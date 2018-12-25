// Package config is for app wide settings that are unmarshalled
// from Viper (see: /cmd)
package config

import (
	"fmt"
	"log"
	"os"
	"path"
	"path/filepath"
	"runtime"
	"sort"

	"github.com/spf13/viper"
)

var (
	// Root is the Root directory of the app (set in init)
	Root = ""

	// singleton is a single settings object
	// used so we only read/unmartial the settings file once
	singleton Config
)

// FragmentConfig settings about fragments
type FragmentConfig struct {
	// the maximum number of fragments in the final assembly
	MaxCount int `mapstructure:"max-count"`

	// the minimum homology between this fragment and the net one
	MinHomology int `mapstructure:"min-homology"`

	// maximum length of homology between two adjacent fragments in bp
	MaxHomology int `mapstructure:" max-homology"`
}

// PCRConfig is settings for PCR
type PCRConfig struct {
	// the cost per bp of primer DNA
	BPCost float64 `mapstructure:"bp-cost"`

	// MinLength is the minimum size of a fragment (used to filter BLAST results)
	MinLength int `mapstructure:"min-length"`

	// the maximum primer3 score allowable
	P3MaxPenalty float64 `mapstructure:"primer3-penalty-max"`

	// the maximum length of a sequence to embed up or downstream of an amplified sequence
	MaxEmbedLength int `mapstructure:"primer-max-embed-length"`

	// MaxOfftargetTm is the maximum tm of an offtarget, above which PCR is abandoned
	MaxOfftargetTm float64 `mapstructure:"primer-max-offtarget-tm"`
}

// SynthCost is meta about the cost of synthesizing DNA up to a certain
// size. Can be fixed (ie everything beneath that limit is the same amount)
// or not (pay by the bp)
type SynthCost struct {
	// whether it's a fixed or variable cost
	Fixed bool `mapstructure:"fixed"`

	// the cost (either per bp or for the whole stretch)
	Dollars float64 `mapstructure:"dollars"`
}

// SynthesisConfig is for synthesis settings
type SynthesisConfig struct {
	// the cost per bp of synthesized DNA (as a step function)
	Cost map[int]SynthCost `mapstructure:"cost"`

	// maximum length of a synthesized piece of DNA
	MaxLength int `mapstructure:"max-length"`

	// minimum length of a synthesized piece of DNA
	MinLength int `mapstructure:"min-length"`
}

// VendorConfig holds the paths to binaries needed for the library:
// blastn, blastdbcmd, and primer3_core
type VendorConfig struct {
	// blastn binary
	Blastn string

	// blastdbcmd binary
	Blastdbcmd string

	// blast io subdirectory
	Blastdir string

	// primer3_core binary
	Primer3core string

	// primer3 config folder, needed for thermodynamic calculations
	Primer3config string

	// primer3 io subdirectory
	Primer3dir string
}

// Config is the Root-level settings struct and is a mix
// of settings available in settings.yaml and those
// available from the command line
type Config struct {
	// paths to the fragment databases
	DBs string

	// whether the user wants to use Addgene as a fragment source
	AddGene bool

	// the cost of a single Addgene vector
	AddGeneVectorCost float64 `mapstructure:"addgene-vector-cost"`

	// Fragment level settings
	Fragments FragmentConfig

	// PCR settings
	PCR PCRConfig

	// Synthesis settings
	Synthesis SynthesisConfig

	// filled is a flag to mark whether it's actually been made yet
	filled bool
}

// New returns a new Config struct populated by settings from
// settings.yaml, in the repo, or some other settings file the user
// points to with the "--config" command
//
// TODO: check for and error out on nonsense config values
// TODO: add back the config file path setting
func New() Config {
	// check if the singleton has been defined, return as is if so
	if singleton.filled {
		return singleton
	}

	// read in intialization files
	if err := viper.ReadInConfig(); err == nil {
		fmt.Println("Using config file: ", viper.ConfigFileUsed())
	} else {
		log.Fatalf("Failed to read in config file: %v", err)
	}

	// move into the singleton Config struct
	if err := viper.Unmarshal(&singleton); err != nil {
		log.Fatalf("Failed to decode settings file %s: %v", viper.ConfigFileUsed(), err)
	}

	// marked as filled to avoid repeating
	singleton.filled = true
	return singleton
}

// SynthCost returns the cost of synthesizing a stretch of DNA
// of the length passed. It's a stand alone function because
// it's highly dependent on the synth cost settings map
//
// find the smallest synth length greater than fragLength
// Ex: a synthesis provider may say it's 32 cents up to 500bp and
// 60 cents up to 2000bp. So, for a 750bp sequence, we want to use
// the 2000bp price
func (c Config) SynthCost(fragLength int) float64 {
	costLengthKeys := []int{}
	for key := range c.Synthesis.Cost {
		costLengthKeys = append(costLengthKeys, key)
	}
	sort.Ints(costLengthKeys)

	synthCostKey := 0
	for _, costLength := range costLengthKeys {
		if costLength > fragLength {
			synthCostKey = costLength
			break
		}
	}

	// find whether this fragment has a fixed or variable cost
	synthCost := c.Synthesis.Cost[synthCostKey]
	if synthCost.Fixed {
		return synthCost.Dollars
	}
	return float64(fragLength) * synthCost.Dollars
}

// Vendors returns a new config for paths to library dependencies. It also creates
// subdirectories for primer3 and blast to store their input+output files in
// (mostly for debugging)
func (c Config) Vendors() VendorConfig {
	blastn := filepath.Join(Root, "vendor", "ncbi-blast-2.7.1+", "bin", "blastn")
	if _, err := os.Stat(blastn); os.IsNotExist(err) {
		log.Fatalf("failed to find a BLAST executable at %s", blastn)
	}

	blastdir := filepath.Join(Root, "bin", "blast")
	if err := os.MkdirAll(blastdir, os.ModePerm); err != nil {
		log.Fatalf("Failed to create a BLAST dir: %v", err)
	}

	blastdbcmd := filepath.Join(Root, "vendor", "ncbi-blast-2.7.1+", "bin", "blastdbcmd")
	if _, err := os.Stat(blastdbcmd); err != nil {
		log.Fatalf("Failed to locate blastdbcmd executable: %v", err)
	}

	p3core := filepath.Join(Root, "vendor", "primer3-2.4.0", "src", "primer3_core")
	if _, err := os.Stat(p3core); err != nil {
		log.Fatalf("Failed to locate primer3 executable: %v", err)
	}

	p3conf := filepath.Join(Root, "vendor", "primer3-2.4.0", "src", "primer3_config") + string(filepath.Separator)
	if _, err := os.Stat(p3conf); err != nil {
		log.Fatalf("Failed to locate primer3 config folder: %v", err)
	}

	p3dir := filepath.Join(Root, "bin", "primer3")
	if err := os.MkdirAll(p3dir, os.ModePerm); err != nil {
		log.Fatalf("Failed to create a primer3 outut dir: %v", err)
	}

	return VendorConfig{
		Blastn:        blastn,
		Blastdir:      blastdir,
		Blastdbcmd:    blastdbcmd,
		Primer3core:   p3core,
		Primer3config: p3conf,
		Primer3dir:    p3dir,
	}
}

// init and set viper's paths to the local config file
func init() {
	_, filename, _, ok := runtime.Caller(0)
	if !ok {
		log.Panicln("No caller information")
	}
	Root, _ = filepath.Abs(path.Join(path.Dir(filename), ".."))

	if confFlag := viper.GetString("config"); confFlag != "" {
		viper.AddConfigPath(confFlag) // settings are in Root of repo
	} else {
		viper.AddConfigPath(Root) // settings are in Root of repo
	}

	viper.SetConfigName("settings") // no yaml needed, just a config file called settings
	viper.AutomaticEnv()            // enviornment variables that match
}
