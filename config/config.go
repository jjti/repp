// Package config is for app wide settings that are unmarshalled
// from Viper (see: /cmd)
package config

import (
	"fmt"
	"log"
	"os/exec"
	"path"
	"sort"

	"github.com/spf13/viper"
)

// FragmentConfig settings about fragments
type FragmentConfig struct {
	// the maximum number of fragments in the final assembly
	MaxCount int `mapstructure:"max-count"`

	// the minimum homology between this fragment and the net one
	MinHomology int `mapstructure:"min-homology"`

	// maximum length of homology between two adjacent fragments in bp
	MaxHomology int `mapstructure:"max-homology"`
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

	// the cost of a single part from the iGEM registry
	IGEMPartCost float64 `mapstructure:"igem-part-cost"`

	// primer3 config folder, needed for thermodynamic calculations
	// created in the settings file during `make install`
	Primer3config string `mapstructure:"primer3_config-path"`

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
func New() *Config {
	if userConfigPath := viper.GetString("config"); userConfigPath != "" {
		viper.AddConfigPath(userConfigPath) // user has specified a path to a settings file
	} else {
		viper.AddConfigPath(path.Join("etc", "defrag")) // settings are /etc/defrag
	}
	viper.SetConfigName("settings") // no yaml needed, just a config file called settings

	// read in intialization files
	if err := viper.ReadInConfig(); err == nil {
		fmt.Println("Using config file: ", viper.ConfigFileUsed())
	} else {
		log.Fatalf("Failed to read in config file: %v", err)
	}

	// make sure all depedencies are available (may belong elsewhere)
	if _, err := exec.LookPath("blastn"); err != nil {
		log.Fatal("no blastn executable available in PATH, try `make install`")
	}

	if _, err := exec.LookPath("blastdbcmd"); err != nil {
		log.Fatal("no blastdbcmd executable available in PATH, try `make install`")
	}

	if _, err := exec.LookPath("primer3_core"); err != nil {
		log.Fatal("no primer3_core executable available in PATH, try `make install`")
	}

	// move into the singleton Config struct
	config := &Config{}
	if err := viper.Unmarshal(&config); err != nil {
		log.Fatalf("Failed to decode settings file %s: %v", viper.ConfigFileUsed(), err)
	}
	return config
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
