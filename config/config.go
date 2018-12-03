// Package config is for app wide settings that are unmarshalled
// from Viper (see: /cmd)
package config

import (
	"fmt"
	"log"
	"path"
	"path/filepath"
	"runtime"
	"sort"
	"strings"

	"github.com/spf13/viper"
)

var (
	// logged holds whether we've already logged the settings file being used
	logged = false

	// singleton is a single settings object
	// used so we only read/unmartial the settings file once
	singleton Config
)

// FragmentConfig settings about fragments
type FragmentConfig struct {
	// the minimum homology between this fragment and the net one
	MinHomology int `mapstructure:"min-homology"`

	// the maximum number of fragments in the final assembly
	MaxCount int `mapstructure:"max-count"`
}

// PCRConfig is settings for PCR
type PCRConfig struct {
	// the cost per bp of primer DNA
	BPCost float32 `mapstructure:"bp-cost"`

	// the maximum primer3 score allowable
	P3MaxPenalty float32 `mapstructure:"primer3-penalty-max"`

	// MinLength is the minimum size of a fragment (used to filter BLAST results)
	MinLength int `mapstructure:"min-length"`
}

// SynthCost is meta about the cost of synthesizing DNA up to a certain
// size. Can be fixed (ie everything beneath that limit is the same amount)
// or not (pay by the bp)
type SynthCost struct {
	// whether it's a fixed or variable cost
	Fixed bool `mapstructure:"fixed"`

	// the cost (either per bp or for the whole stretch)
	Dollars float32 `mapstructure:"dollars"`
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

// Config is the root-level settings struct and is a mix
// of settings available in settings.yaml and those
// available from the command line
type Config struct {
	// path to the root of the repo (hackish)
	Root string

	// paths to the fragment DBs
	DBs string

	// whether the user wants to use Addgene as a fragment source
	AddGene bool

	// Fragment level settings
	Fragments FragmentConfig

	// PCR settings
	PCR PCRConfig

	// Synthesis settings
	Synthesis SynthesisConfig

	// Path to a dir for BLAST io files
	BlastDir string

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
		if !logged {
			fmt.Println("Using config file:", viper.ConfigFileUsed())
			logged = true
		}
	} else {
		log.Fatalf("Failed to read in config file: %v", err)
	}

	// move into the singleton Config struct
	if err := viper.Unmarshal(&singleton); err != nil {
		log.Fatalf("Failed to decode settings file %s: %v", viper.ConfigFileUsed(), err)
	}

	// add on a blast dir path for storing BLAST io files
	singleton.BlastDir = filepath.Join(singleton.Root, "bin", "blast")

	// marked as filled to avoid repeating
	singleton.filled = true
	return singleton
}

// init and set viper's paths to the local config file
func init() {
	// path to the root of the app
	_, filename, _, ok := runtime.Caller(0)
	if !ok {
		log.Panicln("No caller information")
	}
	root, _ := filepath.Abs(path.Join(path.Dir(filename), ".."))
	viper.SetDefault("Root", root)

	// addgene's database

	if confFlag := viper.GetString("config"); confFlag != "" {
		viper.AddConfigPath(confFlag) // settings are in root of repo
	} else {
		viper.AddConfigPath(root) // settings are in root of repo
	}

	viper.SetConfigName("settings") // no yaml needed, just a config file called settings
	viper.AutomaticEnv()            // enviornment variables that match
}

// SynthCost returns the cost of synthesizing a stretch of DNA
// of the length passed. It's a stand alone function because
// it's highly dependent on the synth cost settings map
//
// find the smallest synth length greater than fragLength
// Ex: a synthesis provider may say it's 32 cents up to 500bp and
// 60 cents up to 2000bp. So, for a 750bp sequence, we want to use
// the 2000bp price
func (c Config) SynthCost(fragLength int) float32 {
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
	return float32(fragLength) * synthCost.Dollars
}

// DBList returns a list of absolute paths to BLAST databases used during a given run
func (c Config) DBList() (paths []string, err error) {
	if c.AddGene {
		addgenePath := path.Join(c.Root, "assets", "addgene", "db", "addgene")
		return parseDBs(c.DBs + "," + addgenePath)
	}
	return parseDBs(c.DBs)
}

// parseDBs turns a single string of comma separated BLAST dbs into a
// slice of absolute paths to the BLAST dbs on the local fs
func parseDBs(dbList string) (paths []string, err error) {
	noSpaceDBs := strings.Replace(dbList, " ", "", -1)
	for _, db := range strings.Split(noSpaceDBs, ",") {
		absPath, err := filepath.Abs(db)

		if err != nil {
			return nil, fmt.Errorf("failed to create absolute path: %v", err)
		}

		paths = append(paths, absPath)
	}

	return
}
