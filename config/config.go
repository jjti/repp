// Package config is for app wide settings that are unmarshalled
// from Viper (see: /cmd)
package config

import (
	"fmt"
	"log"
	"path"
	"path/filepath"
	"runtime"

	"github.com/spf13/viper"
)

var (
	// logged holds whether we've already logged the settings file being used
	logged = false
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
// or not-fixed (pay by the bp)
type SynthCost struct {
	// whether it's a fixed or variable cost
	Fixed bool `mapstructure:"fixed"`

	// the cost (either per bp or for the whole stretch)
	Dollars float32 `mapstructure:"dollars"`
}

// SynthesisConfig is for settings involving synthesis
type SynthesisConfig struct {
	// the cost per bp of synthesized DNA (as a step function)
	Cost map[float32]SynthCost `mapstructure:"cost"`

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
	// path to the fragment DB
	DB string
	// Fragment level settings
	Fragments FragmentConfig
	// PCR settings
	PCR PCRConfig
	// Synthesis settings
	Synthesis SynthesisConfig
}

// New returns a new Config struct populated by settings from
// settings.yaml, in the repo, or some other settings file the user
// points to with the "--config" command
//
// TODO: check for and error out on nonsense config values
// TODO: add back the config file path setting
func New() (c Config) {
	// read in intialization files
	if err := viper.ReadInConfig(); err == nil {
		if !logged {
			fmt.Println("Using config file:", viper.ConfigFileUsed())
			logged = true
		}
	} else {
		log.Fatalf("Failed to read in config file: %v", err)
	}

	// move into the new Config struct
	if err := viper.Unmarshal(&c); err != nil {
		log.Fatalf("Failed to decode settings file %s: %v", viper.ConfigFileUsed(), err)
	}

	return
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
	addgeneDB := path.Join(root, "assets", "addgene", "db", "addgene")
	viper.SetDefault("DB", addgeneDB)

	if confFlag := viper.GetString("config"); confFlag != "" {
		viper.AddConfigPath(confFlag) // settings are in root of repo
	} else {
		viper.AddConfigPath(root) // settings are in root of repo
	}

	viper.SetConfigName("settings") // no yaml needed, just a config file called settings
	viper.AutomaticEnv()            // enviornment variables that match
}
