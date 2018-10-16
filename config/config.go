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
	// the maximum number of fragments in the final assembly
	MaxCount int `mapstructure:"max-count"`

	// the minimum homology between fragments in the Gibson assembly
	MinHomology int `mapstructure:"min-homology"`

	// the minimum length of match between a building fragment
	// and the target fragment for it to be considered
	MinMatch int `mapstructure:"min-match"`
}

// PCRConfig is settings for PCR
type PCRConfig struct {
	// the cost per bp of primer DNA
	BPCost float32 `mapstructure:"bp-cost"`

	// the maximum primer3 score allowable
	P3MaxPenalty float32 `mapstructure:"primer3-penalty-max"`
}

// SynthesisConfig is for settings involving synthesis
type SynthesisConfig struct {
	// the cost per bp of synthesized DNA
	BPCost float32 `mapstructure:"bp-cost"`

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
	_, filename, _, ok := runtime.Caller(0)
	if !ok {
		log.Panicln("No caller information")
	}
	// path to the root of the app
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

	viper.SetConfigName("settings") // no yaml needed
	viper.AutomaticEnv()            // enviornment variables that match
}
