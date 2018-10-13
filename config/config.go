// Package config is for app wide settings that are unmarshalled
// from Viper (see: /cmd)
package config

import (
	"fmt"
	"log"
	"path"
	"path/filepath"

	"github.com/spf13/viper"
)

// MakeFlags are those that are passed to the make command
type MakeFlags struct {
	// the path the local target sequence
	TargetPath string `mapstructure:"target"`

	// whether to use the addgene database as a source of building fragments
	Addgene bool `mapstructure:"addgene"`
}

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
}

// Config is the root-level settings struct and is a mix
// of settings available in settings.yaml and those
// available from the command line
type Config struct {
	// path to the fragment DB
	DB string
	// make settings passed thru CLI
	Make MakeFlags
	// Fragment level settings
	Fragments FragmentConfig
	// PCR settings
	PCR PCRConfig
	// Synthesis settings
	Synthesis SynthesisConfig
}

// New returns a new Config struct populated by settings from
// the adjacent settings.yaml
func New() (c Config) {
	viper.AddConfigPath(".")
	viper.SetConfigFile("settings") // no yaml needed
	viper.AutomaticEnv()            // enviornment variables that match

	// read it intialization files
	if err := viper.ReadInConfig(); err == nil {
		fmt.Println("Using config file:", viper.ConfigFileUsed())
	} else {
		log.Fatalf("Failed to read in config file %s: %v", viper.ConfigFileUsed(), err)
	}

	// move into the new Config struct
	err := viper.Unmarshal(&c)
	if err != nil {
		log.Fatalf("Failed to decode settings file %s: %v", viper.ConfigFileUsed(), err)
	}

	// make path to test db
	if c.DB == "" {
		db, _ := filepath.Abs(path.Join("..", "assets", "addgene", "db", "addgene"))
		c.DB = db
	}

	return
}
