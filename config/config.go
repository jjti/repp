// Package config is for app wide settings that are unmarshalled
// from Viper (see: /cmd)
package config

import (
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

// NewConfig returns a new Config struct populated by
// Viper settings (either from the local settings.yaml)
// and/or command line arguments
func NewConfig() Config {
	var c Config

	err := viper.Unmarshal(&c)
	if err != nil {
		log.Fatalf("unable to decode into struct, %v", err)
	}

	// make path to test db
	db, _ := filepath.Abs(path.Join("..", "assets", "addgene", "db", "addgene"))
	c.DB = db

	return c
}
