// Package config is for app wide settings that are unmarshalled
// from Viper (see: /cmd)
package config

import (
	"log"
	"math"
	"os/exec"
	"path/filepath"
	"sort"

	"github.com/mitchellh/go-homedir"
	"github.com/spf13/viper"
)

var (
	home, _ = homedir.Dir()

	// BaseDir is the root directory where defrag settings and database files live
	BaseDir = filepath.Join(home, ".defrag")

	// BaseSettingsFile is the default settings file path for the config package
	BaseSettingsFile = filepath.Join(BaseDir, "config.yaml")

	// IGEMDB is the path to the iGEM db
	IGEMDB = filepath.Join(BaseDir, "igem")

	// AddgeneDB is the path to the Addgene db
	AddgeneDB = filepath.Join(BaseDir, "addgene")

	// FeatureDB is the path to the features db
	FeatureDB = filepath.Join(BaseDir, "features.tsv")
)

// SynthCost contains data of the cost of synthesizing DNA up to a certain
// size. Can be fixed (ie everything beneath that limit is the same amount)
// or not (pay by the bp)
type SynthCost struct {
	// whether it's a fixed or variable cost
	Fixed bool `mapstructure:"fixed"`

	// the cost (either per bp or for the whole stretch)
	Cost float64 `mapstructure:"cost"`
}

// Config is the Root-level settings struct and is a mix
// of settings available in config.yaml and those
// available from the command line
type Config struct {
	// the cost of a single Addgene vector
	AddGeneVectorCost float64 `mapstructure:"addgene-vector-cost"`

	// the cost of a single part from the iGEM registry
	IGEMPartCost float64 `mapstructure:"igem-part-cost"`

	// primer3 config folder, needed for thermodynamic calculations
	// created in the settings file during `make install`
	Primer3Config string `mapstructure:"primer3_config-path"`

	// the maximum number of fragments in the final assembly
	FragmentsMaxCount int `mapstructure:"fragments-max-count"`

	// the minimum homology between this fragment and the net one
	FragmentsMinHomology int `mapstructure:"fragments-min-homology"`

	// maximum length of homology between two adjacent fragments in bp
	FragmentsMaxHomology int `mapstructure:"fragments-max-homology"`

	// the cost per bp of primer DNA
	PCRBPCost float64 `mapstructure:"pcr-bp-cost"`

	// PCRMinLength is the minimum size of a fragment (used to filter BLAST results)
	PCRMinLength int `mapstructure:"pcr-min-length"`

	// the maximum primer3 score allowable
	PCRP3MaxPenalty float64 `mapstructure:"pcr-primer3-penalty-max"`

	// the maximum length of a sequence to embed up or downstream of an amplified sequence
	PCRMaxEmbedLength int `mapstructure:"pcr-primer-max-embed-length"`

	// PCRMaxOfftargetTm is the maximum tm of an offtarget, above which PCR is abandoned
	PCRMaxOfftargetTm float64 `mapstructure:"pcr-primer-max-offtarget-tm"`

	// the cost per bp of synthesized DNA as a fragment (as a step function)
	SynthesisFragmentCost map[int]SynthCost `mapstructure:"synthesis-fragment-cost"`

	// the cost per bp of synthesized clonal DNA  (delivered in a vector)
	SynthesisVectorCost map[int]SynthCost `mapstructure:"synthesis-gene-cost"`

	// maximum length of a synthesized piece of DNA
	SynthesisMaxLength int `mapstructure:"synthesis-max-length"`

	// minimum length of a synthesized piece of DNA
	SynthesisMinLength int `mapstructure:"synthesis-min-length"`
}

// New returns a new Config struct populated by settings from
// config.yaml, in the repo, or some other settings file the user
// points to with the "--config" command
//
// TODO: check for and error out on nonsense config values
// TODO: add back the config file path setting
func New() *Config {
	// read in the default/base settings file first
	viper.SetConfigFile(BaseSettingsFile)
	if err := viper.ReadInConfig(); err != nil {
		log.Fatal(err)
	}

	if userSettings := viper.GetString("settings"); userSettings != "" && userSettings != BaseSettingsFile {
		viper.SetConfigFile(userSettings) // user has specified a new path for a settings file

		// read in user defined settings file
		if err := viper.MergeInConfig(); err == nil {
			// fmt.Println("Using config file: ", viper.ConfigFileUsed())
		} else {
			log.Fatal(err)
		}
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

	// build Config
	config := &Config{}
	if err := viper.Unmarshal(&config); err != nil {
		log.Fatalf("failed to decode settings file %s: %v", viper.ConfigFileUsed(), err)
	}

	return config
}

// SynthFragmentCost returns the cost of synthesizing a linear stretch of DNA
func (c Config) SynthFragmentCost(fragLength int) float64 {
	// by default, we try to synthesize the whole thing in one piece
	// we may optionally need to split it into multiple
	fragCount := math.Ceil(float64(fragLength) / float64(c.SynthesisMaxLength))
	fragLength = int(math.Floor(float64(fragLength) / float64(fragCount)))

	cost := synthCost(fragLength, c.SynthesisFragmentCost)
	if cost.Fixed {
		return fragCount * cost.Cost
	}

	return fragCount * float64(fragLength) * cost.Cost
}

// SynthGeneCost returns the cost of synthesizing the insert and having it delivered in a vector
func (c Config) SynthGeneCost(insertLength int) float64 {
	cost := synthCost(insertLength, c.SynthesisVectorCost)
	if cost.Fixed {
		return cost.Cost
	}

	return float64(insertLength) * cost.Cost
}

// synthCost returns the cost of synthesizing the piece of DNA,
// either linear or circular
//
// find the smallest synth length greater than fragLength
// Ex: a synthesis provider may say it's 32 cents up to 500bp and
// 60 cents up to 2000bp. So, for a 750bp sequence, we want to use
// the 2000bp price
//
// TODO: add error here for if there's no cost for seqLength (too large)
func synthCost(seqLength int, costs map[int]SynthCost) SynthCost {
	costLengthKeys := make([]int, len(costs))
	for key := range costs {
		costLengthKeys = append(costLengthKeys, key)
	}
	sort.Ints(costLengthKeys)

	synthCostKey := 0
	for _, keyLength := range costLengthKeys {
		if keyLength > seqLength {
			synthCostKey = keyLength
			break
		}
	}

	return costs[synthCostKey]
}
