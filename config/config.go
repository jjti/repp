// Package config is for app wide settings
package config

import (
	"log"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"sort"

	"github.com/mitchellh/go-homedir"
	"github.com/mitchellh/mapstructure"
	"github.com/spf13/viper"
	"gopkg.in/yaml.v2"
)

var (
	home, _ = homedir.Dir()

	// rvecDir is the root directory where rvec settings and database files live
	rvecDir = filepath.Join(home, ".rvec")

	// RootSettingsFile is the default settings file path for the config package
	RootSettingsFile = filepath.Join(rvecDir, "config.yaml")

	// Primer3Config is the path to the embedded primer3 config directory
	Primer3Config = filepath.Join(rvecDir, "primer3_config") + string(os.PathSeparator)

	// IGEMDB is the path to the iGEM db
	IGEMDB = filepath.Join(rvecDir, "igem")

	// AddgeneDB is the path to the Addgene db
	AddgeneDB = filepath.Join(rvecDir, "addgene")

	// DNASUDB is the path to the DNASU db
	DNASUDB = filepath.Join(rvecDir, "dnasu")

	// FeatureDB is the path to the features db
	FeatureDB = filepath.Join(rvecDir, "features.tsv")

	// EnzymeDB is the path to the enzymes db file
	EnzymeDB = filepath.Join(rvecDir, "enzymes.tsv")
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
	// Vebose is whether to log debug messages to the stdout
	Verbose bool

	// the cost of a single Addgene plasmid
	CostAddgene float64 `mapstructure:"addgene-cost"`

	// the cost of a single part from the iGEM registry
	CostIGEM float64 `mapstructure:"igem-cost"`

	// the per plasmid cost of DNASU plasmids
	CostDNASU float64 `mapstructure:"dnasu-cost"`

	// the cost per bp of primer DNA
	CostBP float64 `mapstructure:"pcr-bp-cost"`

	// the cost of each PCR reaction
	CostPCR float64 `mapstructure:"pcr-rxn-cost"`

	// the cost of time for each PCR reaction
	CostTimePCR float64 `mapstructure:"pcr-time-cost"`

	// the cost of each Gibson Assembly
	CostGibson float64 `mapstructure:"gibson-assembly-cost"`

	// the cost of time for each Gibson Assembly
	CostTimeGibson float64 `mapstructure:"gibson-assembly-time-cost"`

	// the cost per bp of synthesized DNA as a fragment (as a step function)
	CostSyntheticFragment map[int]SynthCost `mapstructure:"synthetic-fragment-cost"`

	// the cost per bp of synthesized clonal DNA  (delivered in a plasmid)
	CostSynthPlasmid map[int]SynthCost `mapstructure:"synthetic-plasmid-cost"`

	// the maximum number of fragments in the final assembly
	FragmentsMaxCount int `mapstructure:"fragments-max-count"`

	// the minimum homology between this fragment and the net one
	FragmentsMinHomology int `mapstructure:"fragments-min-junction-length"`

	// maximum length of homology between two adjacent fragments in bp
	FragmentsMaxHomology int `mapstructure:"fragments-max-junction-length"`

	// maximum allowable hairpin melting temperature (celcius)
	FragmentsMaxHairpinMelt float64 `mapstructure:"fragments-max-junction-hairpin"`

	// PCRMinLength is the minimum size of a fragment (used to filter BLAST results)
	PCRMinLength int `mapstructure:"pcr-min-length"`

	// the maximum primer3 score allowable
	PCRMaxPenalty float64 `mapstructure:"pcr-primer-max-pair-penalty"`

	// the maximum length of a sequence to embed up or downstream of an amplified sequence
	PCRMaxEmbedLength int `mapstructure:"pcr-primer-max-embed-length"`

	// PCRMaxOfftargetTm is the maximum tm of an offtarget, above which PCR is abandoned
	PCRMaxOfftargetTm float64 `mapstructure:"pcr-primer-max-ectopic-tm"`

	// PCRBufferLength is the length of buffer from the ends of a match in which
	// to allow Primer3 to look for a primer
	PCRBufferLength int `mapstructure:"pcr-buffer-length"`

	// maximum length of a synthesized piece of DNA
	SyntheticMaxLength int `mapstructure:"synthetic-max-length"`

	// minimum length of a synthesized piece of DNA
	SyntheticMinLength int `mapstructure:"synthetic-min-length"`
}

// New returns a new Config struct populated by settings from
// config.yaml, in the repo, or some other settings file the user
// points to with the "--config" command
//
// TODO: check for and error out on nonsense config values
// TODO: add back the config file path setting
func New() *Config {
	// read in the default/base settings file first
	viper.SetConfigType("yaml")
	viper.SetConfigFile(RootSettingsFile)
	if err := viper.ReadInConfig(); err != nil {
		log.Fatal(err)
	}

	if userSettings := viper.GetString("settings"); userSettings != "" && userSettings != RootSettingsFile {
		viper.SetConfigFile(userSettings)             // user has specified a new path for a settings file
		if err := viper.MergeInConfig(); err != nil { // read in user defined settings file
			log.Fatal(err)
		}

		file, _ := os.Open(userSettings)
		userData := make(map[string]interface{})
		if err := yaml.NewDecoder(file).Decode(userData); err != nil {
			log.Fatal(err)
		}

		userConfig := &Config{}
		if err := mapstructure.Decode(userData, userConfig); err != nil {
			log.Fatal(err)
		}

		if userConfig.CostSyntheticFragment != nil {
			viper.Set("synthetic-fragment-cost", userConfig.CostSyntheticFragment)
		}
		if userConfig.CostSynthPlasmid != nil {
			viper.Set("synthetic-plasmid-cost", userConfig.CostSynthPlasmid)
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

	if _, err := exec.LookPath("ntthal"); err != nil {
		log.Fatal("no ntthal executable available in PATH, try `make install`")
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
	fragCount := math.Ceil(float64(fragLength) / float64(c.SyntheticMaxLength))
	fragLength = int(math.Floor(float64(fragLength) / float64(fragCount)))

	cost := synthCost(fragLength, c.CostSyntheticFragment)
	if cost.Fixed {
		return fragCount * cost.Cost
	}

	return fragCount * float64(fragLength) * cost.Cost
}

// SynthPlasmidCost returns the cost of synthesizing the insert and having it delivered in a plasmid
func (c Config) SynthPlasmidCost(insertLength int) float64 {
	cost := synthCost(insertLength, c.CostSynthPlasmid)
	if cost.Fixed {
		return cost.Cost
	}

	return float64(insertLength) * cost.Cost
}

// synthCost returns the cost of synthesizing a piece of DNA
func synthCost(seqLength int, costs map[int]SynthCost) SynthCost {
	// find the smallest synth length greater than fragLength
	// Ex: a synthesis provider may say it's 32 cents up to 500bp and
	// 60 cents up to 2000bp. So, for a 750bp sequence, we want to use
	// the 2000bp price
	// TODO: add error here for if there's no cost for seqLength (too large)
	costLengthKeys := make([]int, len(costs))
	for key := range costs {
		costLengthKeys = append(costLengthKeys, key)
	}
	sort.Ints(costLengthKeys)

	synthCostKey := 0
	for _, keyLength := range costLengthKeys {
		if keyLength >= seqLength {
			synthCostKey = keyLength
			break
		}
	}

	// we're not able to make a fragment/gene this large
	// return an extremely large number
	if synthCostKey == 0 {
		return SynthCost{
			Fixed: true,
			Cost:  math.MaxInt32,
		}
	}

	return costs[synthCostKey]
}
