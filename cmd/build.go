package cmd

import (
	"github.com/jjtimmons/defrag/config"
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var (
	backboneHelp = `backbone to insert the fragments into. Can either be an entry 
in one of the dbs or a file on the local filesystem.`

	enzymeHelp = `enzyme to linearize the backbone with (backbone must be specified).
'defrag enzymes' prints a list of recognized enzymes.`
)

// buildCmd is for finding building a vector from its fragments, features, or sequence
var buildCmd = &cobra.Command{
	Use:                        "build",
	Short:                      "Build a vector from its fragments, features or sequence",
	SuggestionsMinimumDistance: 2,
	Long: `
Find fragments for assembling a vector via Gibson Assembly. Build the vector
against a list of consituent fragment, feature, or a target sequence.`,
	Aliases: []string{"assemble"},
}

// fragmentsCmd is for piecing together a list of input fragments into a vector
var fragmentsCmd = &cobra.Command{
	Use:                        "fragments",
	Short:                      "Prepare fragments for Gibson assembly",
	Run:                        defrag.FragmentsCmd,
	SuggestionsMinimumDistance: 3,
	Long: `
Prepare a list of fragments for assembly via Gibson Assembly. Fragments are
checked for existing homology with their neighbors and are prepared for
assembly with PCR.`,
}

// featuresCmd is for building a vector from its list of contained features
var featuresCmd = &cobra.Command{
	Use:                        "features [feature] ... [featureN]",
	Short:                      "Find or build a vector from its constituent features",
	Run:                        defrag.FeaturesCmd, // TODO
	SuggestionsMinimumDistance: 3,
	Long:                       "",
}

// sequenceCmd is for assembling a vector (single circular sequence) from its target sequence
var sequenceCmd = &cobra.Command{
	Use:                        "sequence",
	Short:                      "Find or build a vector from its target sequence",
	Run:                        defrag.SequenceCmd,
	SuggestionsMinimumDistance: 3,
	Long: `
Build up a vector from its target sequence using a combination of existing and
synthesized fragments.

Solutions have either a minimum fragment count or assembly cost (or both).`,
}

// set flags
func init() {
	// Flags for specifying the paths to the input file, input fragment files, and output file
	fragmentsCmd.Flags().StringP("in", "i", "", "Input file with a list of building fragments <FASTA>")
	fragmentsCmd.Flags().StringP("out", "o", "", "Output file name with fragments and primers <JSON>")
	fragmentsCmd.Flags().StringP("dbs", "d", "", "Comma separated list of databases with building fragments")
	fragmentsCmd.Flags().BoolP("addgene", "a", false, "Whether to use the Addgene repository as a source of building fragments")
	fragmentsCmd.Flags().BoolP("igem", "g", false, "Whether to use the iGEM repository as a source of building fragments")
	fragmentsCmd.Flags().StringP("backbone", "b", "", backboneHelp)
	fragmentsCmd.Flags().StringP("enzyme", "e", "", enzymeHelp)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	featuresCmd.Flags().StringP("out", "o", "", "output file name")
	featuresCmd.Flags().StringP("dbs", "d", "", "list of local fragment databases")
	featuresCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	featuresCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	featuresCmd.Flags().StringP("backbone", "b", "", backboneHelp)
	featuresCmd.Flags().StringP("enzyme", "e", "", enzymeHelp)
	featuresCmd.Flags().StringP("exclude", "x", "", "keywords for excluding fragments")
	featuresCmd.Flags().IntP("identity", "t", 100, "match %-identity threshold (see 'blastn -help')")

	// Flags for specifying the paths to the input file, input fragment files, and output file
	sequenceCmd.Flags().StringP("in", "i", "", "input FASTA with target sequence")
	sequenceCmd.Flags().StringP("out", "o", "", "output file name")
	sequenceCmd.Flags().StringP("dbs", "d", "", "list of local fragment databases")
	sequenceCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	sequenceCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	sequenceCmd.Flags().StringP("backbone", "b", "", backboneHelp)
	sequenceCmd.Flags().StringP("enzyme", "e", "", enzymeHelp)
	sequenceCmd.Flags().StringP("exclude", "x", "", "keywords for excluding fragments")
	sequenceCmd.Flags().IntP("identity", "t", 100, "match %-identity threshold (see 'blastn -help')")

	buildCmd.AddCommand(fragmentsCmd)
	buildCmd.AddCommand(featuresCmd)
	buildCmd.AddCommand(sequenceCmd)

	// settings is an optional parameter for a settings file (that overrides the fields in BaseSettingsFile)
	buildCmd.PersistentFlags().StringP("settings", "s", config.RootSettingsFile, "build settings")
	viper.BindPFlag("settings", buildCmd.PersistentFlags().Lookup("settings"))

	rootCmd.AddCommand(buildCmd)
}
