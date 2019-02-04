package cmd

import (
	"github.com/jjtimmons/defrag/config"
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
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
Build a vectorfrom  its fragments, features, or sequence.`,
	Aliases: []string{"assemble"},
}

// fragmentsCmd is for piecing together a list of input fragments into a vector
// and preparing the fragments to assemble into that vector
var fragmentsCmd = &cobra.Command{
	Use:                        "fragments",
	Short:                      "Assemble a list of fragments",
	Run:                        defrag.FragmentsCmd,
	SuggestionsMinimumDistance: 4,
	Long: `
Prepare a list of fragments for assembly via Gibson Assembly. Fragments are
checked for existing homology with their neighbors and are prepared for
assembly with PCR.`,
}

// featuresCmd is for building a vector from its list of contained features. Also
// includes sub-commands for updating the feature database
var featuresCmd = &cobra.Command{
	Use:                        "features [feature] ... [featureN]",
	Short:                      "Find or build a vector from its features",
	Run:                        defrag.FeaturesCmd, // TODO
	SuggestionsMinimumDistance: 4,
	Long:                       "",
}

// vectorCmd is for assembling a vector (single circular sequence) from its target sequence
var vectorCmd = &cobra.Command{
	Use:                        "vector",
	Short:                      "Build a vector with the target sequence using local and/or remote fragments",
	Run:                        defrag.VectorCmd,
	SuggestionsMinimumDistance: 4,
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
	featuresCmd.Flags().StringP("dbs", "d", "", "delimited list of local fragment databases")
	featuresCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	featuresCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	featuresCmd.Flags().StringP("backbone", "b", "", backboneHelp)
	featuresCmd.Flags().StringP("enzyme", "e", "", enzymeHelp)
	featuresCmd.Flags().StringP("filter", "f", "", "delimited keywords for removing fragments")
	featuresCmd.Flags().IntP("identity", "t", 100, "match %-identity threshold (see 'blastn -help')")

	// Flags for specifying the paths to the input file, input fragment files, and output file
	vectorCmd.Flags().StringP("in", "i", "", "input FASTA with target sequence")
	vectorCmd.Flags().StringP("out", "o", "", "output file name")
	vectorCmd.Flags().StringP("dbs", "d", "", "delimited list of local fragment databases")
	vectorCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	vectorCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	vectorCmd.Flags().StringP("backbone", "b", "", backboneHelp)
	vectorCmd.Flags().StringP("enzyme", "e", "", enzymeHelp)
	vectorCmd.Flags().StringP("filter", "f", "", "delimited keywords for removing fragments")
	vectorCmd.Flags().IntP("identity", "t", 100, "match %-identity threshold (see 'blastn -help')")

	buildCmd.AddCommand(fragmentsCmd)
	buildCmd.AddCommand(featuresCmd)
	buildCmd.AddCommand(vectorCmd)

	// settings is an optional parameter for a settings file (that overrides the fields in BaseSettingsFile)
	buildCmd.PersistentFlags().StringP("settings", "s", config.RootSettingsFile, "build settings")

	rootCmd.AddCommand(buildCmd)
}
