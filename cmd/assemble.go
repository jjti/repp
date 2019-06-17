package cmd

import (
	"github.com/jjtimmons/plade/config"
	"github.com/jjtimmons/plade/internal/plade"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var (
	backboneHelp = `backbone to insert the fragments into. Can either be an entry 
in one of the dbs or a file on the local filesystem.`

	enzymeHelp = `enzyme to linearize the backbone with (backbone must be specified).
'plade ls enzymes' prints a list of recognized enzymes.`
)

// assembleCmd is for finding building a plasmid from its fragments, features, or sequence
var assembleCmd = &cobra.Command{
	Use:                        "assemble",
	Short:                      "Assemble a plasmid from its fragments, features or sequence",
	SuggestionsMinimumDistance: 3,
	Long: `Find fragments for assembling a plasmid via Gibson Assembly. Build the plasmid
against a list of consituent fragment, feature, or a target sequence.`,
	Aliases: []string{"make", "build"},
}

// fragmentsCmd is for piecing together a list of input fragments into a plasmid
var fragmentsCmd = &cobra.Command{
	Use:                        "fragments",
	Short:                      "Build a plasmid from its constituent fragments",
	Run:                        plade.FragmentsCmd,
	SuggestionsMinimumDistance: 3,
	Long: `Prepare a list of fragments for assembly via Gibson Assembly. Fragments are
checked for existing homology with their neighbors and are prepared for
assembly with PCR.`,
}

// featuresCmd is for building a plasmid from its list of contained features
var featuresCmd = &cobra.Command{
	Use:                        "features [feature] ... [featureN]",
	Short:                      "Find or build a plasmid from its constituent features",
	Run:                        plade.FeaturesCmd, // TODO
	SuggestionsMinimumDistance: 3,
}

// sequenceCmd is for assembling a plasmid (single circular sequence) from its target sequence
var sequenceCmd = &cobra.Command{
	Use:                        "sequence",
	Short:                      "Find or build a plasmid from its target sequence",
	Run:                        plade.SequenceCmd,
	SuggestionsMinimumDistance: 2,
	Long: `Build up a plasmid from its target sequence using a combination of existing and
synthesized fragments.

Solutions have either a minimum fragment count or assembly cost (or both).`,
	Aliases: []string{"seq", "plasmid"},
}

// set flags
func init() {
	// Flags for specifying the paths to the input file, input fragment files, and output file
	fragmentsCmd.Flags().StringP("in", "i", "", "input file name")
	fragmentsCmd.Flags().StringP("out", "o", "", "output file name")
	fragmentsCmd.Flags().StringP("dbs", "d", "", "comma separated list of local fragment databases")
	fragmentsCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	fragmentsCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	fragmentsCmd.Flags().BoolP("dnasu", "u", false, "use the DNASU repository")
	fragmentsCmd.Flags().StringP("backbone", "b", "", backboneHelp)
	fragmentsCmd.Flags().StringP("enzyme", "e", "", enzymeHelp)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	featuresCmd.Flags().StringP("out", "o", "", "output file name")
	featuresCmd.Flags().StringP("dbs", "d", "", "comma separated list of local fragment databases")
	featuresCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	featuresCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	featuresCmd.Flags().BoolP("dnasu", "u", false, "use the DNASU repository")
	featuresCmd.Flags().StringP("backbone", "b", "", backboneHelp)
	featuresCmd.Flags().StringP("enzyme", "e", "", enzymeHelp)
	featuresCmd.Flags().StringP("exclude", "x", "", "keywords for excluding fragments")
	featuresCmd.Flags().IntP("identity", "p", 98, "%-identity threshold (see 'blastn -help')")

	// Flags for specifying the paths to the input file, input fragment files, and output file
	sequenceCmd.Flags().StringP("in", "i", "", "input FASTA with target sequence")
	sequenceCmd.Flags().StringP("out", "o", "", "output file name")
	sequenceCmd.Flags().StringP("dbs", "d", "", "list of local fragment databases")
	sequenceCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	sequenceCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	sequenceCmd.Flags().BoolP("dnasu", "u", false, "use the DNASU repository")
	sequenceCmd.Flags().StringP("backbone", "b", "", backboneHelp)
	sequenceCmd.Flags().StringP("enzyme", "e", "", enzymeHelp)
	sequenceCmd.Flags().StringP("exclude", "x", "", "keywords for excluding fragments")
	sequenceCmd.Flags().IntP("identity", "p", 98, "%-identity threshold (see 'blastn -help')")

	assembleCmd.AddCommand(fragmentsCmd)
	assembleCmd.AddCommand(featuresCmd)
	assembleCmd.AddCommand(sequenceCmd)

	// settings is an optional parameter for a settings file (that overrides the fields in BaseSettingsFile)
	assembleCmd.PersistentFlags().StringP("settings", "s", config.RootSettingsFile, "build settings")
	assembleCmd.PersistentFlags().BoolP("verbose", "v", false, "whether to log results to stdout")
	viper.BindPFlag("settings", assembleCmd.PersistentFlags().Lookup("settings"))
	viper.BindPFlag("verbose", assembleCmd.PersistentFlags().Lookup("verbose"))

	rootCmd.AddCommand(assembleCmd)
}
