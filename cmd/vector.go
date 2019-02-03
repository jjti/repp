package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

// vectorCmd is for assembling a vector (single circular sequence) after finding
// building fragments to piece it together efficiently
var vectorCmd = &cobra.Command{
	Use:                        "vector",
	Short:                      "Build a target vector from local and/or remote fragments",
	Run:                        defrag.VectorCmd,
	SuggestionsMinimumDistance: 4,
	Long: `Build up a vector, by its sequence, from some combination of existing and
synthesized fragments.

Combinations of fragments matching the target sequence are checked to find the
one with the fewest fragments and lowest overall assembly cost.`,
}

// set flags
func init() {
	rootCmd.AddCommand(vectorCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	vectorCmd.Flags().StringP("in", "i", "", "Input file with target vector sequence <FASTA>")
	vectorCmd.Flags().StringP("out", "o", "", "Output file name with solutions <JSON>")
	vectorCmd.Flags().StringP("dbs", "d", "", "Comma separated list of local building fragment databases (same machine)")
	vectorCmd.Flags().BoolP("addgene", "a", false, "Whether to use the Addgene repository as a source of building fragments")
	vectorCmd.Flags().BoolP("igem", "g", false, "Whether to use the iGEM repository as a source of building fragments")
	vectorCmd.Flags().StringP("backbone", "b", "", backboneHelp)
	vectorCmd.Flags().StringP("enzyme", "e", "", enzymeHelp)
	vectorCmd.Flags().StringP("filter", "f", "", "Comma separated list of keywords to filter against building fragments")
	vectorCmd.Flags().Float64P("identity", "t", 100, "Percentage identity threshold for identifying building fragments (see 'blastn -help')")
}
