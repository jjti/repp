package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

// VectorCmd is for assembling a vector (single circular sequence) after finding
// building fragments to piece it together efficiently
var VectorCmd = &cobra.Command{
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
	rootCmd.AddCommand(VectorCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	VectorCmd.Flags().StringP("in", "i", "", "Input file with target vector sequence <FASTA>")
	VectorCmd.Flags().StringP("out", "o", "", "Output file name with solutions <JSON>")
	VectorCmd.Flags().StringP("dbs", "d", "", "Comma separated list of databases with building fragments")
	VectorCmd.Flags().BoolP("addgene", "a", false, "Whether to use the Addgene repository as a source of building fragments")
	VectorCmd.Flags().BoolP("igem", "g", false, "Whether to use the iGEM repository as a source of building fragments")
	VectorCmd.Flags().StringP("backbone", "b", "", backboneHelp)
	VectorCmd.Flags().StringP("enzyme", "e", "", enzymeHelp)
}
