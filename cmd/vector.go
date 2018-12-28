package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

// VectorCmd is for assembling a vector (single circular sequence) after finding
// building fragments to piece it together efficiently
var VectorCmd = &cobra.Command{
	Use:   "vector",
	Short: "Build a target vector from local and/or remote fragments",
	Long: `Finds solutions to build an input sequence from database(s) of building fragments
via Gison Assembly
	
All combinations of fragments, that match the target sequence, are checked
to find the one with the fewest fragments and lowest overall
assembly cost`,
	Run: defrag.Vector,
}

// set flags
func init() {
	rootCmd.AddCommand(VectorCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	VectorCmd.Flags().StringP("in", "i", "", "Input file with target vector sequence <FASTA>")
	VectorCmd.Flags().StringP("out", "o", "", "Output file name with solutions <JSON>")
	VectorCmd.Flags().StringP("dbs", "d", "", "Comma separated list of BLAST databases (multi-FASTA)")
	VectorCmd.Flags().BoolP("addgene", "a", false, "Use the Addgene repository as a source of building fragments")
}
