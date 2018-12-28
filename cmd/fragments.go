package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

// FragmentsCmd is for piecing together a list of input fragments into a vector
// and preparing the fragments to assemble into that vector
var FragmentsCmd = &cobra.Command{
	Use:   "fragments",
	Short: "Assemble a list of fragments via Gibson Assembly",
	Long: `Accepts a list of linearized framgents and prepares them for assembly
via Gibson Assembly. Fragments are checked for existing homology for their
neighbors (in the multi-FASTA input file). Fragments lacking the requisite
are prepared for assembly with PCR to add the homology.`,
	Run: defrag.Fragments,
}

// set flags
func init() {
	rootCmd.AddCommand(FragmentsCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	FragmentsCmd.Flags().StringP("in", "i", "", "Input file with a list of building fragments <multi-FASTA>")
	FragmentsCmd.Flags().StringP("out", "o", "", "Output file name with fragments and primers <JSON>")
	// TODO: let the user run their assembly through `defrag vector` with the passed building fragments
	// FragmentsCmd.Flags().StringP("dbs", "d", "", "Comma separated list of BLAST databases (multi-FASTA)")
}
