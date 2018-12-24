package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

// SeqCmd represents the make command
var SeqCmd = &cobra.Command{
	Use:   "seq",
	Short: "Build a target vector (seq) from local and/or remote fragments",
	Long: `Finds solutions to build an input sequence from database(s) of building fragments
via Gison Assembly
	
All combinations of fragments, that match the target sequence, are checked
to find the one with the least assembly fragments and the lowest overall
assembly cost`,
	Run: defrag.Execute,
}

// set flags
func init() {
	rootCmd.AddCommand(SeqCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	SeqCmd.Flags().StringP("in", "i", "", "Input file name of target vector sequence <FASTA>")
	SeqCmd.Flags().StringP("out", "o", "", "Output file name for solutions")
	SeqCmd.Flags().StringP("dbs", "d", "", "Comma/space separated list of BLAST databases (multi-fasta)")
	SeqCmd.Flags().BoolP("addgene", "a", false, "Use the Addgene repository as a source of building fragments")
}
