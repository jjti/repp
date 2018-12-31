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
	Long: `Build up a vector, by its sequencing, from the best combination of existing
	and synthesized fragments.
	
	All combinations of fragments matching the target sequence are checked
	to find the one with the fewest fragments and lowest overall assembly cost.`,
	Run: defrag.Vector,
}

// set flags
func init() {
	rootCmd.AddCommand(VectorCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	VectorCmd.Flags().StringP("in", "i", "", "Input file with target vector sequence <FASTA>")
	VectorCmd.Flags().StringP("out", "o", "", "Output file name with solutions <JSON>")
	VectorCmd.Flags().StringP("dbs", "d", "", "Comma separated list of databases with building fragments")
	VectorCmd.Flags().BoolP("addgene", "a", false, "Whether to use the Addgene repository as a source of building fragments")
	FragmentsCmd.Flags().StringP("backbone", "b", "",
		`Entry name of a backbone to insert the fragments into. Must match an entry 
	in one of the dbs (either those passed manually or AddGene, iGEM). If an
	enzyme chosen, the backbone will be linearized with that enzyme and
	the largest resulting fragment will be used as the backbone. If no enzyme
	is specified, defrag will chose one nearest the first bp of the backbone with a
	single recognition site`)
	FragmentsCmd.Flags().StringP("enzyme", "e", "",
		`Name of a enzyme to linearize the backbone with (backbone must be specified).
	The enzyme's name must be recognizable by defrag. Use 'defrag enzymes' for a list
	of recognized enzyme names and their recognition sequences`)
}
