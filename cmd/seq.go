package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

// SeqCmd represents the make command
var SeqCmd = &cobra.Command{
	Use:   "seq",
	Short: "Build a target vector (seq) from local and/or remote fragments",
	Long: `Build a vector from its target sequence and database(s) of building fragments
	
The fragments are prepared for Gibson Assembly in the most efficient 
way possible

It's a declaractive approach to vector design. This means that, rather than 
telling make which DNA fragments go together, and in what order, "defrag seq" simply builds the 
vector it's told to. It does this by:

1. Screening fragments in BLAST databases, "dbs," to find those that match
   portions of the target vector sequence
2. Creating a list of possible assemblies using fragments from 1 and ranking them by their
   estimated assembly cost
3. Finding the cheapest assembly-plan with the best Gibson Assembly fragments possible`,
	Run: defrag.Execute,
}

// set flags
func init() {
	rootCmd.AddCommand(SeqCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	SeqCmd.Flags().StringP("in", "i", "", "Input file name of target vector sequence <FASTA>")
	SeqCmd.Flags().StringP("out", "o", "", "Output file name for Gibson assemblies")
	SeqCmd.Flags().StringP("dbs", "d", "", "Comma/space separated list of BLAST databases (multi-fasta)")
	SeqCmd.Flags().BoolP("addgene", "a", false, "Use the Addgene repository as a source of building fragments")
}
