package cmd

import (
	"github.com/jjtimmons/defrag/internal/exec"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

// makeCmd represents the make command
var makeCmd = &cobra.Command{
	Use:   "make",
	Short: "Make a target vector from a database of building fragments",
	Long: `Make a vector from its target sequence and a database of building fragments
	
"defrag make" is for assembling a vector, using Gibson Assembly, in the most efficient 
way possible. It's a declaractive approach to vector design. This means that, rather than 
telling make which DNA fragments go together, and in what order, make simply builds the 
vector it's told to. It does this by:

1. Screening fragments in a BLAST database, at "db-fasta," to find those that match
   portions of the target vector sequence
2. Creating a list of possible assemblies using fragments from 1 and ranking them by their
   estimated assembly cost
3. Finding the cheapest assembly-plan with the best Gibson Assembly fragments possible`,
	Run: exec.Execute,
}

// set flags
func init() {
	rootCmd.AddCommand(makeCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	makeCmd.PersistentFlags().StringP("target", "t", "", "Input file name of target vector sequence <FASTA>")
	makeCmd.PersistentFlags().StringP("dbs", "d", "", "Comma separated list of BLAST databases (multi-fasta)")
	makeCmd.PersistentFlags().StringP("out", "o", "", "Output file name")
	makeCmd.PersistentFlags().BoolP("addgene", "a", false, "Use the Addgene repository as a source of building fragments")

	// db path is needed globally
	viper.BindPFlag("dbs", makeCmd.PersistentFlags().Lookup("dbs"))
}

func makeExec(cmd *cobra.Command, args []string) {

	// os.Exit(0)
}
