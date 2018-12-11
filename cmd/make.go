package cmd

import (
	"github.com/jjtimmons/defrag/internal/exec"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

// MakeCmd represents the make command
var MakeCmd = &cobra.Command{
	Use:   "make",
	Short: "Make a target vector from a database of building fragments",
	Long: `Make a vector from its target sequence and a database of building fragments
	
"defrag make" is for assembling a vector, using Gibson Assembly, in the most efficient 
way possible. It's a declaractive approach to vector design. This means that, rather than 
telling make which DNA fragments go together, and in what order, make simply builds the 
vector it's told to. It does this by:

1. Screening fragments in BLAST databases, "dbs," to find those that match
   portions of the target vector sequence
2. Creating a list of possible assemblies using fragments from 1 and ranking them by their
   estimated assembly cost
3. Finding the cheapest assembly-plan with the best Gibson Assembly fragments possible`,
	Run: exec.Execute,
}

// set flags
func init() {
	rootCmd.AddCommand(MakeCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	MakeCmd.PersistentFlags().StringP("target", "t", "", "Input file name of target vector sequence <FASTA>")
	MakeCmd.PersistentFlags().StringP("dbs", "d", "", "Comma separated list of BLAST databases (multi-fasta)")
	MakeCmd.PersistentFlags().StringP("out", "o", "", "Output file name")
	MakeCmd.PersistentFlags().BoolP("addgene", "a", false, "Use the Addgene repository as a source of building fragments")

	// db path is needed globally
	viper.BindPFlag("dbs", MakeCmd.PersistentFlags().Lookup("dbs"))
}
