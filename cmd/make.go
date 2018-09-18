package cmd

import (
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var targetPath string
var inputFastaPath string

// makeCmd represents the make command
var makeCmd = &cobra.Command{
	Use:   "make",
	Short: "Make a vector from its target sequence and a database of building fragments",
	Long: `Make a vector from its target sequence and a database of building fragments
	
"decvec make" is for assembling a vector, using Gibson Assembly, in the most efficient 
way possible. It's a declaractive approach to vector design. This means that, rather than 
telling make which DNA fragments go together, and in what order, make simply builds the 
vector it's told to. It does this by:

1. Screening fragments in a BLAST database, at "db-fasta," to find those that match
   portions of the target vector sequence
2. Creating a list of possible assemblies using fragments from 1 and ranking them by their
   estimated assembly cost
3. Finding the cheapest assembly-plan with the best Gibson Assembly fragments possible`,
	Run: makeExec,
}

// set flags
func init() {
	rootCmd.AddCommand(makeCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	makeCmd.Flags().StringVarP(&targetPath, "target", "t", "", "path to a FASTA file with the target vector")
	makeCmd.Flags().StringVarP(&inputFastaPath, "fragment-db", "f", "", "path to a BLAST database with possible building fragments")

	// Mark required flags
	makeCmd.MarkFlagRequired("target")
	makeCmd.MarkFlagRequired("fragment-db")

	// Bind the paramters to viper
	viper.BindPFlag("make.target", makeCmd.Flags().Lookup("target"))
	viper.BindPFlag("make.fragment-db", makeCmd.Flags().Lookup("fragment-db"))
}

// makeExec is the root of the make functionality
func makeExec(cmd *cobra.Command, args []string) {
	// _ := config.NewConfig()

	println("called make")
}
