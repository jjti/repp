package cmd

import (
	"log"

	"github.com/spf13/viper"

	"github.com/jjtimmons/decvec/config"
	"github.com/jjtimmons/decvec/internal/dvec"
	"github.com/spf13/cobra"
)

// makeCmd represents the make command
var makeCmd = &cobra.Command{
	Use:   "make",
	Short: "Make a target vector from a database of building fragments",
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
	makeCmd.PersistentFlags().StringP("target", "t", "", "Path to FASTA file with target vector sequence")
	makeCmd.PersistentFlags().StringP("db", "d", "", "Database of building fragments")

	viper.BindPFlag("db", makeCmd.PersistentFlags().Lookup("db"))

}

// makeExec is the root of the make functionality
//
// the goal is to find an "optimal" assembly vector with:
// 	1. fewest number of fragments
// 	2. cheapest assembly cost ($)
// and, secondarily:
//	3. no duplicate end regions between Gibson fragments
// 	4. no inverted repeats in the junctions
// 	5. no off-target binding sites in the parent vectors
//	6. low primer3 penalty scores
func makeExec(cmd *cobra.Command, args []string) {
	target, err := cmd.Flags().GetString("target")
	if err != nil {
		log.Fatalf("Cannot get target from arguments")
	}

	c := config.New()

	// no path to input file
	if target == "" {
		log.Fatalf("Failed, no target fragment path set in config %+v", c)
	}

	// read in fragments
	fragments, err := dvec.Read(target)
	if err != nil {
		log.Fatalf("Failed to read in fasta files at %s: %v", target, err)
	}

	// set target fragment
	if len(fragments) > 1 {
		println(
			"Warning: %d building fragments were in %s. Only targeting the first: %s",
			len(fragments),
			target,
			fragments[0].ID,
		)
	}

	// target := fragments[0]

	// dbPath, err := filepath.Abs(path.Join("..", "..", "assets", "addgene", "db", "addgene"))
	// handle(err)
}
