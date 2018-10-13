package cmd

import (
	"log"

	"github.com/jjtimmons/decvec/config"
	"github.com/jjtimmons/decvec/internal/dvec"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

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

var targetPath string
var inputFastaPath string
var useAddgene bool

// set flags
func init() {
	rootCmd.AddCommand(makeCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	makeCmd.Flags().StringVarP(&targetPath, "target", "t", "", "path to a FASTA file with the target vector")
	makeCmd.Flags().BoolVarP(&useAddgene, "addgene", "a", false, "flag signalling we should use addgene's vector db")

	// Mark required flags
	makeCmd.MarkFlagRequired("target")
	makeCmd.MarkFlagRequired("addgene")

	// Bind the paramters to viper
	viper.BindPFlag("make.target", makeCmd.Flags().Lookup("target"))
	viper.BindPFlag("make.addgene", makeCmd.Flags().Lookup("addgene"))
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
	c := config.New()

	// no path to input file
	if c.Make.TargetPath == "" {
		log.Fatal("Failed, no target fragment path set")
	}

	// read in fragments
	fragments, err := dvec.Read(c.Make.TargetPath)
	if err != nil {
		log.Fatalf("Failed to read in fasta files at %s: %v", c.Make.TargetPath, err)
	}

	// set target fragment
	if len(fragments) > 1 {
		println(
			"Warning: %d building fragments were in %s. Only targeting the first: %s",
			len(fragments),
			c.Make.TargetPath,
			fragments[0].ID,
		)
	}

	// target := fragments[0]

	// dbPath, err := filepath.Abs(path.Join("..", "..", "assets", "addgene", "db", "addgene"))
	// handle(err)
}
