package cmd

import (
	"log"
	"path/filepath"

	"github.com/jjtimmons/defrag/config"
	"github.com/jjtimmons/defrag/internal/assemble"
	"github.com/jjtimmons/defrag/internal/blast"
	"github.com/jjtimmons/defrag/internal/io"

	"github.com/spf13/viper"

	"github.com/spf13/cobra"
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
	Run: makeExec,
}

// set flags
func init() {
	rootCmd.AddCommand(makeCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	makeCmd.PersistentFlags().StringP("target", "t", "", "Input file name of target vector sequence <FASTA>")
	makeCmd.PersistentFlags().StringP("dbs", "d", "", "Comma separated list of building fragment databases (multi-fasta)")
	makeCmd.PersistentFlags().StringP("out", "o", "", "Output file name")
	makeCmd.PersistentFlags().BoolP("addgene", "a", false, "Use the Addgene repository as a source of building fragments")

	// db path is needed globally
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
	conf := config.New()

	target, err := cmd.PersistentFlags().GetString("target")
	if err != nil {
		log.Fatalf("Cannot get target from arguments: %v", err)
	}

	// no path to input file
	if target == "" {
		log.Fatal("Failed, no target fragment path set")
	}

	output, err := cmd.PersistentFlags().GetString("out")
	if err != nil {
		log.Fatalf("Cannot find the output path: %v", err)
	}

	// read in fragments, the first is the target sequence
	fragments, err := io.Read(target)
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
	targetFrag := fragments[0]

	// get all the matches against the fragment
	matches, err := blast.BLAST(&targetFrag, conf.DB, conf.BlastDir, conf.PCR.MinLength)
	if err != nil {
		log.Fatalf("Failed to blast %s against the BLAST DB: %v", targetFrag.ID, err)
	}

	// build up the assemblies
	builds := assemble.Assemble(matches, targetFrag.Seq)

	// try to write the JSON to the filepath
	if !filepath.IsAbs(output) {
		output, err = filepath.Abs(output)
		if err != nil {
			log.Fatalf("Failed to make output path absolute: %v", err)
		}
	}
	io.Write(output, targetFrag, builds)

	// os.Exit(0)
}
