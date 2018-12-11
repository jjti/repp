package exec

import (
	"log"
	"path/filepath"

	"github.com/jjtimmons/defrag/config"
	"github.com/jjtimmons/defrag/internal/assemble"
	"github.com/jjtimmons/defrag/internal/blast"
	"github.com/jjtimmons/defrag/internal/io"
	"github.com/spf13/cobra"
)

// Execute is the root of the make functionality
//
// the goal is to find an "optimal" assembly vector with:
// 	1. fewest fragments
// 	2. lowest overall assembly cost ($)
// and, secondarily:
//	3. no duplicate end regions between Gibson fragments
// 	4. no inverted repeats in the junctions
// 	5. no off-target binding sites in the parent vectors
//	6. low primer3 penalty scores
func Execute(cmd *cobra.Command, args []string) {
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

	// read in the BLAST DB paths from config
	dbPaths, err := conf.DBList()
	if err != nil {
		log.Fatalf("Failed to find a BLAST database: %v", err)
	}

	// get all the matches against the fragment
	matches, err := blast.BLAST(&targetFrag, dbPaths, conf.Fragments.MinHomology, conf.Vendors())
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
}
