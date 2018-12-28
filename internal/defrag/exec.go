package defrag

import (
	"fmt"
	"io/ioutil"
	"log"
	"path"
	"path/filepath"
	"regexp"
	"strings"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
)

// Vector is the root of the `defrag vector` functionality
//
// the goal is to find an "optimal" assembly vector with:
// 	1. the fewest fragments
// 	2. the lowest overall assembly cost ($)
// and, secondarily:
//	3. no duplicate end regions between Gibson fragments
// 	4. no inverted repeats in the junctions
// 	5. no off-target binding sites in the parent vectors
//	6. low primer3 penalty scores
func Vector(cmd *cobra.Command, args []string) {
	in, err := cmd.Flags().GetString("in")
	if err != nil {
		in = guessInput()
	}

	out, err := cmd.Flags().GetString("out")
	if err != nil {
		out = guessOutput(in)
	}

	addgene, err := cmd.Flags().GetBool("addgene")
	if err != nil {
		log.Fatal(err)
	}

	dbs, err := cmd.Flags().GetString("dbs")
	if err != nil && !addgene {
		log.Fatalf("failed to find dbs of building fragments: %v", err)
	}

	execVector(in, out, dbs, addgene)
}

// execVector accepts an input path with a vector to build, an output path
// to write solutions to, dbs of possible building fragments and
// a flag for whether to use addgene as a building fragment source
func execVector(in, out, dbs string, addgene bool) [][]Fragment {
	conf := config.New()

	// no path to input file
	if in == "" {
		log.Fatal("failed without an input file argument")
	}

	// read the target sequence (the first in the slice is used)
	fragments, err := read(in)
	if err != nil {
		log.Fatalf("failed to read in fasta files at %s: %v", in, err)
	}

	// set target fragment
	if len(fragments) > 1 {
		log.Printf(
			"warning: %d building fragments were in %s. Only targeting the first: %s",
			len(fragments),
			in,
			fragments[0].ID,
		)
	}
	targetFrag := fragments[0]

	// read in the BLAST DB paths
	dbPaths, err := getDBs(dbs, config.Root, addgene)
	if err != nil {
		log.Fatalf("failed to find a fragment database: %v", err)
	}

	// get all the matches against the fragment
	matches, err := blast(&targetFrag, dbPaths, conf.Fragments.MinHomology, conf.Vendors())
	if err != nil {
		log.Fatalf("failed to blast %s against the BLAST DB: %v", targetFrag.ID, err)
	}

	// build up the assemblies
	builds := assembleREV(matches, targetFrag.Seq, &conf)

	// try to write the JSON to the output file path
	if !filepath.IsAbs(out) {
		if out, err = filepath.Abs(out); err != nil {
			log.Fatalf("failed to make output path absolute: %v", err)
		}
	}
	write(out, targetFrag, builds)

	// return the builds (for e2e testing)
	return builds
}

// Fragments accepts a cobra.Command with flags for assembling a list of
// fragments together into a vector (in the order specified). Fragments
// without junctions for their neighbors are prepared via PCR
func Fragments(cmd *cobra.Command, args []string) {
	in, err := cmd.Flags().GetString("in")
	if err != nil {
		in = guessInput()
	}

	out, err := cmd.Flags().GetString("out")
	if err != nil {
		out = guessOutput(in)
	}

	execFragments(in, out)
}

// execFragments takes an input file with a list of fragments to assemble,
// pieces them together (preparing them if necessary) and writing the results
// with primers to the out file
func execFragments(in, out string) {
	c := config.New()

	// read the target sequence (the first in the slice is used)
	inputFragments, err := read(in)
	if err != nil {
		log.Fatalf("failed to read in fasta files at %s: %v", in, err)
	}

	// try to find the target vector (sequence) and prepare the fragments to
	// assemble it
	target, fragments := assembleFWD(inputFragments, &c)

	// write the single list of fragments as a possible solution to the output file
	write(out, target, [][]Fragment{fragments})
}

// guessInput returns the first fasta file in the current directory. Is used
// if the user hasn't specified an input file
func guessInput() (in string) {
	dir, _ := filepath.Abs(".")
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		log.Fatal(err)
	}
	for _, file := range files {
		if file.IsDir() {
			continue
		}

		ext := strings.ToLower(filepath.Ext(file.Name()))
		if ext == ".fa" || ext == ".fasta" {
			return file.Name()
		}
	}

	log.Fatalf("failed: no input argument set and no fasta file found in %s", dir)
	return
}

// guessOutput gets an outpath path from an input path (if no output path is
// specified). It uses the same name as the input path to create an output
func guessOutput(in string) (out string) {
	ext := filepath.Ext(in)
	noExt := in[0 : len(in)-len(ext)]
	return noExt + ".defrag.json"
}

// getDBs returns a list of absolute paths to BLAST databases
func getDBs(dbsInput, root string, addgene bool) (paths []string, err error) {
	if addgene {
		addgenePath := path.Join(root, "assets", "addgene", "db", "addgene")
		return parseDBs(dbsInput + "," + addgenePath)
	}

	return parseDBs(dbsInput)
}

// parseDBs turns a single string of comma separated BLAST dbs into a
// slice of absolute paths to the BLAST dbs on the local fs
func parseDBs(dbList string) (paths []string, err error) {
	re := regexp.MustCompile(",")

	for _, db := range re.Split(dbList, -1) {
		dbTrimmed := strings.Trim(db, " ,")
		if dbTrimmed == "" {
			continue
		}

		absPath, err := filepath.Abs(dbTrimmed)

		if err != nil {
			return nil, fmt.Errorf("failed to create absolute path: %v", err)
		}
		paths = append(paths, absPath)
	}
	return
}
