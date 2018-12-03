// Package blast finds Matches between a Fragment and the Fragment databases.
package blast

import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"

	"github.com/jjtimmons/defrag/config"
	"github.com/jjtimmons/defrag/internal/defrag"
)

var (
	// config object
	conf = config.New()

	// path to the blast executable
	blast = path.Join(conf.Root, "vendor", "ncbi-blast-2.7.1+", "bin", "blastn")
)

// blastExec is a small utility function for executing BLAST
// on a fragment.
type blastExec struct {
	// the fragment we're BLASTing
	f *defrag.Fragment

	// the path to the database we're BLASTing against
	db string

	// the path to the input BLAST file
	in string

	// the path for the BLAST output
	out string

	// optional path to a FASTA file with a subject FASTA sequence
	subject string
}

// BLAST the passed Fragment against a set from the command line and create
// matches for those that are long enough
//
// Accepts a fragment to BLAST against
func BLAST(f *defrag.Fragment, dbs, dir string, minLength int) (matches []defrag.Match, err error) {
	parsedDBs, err := parseDBs(dbs)
	if err != nil {
		return
	}

	for _, db := range parsedDBs {
		b := &blastExec{
			f:   f,
			db:  db,
			in:  path.Join(dir, f.ID+".input.fa"),
			out: path.Join(dir, f.ID+".output"),
		}

		// make sure the db exists
		if _, err := os.Stat(db); os.IsNotExist(err) {
			return nil, fmt.Errorf("failed to find an Addgene database at %s", db)
		}

		// make sure the blast executable exists
		if _, err := os.Stat(blast); os.IsNotExist(err) {
			return nil, fmt.Errorf("failed to find a BLAST executable at %s", blast)
		}

		// create the input file
		if err := b.create(); err != nil {
			return nil, fmt.Errorf("failed at creating BLAST input file at %s: %v", b.in, err)
		}

		// execute BLAST
		if err := b.run(); err != nil {
			return nil, fmt.Errorf("failed at executing BLAST: %v", err)
		}

		// parse the output file to Matches against the Fragment
		dbMatches, err := b.parse()
		if err != nil {
			return nil, fmt.Errorf("failed to parse BLAST output: %v", err)
		}

		// add these matches against the growing list of matches
		matches = append(matches, dbMatches...)
	}

	// keep only "proper" arcs (non-self-contained)
	matches = filter(matches, minLength)
	if len(matches) < 1 {
		return nil, fmt.Errorf("did not find any matches for %s", f.ID)
	}
	return matches, err
}

// parseDBs turns a single string of comma separated BLAST dbs into a
// slice of absolute paths to the BLAST dbs on the local fs
func parseDBs(dbList string) (paths []string, err error) {
	noSpaceDBs := strings.Replace(dbList, " ", "", -1)
	for _, db := range strings.Split(noSpaceDBs, ",") {
		absPath, err := filepath.Abs(db)

		if err != nil {
			return nil, fmt.Errorf("failed to create absolute path: %v", err)
		}

		paths = append(paths, absPath)
	}

	return
}

// input creates an input file for BLAST
// return the path to the file and an error if there was one
func (b *blastExec) create() error {
	// create the query sequence file.
	// add the sequence to itself because it's circular
	// and we want to find matches across the zero-index.
	file := fmt.Sprintf(">%s\n%s\n", b.f.ID, b.f.Seq+b.f.Seq)
	return ioutil.WriteFile(b.in, []byte(file), 0666)
}

// run calls the external blastn binary
func (b *blastExec) run() error {
	threads := runtime.NumCPU() - 1
	if threads < 1 {
		threads = 1
	}

	// create the blast command
	// https://www.ncbi.nlm.nih.gov/books/NBK279682/
	blastCmd := exec.Command(
		blast,
		"-task", "blastn",
		"-db", b.db,
		"-query", b.in,
		"-out", b.out,
		"-outfmt", "7 sseqid qstart qend sstart send sseq mismatch",
		"-perc_identity", "100",
		"-num_threads", strconv.Itoa(threads),
	)

	// execute BLAST and wait on it to finish
	if output, err := blastCmd.CombinedOutput(); err != nil {
		log.Fatalf("Failed to execute BLAST: %v: %s", err, string(output))
		return err
	}
	return nil
}

// runs blast on the query file against another subject file
// (rather than the blastdb)
func (b *blastExec) runAgainst() error {
	// create the blast command
	// https://www.ncbi.nlm.nih.gov/books/NBK279682/
	blastCmd := exec.Command(
		blast,
		"-task", "blastn",
		"-query", b.in,
		"-subject", b.subject,
		"-out", b.out,
		"-outfmt", "7 sseqid qstart qend sstart send sseq mismatch",
	)

	// execute BLAST and wait on it to finish
	if output, err := blastCmd.CombinedOutput(); err != nil {
		log.Fatalf("Failed to execute BLAST: %v: %s", err, string(output))
		return err
	}
	return nil
}

// parse reads the output file into Matches on the Fragment
// returns a slice of Matches for the blasted fragment
func (b *blastExec) parse() (matches []defrag.Match, err error) {
	// read in the results
	file, err := ioutil.ReadFile(b.out)
	if err != nil {
		return
	}
	fileS := string(file)

	// read it into Matches
	var ms []defrag.Match
	for _, line := range strings.Split(fileS, "\n") {
		// comment lines start with a #
		if strings.HasPrefix(line, "#") {
			continue
		}

		// split on white space
		cols := strings.Fields(line)
		if len(cols) < 6 {
			continue
		}

		// the full id of the entry in the db
		id := strings.Replace(cols[0], ">", "", -1)

		start, _ := strconv.Atoi(cols[1])
		end, _ := strconv.Atoi(cols[2])
		seq := cols[5]
		mismatch, _ := strconv.Atoi(cols[6])

		// direction not gaurenteed
		// TODO: how is guarentee spelled?...
		if start > end {
			start, end = end, start
		}

		// create and append the new match
		ms = append(ms, defrag.Match{
			// for later querying when checking for off-targets
			Entry: id,
			Seq:   strings.Replace(seq, "-", "", -1),
			// convert 1-based numbers to 0-based
			Start: start - 1,
			End:   end - 1,
			// brittle, but checking for circular in entry's id
			Circular: strings.Contains(id, "(circular)"),
			Mismatch: mismatch,
		})
	}
	return ms, nil
}

// init ensures there's a blast subdirectory in the binary's execution enviornment
// for the results this is about to create
func init() {
	err := os.MkdirAll(conf.BlastDir, os.ModePerm)
	if err != nil {
		log.Fatalf("Failed to create a BLAST dir: %v", err)
	}
}
