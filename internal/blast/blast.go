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
	"strconv"
	"strings"

	"github.com/jjtimmons/decvec/config"
	"github.com/jjtimmons/decvec/internal/dvec"
)

var (
	// config object
	conf = config.New()

	// path to the blast directory for putting results into
	blastDir = filepath.Join(conf.Root, "bin", "blast")

	// path to the blast executable
	blast = path.Join(conf.Root, "vendor", "ncbi-blast-2.7.1+", "bin", "blastn")
)

// blastExec is a small utility function for executing BLAST
// on a fragment.
type blastExec struct {
	// the fragment we're BLASTing
	f *dvec.Fragment

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
// Accepts a fragment to blast against
func BLAST(f *dvec.Fragment) (matches []dvec.Match, err error) {
	b := &blastExec{
		f:   f,
		in:  path.Join(blastDir, f.ID+".input.fa"),
		out: path.Join(blastDir, f.ID+".output"),
	}

	// make sure the addgene db is there
	if _, err := os.Stat(conf.DB); os.IsNotExist(err) {
		return nil, fmt.Errorf("failed to find an Addgene database at %s", conf.DB)
	}

	// make sure the blast executable is there
	if _, err := os.Stat(blast); os.IsNotExist(err) {
		return nil, fmt.Errorf("failed to find a BLAST executable at %s", blast)
	}

	// create the input file
	if err := b.create(); err != nil {
		return nil, fmt.Errorf("failed at creating BLAST input file at %s: %v", b.in, err)
	}

	// execute BLAST on it
	if err := b.run(); err != nil {
		return nil, fmt.Errorf("failed at executing BLAST: %v", err)
	}

	// parse the BLAST output file to Matches for the Fragment
	matches, err = b.parse()
	if err != nil {
		return nil, fmt.Errorf("failed to parse BLAST output: %v", err)
	} else if len(matches) < 1 {
		return nil, fmt.Errorf("did not find any matches for %s", b.f.ID)
	}
	return matches, err
}

// input creates an input file for BLAST
// return the path to the file and an error if there was one
func (b *blastExec) create() error {
	// create the query sequence file.
	// add the sequence to itself because it's circular
	// and we want to find matches across the zero-index.
	file := fmt.Sprintf(">%s\n%s%s\n", b.f.ID, b.f.Seq, b.f.Seq)
	return ioutil.WriteFile(b.in, []byte(file), 0666)
}

// run calls the external blastn binary
func (b *blastExec) run() error {
	// create the blast command
	// https://www.ncbi.nlm.nih.gov/books/NBK279682/
	blastCmd := exec.Command(
		blast,
		"-task", "blastn",
		"-db", conf.DB,
		"-query", b.in,
		"-out", b.out,
		"-outfmt", "7 sseqid qstart qend sstart send sseq mismatch",
		"-ungapped",
		"-perc_identity", "100",
		"-max_target_seqs", "100000",
	)

	// execute BLAST and wait on it to finish
	if output, err := blastCmd.CombinedOutput(); err != nil {
		fmt.Fprintln(os.Stderr, err, string(output))
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
		"-ungapped",
	)

	// execute BLAST and wait on it to finish
	if output, err := blastCmd.CombinedOutput(); err != nil {
		fmt.Fprintln(os.Stderr, err, string(output))
		return err
	}
	return nil
}

// parse reads the output file into Matches on the Fragment
// returns a slice of Matches for the blasted fragment
func (b *blastExec) parse() (matches []dvec.Match, err error) {
	// read in the results
	file, err := ioutil.ReadFile(b.out)
	if err != nil {
		return
	}
	fileS := string(file)

	// read it into Matches
	var ms []dvec.Match
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

		// direction not guarenteed
		if start > end {
			start, end = end, start
		}

		// create and append the new match
		ms = append(ms, dvec.Match{
			// for later querying when checking for off-targets
			Entry: id,
			Seq:   seq,
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

	err := os.MkdirAll(blastDir, os.ModePerm)
	if err != nil {
		log.Fatalf("failed to create a BLAST dir: %v", err)
	}
}
