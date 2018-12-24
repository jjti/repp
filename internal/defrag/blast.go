package defrag

import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"runtime"
	"sort"
	"strconv"
	"strings"

	"github.com/jjtimmons/defrag/config"
)

// blastExec is a small utility function for executing BLAST on a fragment
type blastExec struct {
	// the fragment we're BLASTing
	f *Fragment

	// the path to the database we're BLASTing against
	db string

	// the path to the input BLAST file
	in string

	// the path for the BLAST output
	out string

	// optional path to a FASTA file with a subject FASTA sequence
	subject string

	// path to the blastn executable
	blastn string
}

// blast the passed Fragment against a set from the command line and create
// matches for those that are long enough
//
// Accepts a fragment to BLAST against, a list of dbs to BLAST it against,
// a minLength for a match, and settings around blastn location, output dir, etc
func blast(f *Fragment, dbs []string, minLength int, v config.VendorConfig) (matches []Match, err error) {
	for _, db := range dbs {
		b := &blastExec{
			f:      f,
			db:     db,
			in:     path.Join(v.Blastdir, f.ID+".input.fa"),
			out:    path.Join(v.Blastdir, f.ID+".output"),
			blastn: v.Blastn,
		}

		// make sure the db exists
		if _, err := os.Stat(db); os.IsNotExist(err) {
			return nil, fmt.Errorf("failed to find an Addgene database at %s", db)
		}

		// create the input file
		if err := b.create(); err != nil {
			return nil, fmt.Errorf("failed at creating BLAST input file at %s: %v", b.in, err)
		}

		// execute BLAST
		if err := b.run(); err != nil {
			return nil, fmt.Errorf("failed executing BLAST: %v", err)
		}

		// parse the output file to Matches against the Fragment
		dbMatches, err := b.parse()
		if err != nil {
			return nil, fmt.Errorf("failed to parse BLAST output: %v", err)
		}

		log.Printf("%d matches in %s\n", len(dbMatches), db)

		// add these matches against the growing list of matches
		matches = append(matches, dbMatches...)
	}

	// keep only "proper" arcs (non-self-contained)
	matches = filter(matches, minLength)
	if len(matches) < 1 {
		return nil, fmt.Errorf("did not find any matches for %s", f.ID)
	}

	for _, m := range matches {
		fmt.Println(m.Entry, m.Start, m.End)
	}

	return matches, err
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

// run calls the external blastn binary on the input library
func (b *blastExec) run() error {
	threads := runtime.NumCPU() - 1
	if threads < 1 {
		threads = 1
	}

	// create the blast command
	// https://www.ncbi.nlm.nih.gov/books/NBK279682/
	blastCmd := exec.Command(
		b.blastn,
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
		log.Fatalf("failed to execute blastn against db, %s: %v: %s", b.db, err, string(output))
	}

	return nil
}

// runs blast on the query file against another subject file (rather than blastdb)
func (b *blastExec) runAgainst() error {
	// create the blast command
	// https://www.ncbi.nlm.nih.gov/books/NBK279682/
	blastCmd := exec.Command(
		b.blastn,
		"-task", "blastn",
		"-query", b.in,
		"-subject", b.subject,
		"-out", b.out,
		"-outfmt", "7 sseqid qstart qend sstart send sseq mismatch",
	)

	// execute BLAST and wait on it to finish
	if output, err := blastCmd.CombinedOutput(); err != nil {
		log.Fatalf("failed to execute blastn against db, %s: %v: %s", b.db, err, string(output))
		return err
	}

	return nil
}

// parse reads the output file into Matches on the Fragment
// returns a slice of Matches for the blasted fragment
func (b *blastExec) parse() (matches []Match, err error) {
	// read in the results
	file, err := ioutil.ReadFile(b.out)
	if err != nil {
		return
	}
	fileS := string(file)

	// read it into Matches
	var ms []Match
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
		ms = append(ms, Match{
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

// filter "proper-izes" the matches from BLAST
//
// TODO: filter on multiple BLAST databases. Keep ALL fragments from non-remote repositories (ie free
// fragments) but still properize out the remote fragments. So if a large fragment in the local repository
// completely emcompasses a smaller fragment in Addgene, remove the fragment from the Addgene database.
// Cannot do the same if there's a large fragment in the remote but not local db
//
// proper-izing fragment matches means removing those that are completely
// self-contained in other fragments: the larger of the available fragments
// will be the better one, since it covers a greater region and will almost
// always be preferable to the smaller one
//
// Circular-arc graph: https://en.wikipedia.org/wiki/Circular-arc_graph
//
// also remove small fragments here, that are too small to be useful during assembly
func filter(matches []Match, minSize int) (properized []Match) {
	properized = []Match{}

	// remove fragments that are shorter the minimum cut off size
	var largeEnough []Match
	for _, m := range matches {
		if m.Length() > minSize {
			largeEnough = append(largeEnough, m)
		}
	}

	// sort largeEnough by their start index
	// for fragments with equivelant starting indexes, put the larger one first
	sort.Slice(largeEnough, func(i, j int) bool {
		if largeEnough[i].Start != largeEnough[j].Start {
			return largeEnough[i].Start < largeEnough[j].Start
		}
		return largeEnough[i].Length() > largeEnough[j].Length()
	})

	// only include those that aren't encompassed by the one before it
	for _, m := range largeEnough {
		lastMatch := len(properized) - 1
		if lastMatch < 0 || m.End > properized[lastMatch].End {
			properized = append(properized, m)
		}
	}

	return properized
}
