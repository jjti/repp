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

// match is a blast "hit" in the blastdb
type match struct {
	// entry of the matched fragment in the database
	entry string

	// seq of the match on the target vector
	seq string

	// start of the fragment (0-indexed)
	start int

	// end of the fragment (0-indexed)
	end int

	// circular if it's a circular fragment (vector, plasmid, etc)
	circular bool

	// mismatching number of bps in the match (for primer off-targets)
	mismatching int

	// internal if the fragment doesn't have to be procured from a remote repository (eg Addgene, iGEM)
	internal bool
}

// Length returns the length of the match on the target fragment
func (m *match) Length() int {
	return m.end - m.start + 1 // it's inclusive
}

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

	// internal if the db is a local/user owned list of fragments (ie free)
	internal bool
}

// blast the passed Fragment against a set from the command line and create
// matches for those that are long enough
//
// Accepts a fragment to BLAST against, a list of dbs to BLAST it against,
// a minLength for a match, and settings around blastn location, output dir, etc
func blast(f *Fragment, dbs []string, minLength int, v config.VendorConfig) (matches []match, err error) {
	for _, db := range dbs {
		internal := true
		if strings.Contains(db, "addgene") {
			internal = false
		}

		b := &blastExec{
			f:        f,
			db:       db,
			in:       path.Join(v.Blastdir, f.ID+".input.fa"),
			out:      path.Join(v.Blastdir, f.ID+".output"),
			blastn:   v.Blastn,
			internal: internal,
		}

		// make sure the db exists
		if _, err := os.Stat(db); os.IsNotExist(err) {
			return nil, fmt.Errorf("failed to find a BLAST database at %s", db)
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
func (b *blastExec) parse() (matches []match, err error) {
	// read in the results
	file, err := ioutil.ReadFile(b.out)
	if err != nil {
		return
	}
	fileS := string(file)

	// read it into Matches
	var ms []match
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
		ms = append(ms, match{
			// for later querying when checking for off-targets
			entry: id,
			seq:   strings.Replace(seq, "-", "", -1),
			// convert 1-based numbers to 0-based
			start: start - 1,
			end:   end - 1,
			// brittle, but checking for circular in entry's id
			circular:    strings.Contains(id, "(circular)"),
			mismatching: mismatch,
			internal:    b.internal,
		})
	}
	return ms, nil
}

// filter "proper-izes" the matches from BLAST
//
// TODO: filter further here, can remove external matches that are
// entirely contained by internal matches but am not doing that here
//
// proper-izing fragment matches means removing those that are completely
// self-contained in other fragments: the larger of the available fragments
// will be the better one, since it covers a greater region and will almost
// always be preferable to the smaller one
//
// Circular-arc graph: https://en.wikipedia.org/wiki/Circular-arc_graph
//
// also remove small fragments here, that are too small to be useful during assembly
func filter(matches []match, minSize int) (properized []match) {
	properized = []match{}

	// remove fragments that are shorter the minimum cut off size
	// separate the internal and external fragments. the internal
	// ones should not be removed just if they're self-contained
	// in another, because they may be cheaper to assemble
	var internal []match
	var external []match
	for _, m := range matches {
		if m.Length() > minSize {
			if m.internal {
				internal = append(internal, m)
			} else {
				external = append(external, m)
			}
		}
	}

	return append(properize(internal), properize(external)...)
}

// properize remove matches that are entirely contained within others
func properize(matches []match) (properized []match) {
	// sort matches by their start index
	// for fragments with equivelant starting indexes, put the larger one first
	sort.Slice(matches, func(i, j int) bool {
		if matches[i].start != matches[j].start {
			return matches[i].start < matches[j].start
		}
		return matches[i].Length() > matches[j].Length()
	})

	// only include those that aren't encompassed by the one before it
	for _, m := range matches {
		lastMatch := len(properized) - 1
		if lastMatch < 0 || m.end > properized[lastMatch].end {
			properized = append(properized, m)
		}
	}
	return
}
