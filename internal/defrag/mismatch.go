package defrag

import (
	"fmt"
	"io/ioutil"
	"log"
	"os/exec"
	"path"
	"path/filepath"
	"strings"

	"github.com/jjtimmons/defrag/config"
)

// parentMismatch both searches for a the parent fragment in its source DB and queries for
// any mismatches in the seq before returning
func parentMismatch(primers []Primer, parent, db string, conf *config.Config) (wasMismatch bool, m match, err error) {
	// try and query for the parent in the source DB and write to a file
	parentFile := blastDBCmd(parent, db, conf)
	if parentFile == "" {
		return false, match{}, nil
	}

	// check each primer for mismatches
	for _, primer := range primers {
		wasMismatch, m, err = mismatch(primer.Seq, parentFile, conf)
		if wasMismatch || err != nil {
			return
		}
	}
	return
}

// seqMismatch queries for any mismatching primer locations in the parent sequence
// unlike parentMismatch, it doesn't first find the parent fragment from the db it came from
func seqMismatch(primers []Primer, parentID, parentSeq string, conf *config.Config) (wasMismatch bool, m match, err error) {
	v := conf.Vendors()

	// write a parent seq file to BLAST against
	parentFile, _ := filepath.Abs(path.Join(v.Blastdir, parentID+".parent"))
	inContent := fmt.Sprintf(">%s\n%s\n", parentID, parentSeq)
	if err = ioutil.WriteFile(parentFile, []byte(inContent), 0666); err != nil {
		return false, m, fmt.Errorf("failed to write primer sequence to query FASTA file: %v", err)
	}

	// check each primer for mismatches
	for _, primer := range primers {
		wasMismatch, m, err = mismatch(primer.Seq, parentFile, conf)
		if wasMismatch || err != nil {
			return
		}
	}
	return
}

// blastDBCmd gets the full parent's sequence from the DB that it came from
// returns an empty string if we can't find the sequence
//
// parent here is the id that's associated with the fragment in its source DB
func blastDBCmd(parent, db string, c *config.Config) (parentFile string) {
	v := c.Vendors()

	// path to the entry batch file to hold the parent entry accession
	entry, _ := filepath.Abs(path.Join(v.Blastdir, parent+".parent.entry"))

	// path to the output sequence file from querying the parent's sequence from the BLAST db
	parentPath, _ := filepath.Abs(path.Join(v.Blastdir, parent+".parent"))

	// write entry to file
	// this was a 2-day bug I couldn't resolve...
	// I was using the "-entry" flag on exec.Command, but have since
	// switched to the simpler -entry_batch command (on a file) that resolves the issue
	if err := ioutil.WriteFile(entry, []byte(parent), 0666); err != nil {
		log.Fatalf("failed to write batch entry list: %v", err)
		return
	}

	// make a blastdbcmd command (for querying a DB, very different from blastn)
	queryCmd := exec.Command(
		v.Blastdbcmd,
		"-db", db,
		"-dbtype", "nucl",
		"-entry_batch", entry,
		"-out", parentPath,
		"-outfmt", "%f", // fasta format
	)
	if _, err := queryCmd.CombinedOutput(); err != nil {
		log.Printf("warning: failed to query %s from %s", parent, db)
		return "" // failed to get parent sequence from the db
	}

	// read in the results as a fragment and return just the seq
	fragments, err := read(parentPath)
	if err != nil && len(fragments) >= 1 {
		return parentPath
	}
	return
}

// mismatch finds mismatching sequences between the query sequence and
// the parent sequence (in the parent file)
//
// The fragment to query against is stored in a file ".parent"
// db is passed as the path to the db we're blasting against
func mismatch(primer, parentFile string, c *config.Config) (wasMismatch bool, m match, err error) {
	v := c.Vendors()

	// path the query sequence input file
	in, _ := filepath.Abs(path.Join(v.Blastdir, parentFile+".primer.query"))

	// path to the blastOutput file
	out, _ := filepath.Abs(path.Join(v.Blastdir, parentFile+".blast"))

	// create blast input file
	inContent := fmt.Sprintf(">primer\n%s\n", primer)
	if err = ioutil.WriteFile(in, []byte(inContent), 0666); err != nil {
		return false, m, fmt.Errorf("failed to write primer sequence to query FASTA file: %v", err)
	}

	// blast the query sequence against the parentFile sequence
	b := blastExec{
		in:      in,
		out:     out,
		subject: parentFile,
		blastn:  v.Blastn,
	}

	// execute blast
	if err = b.runAgainst(); err != nil {
		return false, m, fmt.Errorf("failed to run blast against parent: %v", err)
	}

	// get the BLAST matches
	matches, err := b.parse()
	if err != nil {
		return false, match{}, fmt.Errorf("failed to parse matches from %s: %v", out, err)
	}

	// parse the results and check whether any are cause for concern (by Tm)
	primerCount := 1 // times we expect to see the primer itself
	for i, m := range matches {
		if i == 0 && m.circular {
			// if the match is against a circular fragment, we might expect to see
			// the primer's sequence twice, rather than just once
			primerCount++
		}

		// one of the matches will, of course, be against the primer itself
		// and we don't want to double count it
		if primerCount > 0 && m.seq == primer {
			primerCount--
			continue
		} else if isMismatch(m, c) {
			return true, m, nil
		}
	}

	return false, match{}, nil
}

// isMismatch reutns whether the match constitutes a mismatch
// between it and the would be primer sequence
//
// source: http://depts.washington.edu/bakerpg/primertemp/
//
// The equation used for the melting temperature is:
// Tm = 81.5 + 0.41(%GC) - 675/N - % mismatch, where N = total number of bases.
func isMismatch(m match, c *config.Config) bool {
	primer := strings.ToLower(m.seq)
	primerL := float64(len(primer))

	noA := strings.Replace(primer, "a", "", -1)
	noT := strings.Replace(noA, "t", "", -1)
	gcPerc := float64(len(noT)) / primerL
	tmNoMismatch := 81.5 + 0.41*gcPerc - 675/float64(len(primer))
	tmWithMismatch := tmNoMismatch - float64(m.mismatching)/primerL

	return tmWithMismatch > c.PCR.MaxOfftargetTm
}
