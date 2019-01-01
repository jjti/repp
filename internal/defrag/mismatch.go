package defrag

import (
	"fmt"
	"io/ioutil"
	"log"
	"path"
	"path/filepath"
	"strings"

	"github.com/jjtimmons/defrag/config"
)

// parentMismatch both searches for a the parent fragment in its source DB and queries for
// any mismatches in the seq before returning
func parentMismatch(primers []Primer, parent, db string, conf *config.Config) (wasMismatch bool, m match, err error) {
	// try and query for the parent in the source DB and write to a file
	parentFile, err := blastdbcmd(parent, db, conf)

	// ugly check here for whether we just failed to get the parent entry from a db
	// which isn't a huge deal (shouldn't be flagged as a mismatch)
	// this is similar to what io.IsNotExist does
	if err != nil {
		if strings.Contains(err.Error(), "failed to query") {
			log.Println(err) // just write the error
			// TODO: if we fail to find the parent, query the fullSeq as it was sent
			return false, match{}, nil
		}
		return false, match{}, err
	}

	// check each primer for mismatches
	if parentFile != "" {
		for _, primer := range primers {
			wasMismatch, m, err = mismatch(primer.Seq, parentFile, conf)
			if wasMismatch || err != nil {
				return
			}
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

	return false, match{}, nil
}

// mismatch finds mismatching sequences between the query sequence and
// the parent sequence (in the parent file)
//
// The fragment to query against is stored in a file ".parent"
// db is passed as the path to the db we're blasting against
func mismatch(primer, parentFile string, c *config.Config) (wasMismatch bool, m match, err error) {
	v := c.Vendors()

	// path the query sequence input file
	in := parentFile + ".primer.query"

	// path to the blastOutput file
	out := parentFile + ".blast"

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
