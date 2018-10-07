package blast

import (
	"fmt"
	"io/ioutil"
	"os/exec"
	"path"
	"path/filepath"
	"strings"

	"github.com/jjtimmons/decvec/internal/frag"
)

// isMismatch reutns whether the match constitutes a mismatch
// between it and the would be primer sequence
//
// source: http://depts.washington.edu/bakerpg/primertemp/
//
// The equation used for the melting temperature is:
// Tm = 81.5 + 0.41(%GC) - 675/N - % mismatch, where N = total number of bases.
func isMismatch(match frag.Match) bool {
	primer := strings.ToLower(match.Seq)
	primerL := float32(len(primer))

	noA := strings.Replace(primer, "a", "", -1)
	noT := strings.Replace(noA, "t", "", -1)
	gcPerc := float32(len(noT)) / primerL
	tmNoMismatch := 81.5 + 0.41*gcPerc - 675/float32(len(primer))
	tmWithMismatch := tmNoMismatch - float32(match.Mismatch)/primerL

	return tmWithMismatch > 45 // TODO: move to settings
}

// Mismatch finds mismatching sequences between the query sequence and
// the parent sequence
//
// The parent sequence is passed as the entry id as it exists in the
// blast database
func Mismatch(primer, parent string) (mismatch bool, match frag.Match, err error) {
	// path to the blastdbcmd binary
	blastcmd, _ := filepath.Abs(path.Join("..", "..", "vendor", "ncbi-blast-2.7.1+", "bin", "blastdbcmd"))
	// path to the entry batch file to hold the parent entry accession
	entry, _ := filepath.Abs(path.Join(blastDir, parent+".entry"))
	// path to the output sequence file
	entryOutput, _ := filepath.Abs(path.Join(blastDir, parent+".out"))
	// path the query sequence input file
	in, _ := filepath.Abs(path.Join(blastDir, primer+"-"+parent+".query"))
	// path to the blastOutput file
	out, _ := filepath.Abs(path.Join(blastDir, parent+".blast"))

	// write entry to file
	// this was a 2-day bug I couldn't resolve...
	// I was using the "-entry" flag on exec.Command, but have since
	// switched to the far simpler -entry_batch command that fixes the issue
	err = ioutil.WriteFile(entry, []byte(parent), 0666)
	if err != nil {
		return false, match, fmt.Errorf("failed to write batch entry list: %v", err)
	}

	// make a blast command
	queryCmd := exec.Command(
		blastcmd,
		"-db", db,
		"-dbtype", "nucl",
		"-entry_batch", entry,
		"-out", entryOutput,
		"-outfmt", "%f", // fasta format
	)
	if output, err := queryCmd.CombinedOutput(); err != nil {
		return false, match, fmt.Errorf("failed to run blastdbcmd %s, %v, %s", entry, err, string(output))
	}

	// create blast input file
	inContent := fmt.Sprintf(">primer\n%s\n", primer)
	err = ioutil.WriteFile(in, []byte(inContent), 0666)
	if err != nil {
		return false, match, fmt.Errorf("failed to write primer sequence to query FASTA file: %v", err)
	}

	// blast the query sequence against the parent sequence
	b := blastExec{
		in:      in,
		out:     out,
		subject: entryOutput,
	}

	// execute blast
	if err = b.runAgainst(); err != nil {
		return false, match, fmt.Errorf("failed to run blast against parent: %v", err)
	}

	// get the BLAST matches
	matches, err := b.parse()
	if err != nil {
		return false, match, fmt.Errorf("failed to parse matches from %s: %v", out, err)
	}

	// parse the results and check whether any are cause for concern (by Tm)
	primerSeen := false
	for _, m := range matches {
		// one of the matches will, of course, be against the primer itself
		// and we don't want to double count it
		if !primerSeen && m.Seq == primer {
			primerSeen = true
			continue
		}

		if isMismatch(m) {
			return true, m, nil
		}
	}

	return false, match, nil
}