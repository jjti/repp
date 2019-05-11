package rvec

import (
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"text/tabwriter"

	"github.com/jjtimmons/rvec/config"
)

var (
	// blastnDir is a temporary directory for all blastn output
	blastnDir = ""

	// blastdbcmd is a temporary directory for all blastdbcmd output
	blastdbcmdDir = ""
)

// match is a blast "hit" in the blastdb.
type match struct {
	// entry of the matched building fragment in the database
	entry string

	// unique id for the match (entry name + start index % seqL). also used by fragments
	uniqueID string

	// seq of the match on the subject match
	querySeq string

	// queryStart of the queried seq match (0-indexed)
	queryStart int

	// queryEnd of the queried seq match (0-indexed)
	queryEnd int

	// sequence of the match in the subject sequence
	seq string

	// start index of the match on the subject fragment (entry) (0-indexed)
	subjectStart int

	// end index of the match on the subject fragment (entry) (0-indexed)
	subjectEnd int

	// the db that was BLASTed against (used later for checking off-targets in parents)
	db string

	// titles from the db. eg: year it was created
	title string

	// circular if it's a circular fragment in the db (vector, plasmid, etc)
	circular bool

	// mismatching number of bps in the match (for primer off-targets)
	mismatching int

	// internal if the fragment doesn't have to be procured from a remote repository (eg Addgene, iGEM)
	internal bool

	// forward if the match is along the sequence strand versus the reverse complement strand
	forward bool
}

// blastExec is a small utility object for executing BLAST.
type blastExec struct {
	// the name of the query
	name string

	// the sequence of the query
	seq string

	// whether to circularize the queries sequence in the input file
	circular bool

	// the path to the database we're BLASTing against
	db string

	// the input BLAST file
	in *os.File

	// the output BLAST file
	out *os.File

	// optional path to a FASTA file with a subject FASTA sequence
	subject string

	// internal if the db is a local/user owned list of fragments (ie free)
	internal bool

	// the percentage identity for BLAST queries
	identity int

	// the expect value of a BLAST query (defaults to 10)
	evalue int
}

// length returns the length of the match on the queried fragment.
func (m *match) length() int {
	queryLength := m.queryEnd - m.queryStart + 1
	subjectLength := m.subjectEnd - m.subjectStart + 1

	if queryLength > subjectLength {
		return queryLength
	}

	return subjectLength
}

// copyWithQueryRange returns a new match with the new start, end.
func (m *match) copyWithQueryRange(start, end int) match {
	return match{
		entry:        m.entry,
		uniqueID:     m.uniqueID,
		seq:          m.seq,
		queryStart:   start,
		queryEnd:     end,
		subjectStart: m.subjectStart,
		subjectEnd:   m.subjectEnd,
		db:           m.db,
		circular:     m.circular,
		mismatching:  m.mismatching,
		internal:     m.internal,
	}
}

// log the match in string format
func (m *match) log() {
	fmt.Printf("%s %d %d\n", m.entry, m.queryStart, m.queryEnd)
}

// blast the seq against all dbs and acculate matches.
func blast(
	name, seq string,
	circular bool,
	dbs, filters []string,
	identity int,
	tw *tabwriter.Writer,
) (matches []match, err error) {
	in, err := ioutil.TempFile(blastnDir, name+"in-*")
	if err != nil {
		return nil, err
	}
	defer os.Remove(in.Name())

	out, err := ioutil.TempFile(blastnDir, name+"out-*")
	if err != nil {
		return nil, err
	}
	defer os.Remove(out.Name())

	for _, db := range dbs {
		internal := true
		if strings.Contains(db, "addgene") || strings.Contains(db, "igem") {
			internal = false
		}

		b := &blastExec{
			name:     name,
			seq:      seq,
			circular: circular,
			db:       db,
			in:       in,
			out:      out,
			internal: internal,
			identity: identity,
		}

		// make sure the db exists
		if _, err := os.Stat(db); os.IsNotExist(err) {
			return nil, fmt.Errorf("failed to find a BLAST database at %s", db)
		}

		// create the input file
		if err := b.input(); err != nil {
			return nil, fmt.Errorf("failed to write a BLAST input file at %s: %v", b.in.Name(), err)
		}

		// execute BLAST
		if err := b.run(); err != nil {
			return nil, fmt.Errorf("failed executing BLAST: %v", err)
		}

		// parse the output file to Matches against the Frag
		dbMatches, err := b.parse(filters)
		if err != nil {
			return nil, fmt.Errorf("failed to parse BLAST output: %v", err)
		}

		fmt.Fprintf(tw, "%s\t%d\t%s\n", name, len(dbMatches), db)

		// add these matches against the growing list of matches
		matches = append(matches, dbMatches...)
	}

	return matches, nil
}

// blastWriter returns a new tabwriter specifically for blast database calls.
func blastWriter() *tabwriter.Writer {
	tw := tabwriter.NewWriter(os.Stdout, 0, 4, 3, ' ', 0)
	fmt.Fprintf(tw, "entry\tmatches\tdatabase\t\n")

	return tw
}

// input creates an input query file (FASTA) for blastn.
func (b *blastExec) input() error {
	// create the query sequence file.
	// if circular, add the sequence to itself because it's circular
	// and we want to find matches across the zero-index
	querySeq := b.seq
	if b.circular {
		querySeq = querySeq + b.seq
	}

	_, err := b.in.WriteString(fmt.Sprintf(">%s\n%s\n", b.name, querySeq))

	return err
}

// run calls the external blastn binary on the input file.
func (b *blastExec) run() (err error) {
	threads := runtime.NumCPU() - 1
	if threads < 1 {
		threads = 1
	}

	flags := []string{
		"-task", "blastn",
		"-db", b.db,
		"-query", b.in.Name(),
		"-out", b.out.Name(),
		"-outfmt", "7 sseqid qstart qend sstart send sseq mismatch gaps stitle",
		"-perc_identity", fmt.Sprintf("%d", b.identity),
		"-culling_limit", "50", // "If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit"
		"-num_threads", strconv.Itoa(threads),
	}

	if b.identity > 99 {
		flags = append(flags,
			"-reward", "1", // most negative penalty I could find in blast_stat.c
			"-penalty", "-5", // needed because mismatches were being included in the end of pSB1A3 matches
			"-gapopen", "6",
			"-gapextend", "6",
		)
	} else if b.identity >= 98 {
		// Table D1 https://www.ncbi.nlm.nih.gov/books/NBK279684/
		flags = append(flags,
			"-reward", "1",
			"-penalty", "-3",
			"-gapopen", "3",
			"-gapextend", "3",
		)
	} else if b.identity >= 90 {
		// Table D1 https://www.ncbi.nlm.nih.gov/books/NBK279684/
		flags = append(flags,
			"-reward", "1",
			"-penalty", "-2",
			"-gapopen", "1",
			"-gapextend", "2",
		)
	} else {
		// Table D1 https://www.ncbi.nlm.nih.gov/books/NBK279684/
		flags = append(flags,
			"-reward", "1",
			"-penalty", "-1",
			"-gapopen", "1",
			"-gapextend", "2",
		)
	}

	if b.evalue != 0 {
		flags = append(flags, "-evalue", strconv.Itoa(b.evalue))
	} else if b.identity < 90 {
		flags = append(flags, "-evalue", "1000")
	} else if b.identity < 98 {
		flags = append(flags, "-evalue", "500")
	}

	// https://www.ncbi.nlm.nih.gov/books/NBK279682/
	blastCmd := exec.Command("blastn", flags...)

	// execute BLAST and wait on it to finish
	if output, err := blastCmd.CombinedOutput(); err != nil {
		return fmt.Errorf("failed to execute blastn against %s: %v: %s", b.db, err, string(output))
	}

	return
}

// parse reads the output of blastn into matches.
func (b *blastExec) parse(filters []string) (matches []match, err error) {
	// read in the results
	file, err := ioutil.ReadFile(b.out.Name())
	if err != nil {
		return
	}
	fileS := string(file)

	fullQuery := b.seq + b.seq
	identityThreshold := float32(b.identity)/100.0 - 0.0001

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

		entry := strings.Replace(cols[0], ">", "", -1)
		queryStart, _ := strconv.Atoi(cols[1])
		queryEnd, _ := strconv.Atoi(cols[2])
		subjectStart, _ := strconv.Atoi(cols[3])
		subjectEnd, _ := strconv.Atoi(cols[4])
		seq := cols[5]                          // subject sequence
		mismatching, _ := strconv.Atoi(cols[6]) // mismatch count
		gaps, _ := strconv.Atoi(cols[7])        // gap count
		titles := cols[8]                       // salltitles, eg: "fwd-terminator-2011"
		forward := true

		// check whether the mismatch ratio is less than the set limit
		matchRatio := float32(len(seq)-(mismatching+gaps)) / float32(len(seq))
		if matchRatio < identityThreshold {
			// fmt.Printf("getting rid of %s", entry)
			continue
		}

		seq = strings.Replace(seq, "-", "", -1) // remove gap markers
		queryStart--                            // convert from 1-based to 0-based
		queryEnd--
		subjectStart--
		subjectEnd--

		// bug where titles are being included in the entry
		entryCols := strings.Fields(entry)
		if len(entryCols) > 1 {
			entry = entryCols[0]
			titles = entryCols[1] + titles
		}

		// flip if blast is reading right to left
		if queryStart > queryEnd {
			queryStart, queryEnd = queryEnd, queryStart
			forward = !forward
		}
		if subjectStart > subjectEnd {
			subjectStart, subjectEnd = subjectEnd, subjectStart
			forward = !forward
		}

		// filter on titles
		matchesFilter := false
		titles += entry
		titles = strings.ToUpper(titles)
		for _, f := range filters {
			if strings.Contains(titles, f) {
				matchesFilter = true
				break
			}
		}
		if matchesFilter {
			continue // has been filtered out because of the "exclude" CLI flag
		}

		// get a unique identifier to distinguish this match/fragment from the others
		uniqueID := entry + strconv.Itoa(queryStart%len(b.seq))

		// gather the query sequence
		querySeq := fullQuery[queryStart : queryEnd+1]

		// create and append the new match
		ms = append(ms, match{
			entry:        entry,
			uniqueID:     uniqueID,
			querySeq:     querySeq,
			queryStart:   queryStart,
			queryEnd:     queryEnd,
			seq:          seq,
			subjectStart: subjectStart,
			subjectEnd:   subjectEnd,
			circular:     strings.Contains(entry+titles, "circular"),
			mismatching:  mismatching + gaps,
			internal:     b.internal,
			db:           b.db,
			title:        titles,
			forward:      forward,
		})
	}

	return ms, nil
}

// culling removes matches that are engulfed in others
//
// culling fragment matches means removing those that are completely
// self-contained in other fragments: the larger of the available fragments
// will be the better one, since it covers a greater region and will almost
// always be preferable to the smaller one
//
// also remove small fragments here that are too small to be useful during assembly
func cull(matches []match, targetLength, minSize int) (culled []match) {
	culled = []match{}

	// remove fragments that are shorter the minimum cut off size
	// separate the internal and external fragments. the internal
	// ones should not be removed just if they're self-contained
	// in another, because they may be cheaper to assemble
	var internal []match
	var external []match
	for _, m := range matches {
		if minSize > 0 && m.length() < minSize {
			continue // too short
		}

		if m.internal {
			internal = append(internal, m)
		} else {
			external = append(external, m)
		}
	}

	// create culled matches (non-self contained)
	culled = append(properize(internal), properize(external)...)

	// because we culled the matches, we may have removed a match from the
	// start or the end. right now, a match showing up twice in the vector
	// is how we circularize, so have to add back matches to the start or end
	matchCount := make(map[string]int)
	for _, m := range culled {
		if _, counted := matchCount[m.uniqueID]; counted {
			matchCount[m.uniqueID]++
		} else {
			matchCount[m.uniqueID] = 1
		}
	}

	// sort again now that we added copied matches
	sortMatches(culled)

	// fmt.Println("matches")
	// for _, m := range matches {
	// 	if m.queryStart < targetLength {
	// 		m.log()
	// 	}
	// }

	// fmt.Println("after culling")
	// for _, m := range culled {
	// 	if m.queryStart < targetLength {
	// 		m.log()
	// 	}
	// }

	return culled
}

// properize remove matches that are entirely contained within others
func properize(matches []match) (culled []match) {
	sortMatches(matches)

	// only include those that aren't encompassed by the one before it
	for _, m := range matches {
		lastMatch := len(culled) - 1
		if lastMatch < 0 || m.queryEnd > culled[lastMatch].queryEnd {
			culled = append(culled, m)
		}
	}

	return
}

// sortMatches sorts matches by their start index
// for fragments with equivelant starting indexes, put the larger one first
func sortMatches(matches []match) {
	sort.Slice(matches, func(i, j int) bool {
		if matches[i].queryStart != matches[j].queryStart {
			return matches[i].queryStart < matches[j].queryStart
		} else if matches[i].length() != matches[j].length() {
			return matches[i].length() > matches[j].length()
		}
		return matches[i].entry > matches[j].entry
	})
}

// queryDatabases is for finding a fragment/vector with the entry name in one of the dbs
func queryDatabases(entry string, dbs []string) (f *Frag, err error) {
	// first try to get the entry out of a local file
	if frags, err := read(entry, false); err == nil && len(frags) > 0 {
		return frags[0], nil // it was a local file
	}

	// channel that returns filename to an output result from blastdbcmd
	outFileCh := make(chan string, len(dbs))
	dbSourceCh := make(chan string, len(dbs))

	// move through each db and see if it contains the entry
	for _, db := range dbs {
		go func(db string) {
			// if outFile is defined here we managed to query the entry from the db
			outFile, err := blastdbcmd(entry, db)
			if err == nil && outFile != nil {
				outFileCh <- outFile.Name() // "" if not found
				dbSourceCh <- db
			} else {
				outFileCh <- ""
				dbSourceCh <- ""
			}
		}(db)
	}

	// try and return a read fragment file
	for i := 0; i < len(dbs); i++ {
		outFile := <-outFileCh
		dbSource := <-dbSourceCh
		if outFile == "" {
			continue // failed to query from this DB
		}

		defer os.Remove(outFile)

		if frags, err := read(outFile, false); err == nil {
			targetFrag := frags[0]

			// fix the ID, don't want titles in the ID (bug)
			idSplit := strings.Fields(targetFrag.ID)
			if len(idSplit) > 1 {
				targetFrag.ID = idSplit[0]
			}

			targetFrag.db = dbSource
			return targetFrag, nil
		}

		return &Frag{}, err
	}

	close(outFileCh)
	close(dbSourceCh)

	sep := "\n\t"

	return &Frag{}, fmt.Errorf("failed to find %s in any of:%s", entry, sep+strings.Join(dbs, sep))
}

// seqMismatch queries for any mismatching primer locations in the parent sequence
// unlike parentMismatch, it doesn't first find the parent fragment from the db it came from
// the sequence is passed directly as parentSeq
func seqMismatch(primers []Primer, parentID, parentSeq string, conf *config.Config) (wasMismatch bool, m match, err error) {
	parentFile, err := ioutil.TempFile(blastnDir, "parent-*")
	if err != nil {
		return false, match{}, err
	}
	defer os.Remove(parentFile.Name())

	if parentID == "" {
		parentID = "parent"
	}
	inContent := fmt.Sprintf(">%s\n%s\n", parentID, parentSeq)
	if _, err = parentFile.WriteString(inContent); err != nil {
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

// parentMismatch both searches for a the parent fragment in its source DB and queries for
// any mismatches in the seq before returning
func parentMismatch(primers []Primer, parent, db string, conf *config.Config) (wasMismatch bool, m match, err error) {
	// try and query for the parent in the source DB and write to a file
	parentFile, err := blastdbcmd(parent, db)

	// ugly check here for whether we just failed to get the parent entry from a db
	// which isn't a huge deal (shouldn't be flagged as a mismatch)
	// this is similar to what io.IsNotExist does
	if err != nil {
		if strings.Contains(err.Error(), "failed to query") {
			stderr.Println(err) // just write the error
			// TODO: if we fail to find the parent, query the fullSeq as it was sent
			return false, match{}, nil
		}
		return false, match{}, err
	}

	// check each primer for mismatches
	if parentFile.Name() != "" {
		defer os.Remove(parentFile.Name())

		for _, primer := range primers {
			wasMismatch, m, err = mismatch(primer.Seq, parentFile, conf)
			if wasMismatch || err != nil {
				return
			}
		}
	}

	return
}

// blastdbcmd queries a fragment/vector by its FASTA entry name (entry) and writes the
// results to a temporary file (to be BLAST'ed against)
//
// entry here is the ID that's associated with the fragment in its source DB (db)
func blastdbcmd(entry, db string) (output *os.File, err error) {
	// path to the entry batch file to hold the entry accession
	entryFile, err := ioutil.TempFile(blastdbcmdDir, "in-*")
	if err != nil {
		return nil, err
	}
	defer os.Remove(entryFile.Name())

	// path to the output sequence file from querying the entry's sequence from the BLAST db
	output, err = ioutil.TempFile(blastdbcmdDir, "out-*")
	if err != nil {
		return nil, err
	}

	// write entry to file
	// this was a 2-day issue I couldn't resolve...
	// I was using the "-entry" flag on exec.Command, but have since
	// switched to the simpler -entry_batch command (on a file) that resolves the issue
	if _, err := entryFile.WriteString(entry); err != nil {
		return nil, fmt.Errorf("failed to write blastdbcmd entry file at %s: %v", entryFile.Name(), err)
	}

	// make a blastdbcmd command (for querying a DB, very different from blastn)
	queryCmd := exec.Command(
		"blastdbcmd",
		"-db", db,
		"-dbtype", "nucl",
		"-entry_batch", entryFile.Name(),
		"-out", output.Name(),
		"-outfmt", "%f ", // fasta format
	)

	// execute
	if _, err := queryCmd.CombinedOutput(); err != nil {
		return nil, fmt.Errorf("warning: failed to query %s from %s\n\t%s", entry, db, err.Error())
	}

	// read in the results as fragments. set their sequence to the full one returned from blastdbcmd
	fragments, err := read(output.Name(), false)
	if err == nil && len(fragments) >= 1 {
		for _, f := range fragments {
			f.fullSeq = f.Seq // set fullSeq, faster to check for primer off-targets later
		}
		return output, nil
	}

	return nil, fmt.Errorf("warning: failed to query %s from %s", entry, db)
}

// mismatch finds mismatching sequences between the query sequence and
// the parent sequence (in the parent file)
//
// The fragment to query against is stored in parentFile
func mismatch(primer string, parentFile *os.File, c *config.Config) (wasMismatch bool, m match, err error) {
	// path to the entry batch file to hold the entry accession
	in, err := ioutil.TempFile(blastnDir, "primer.in-*")
	if err != nil {
		return false, match{}, err
	}
	defer os.Remove(in.Name())

	// path to the output sequence file from querying the entry's sequence from the BLAST db
	out, err := ioutil.TempFile(blastnDir, "primer.out-*")
	if err != nil {
		return false, match{}, err
	}
	defer os.Remove(out.Name())

	// create input file
	inContent := fmt.Sprintf(">primer\n%s\n", primer)
	if _, err = in.WriteString(inContent); err != nil {
		return false, m, fmt.Errorf("failed to write primer sequence to query FASTA file: %v", err)
	}

	// BLAST the query sequence against the parentFile sequence
	b := blastExec{
		in:       in,
		out:      out,
		subject:  parentFile.Name(),
		seq:      primer,
		identity: 65,    // see Primer-BLAST https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3412702/
		evalue:   30000, // see Primer-BLAST
	}

	// execute BLAST
	if err = b.runAgainst(); err != nil {
		return false, m, fmt.Errorf("failed to run blast against parent: %v", err)
	}

	// get the BLAST matches
	matches, err := b.parse([]string{})
	if err != nil {
		return false, match{}, fmt.Errorf("failed to parse matches from %s: %v", out.Name(), err)
	}

	// parse the results and check whether any are cause for concern (by Tm)
	primerCount := 1 // number of times we expect to see the primer itself
	parentFileContents, err := ioutil.ReadFile(parentFile.Name())
	if err != nil {
		return false, match{}, err
	}

	if strings.Contains(string(parentFileContents), "circular") {
		// if the match is against a circular fragment, we expect to see the primer's binding location
		// twice because circular fragments' sequences are doubled in the DBs
		primerCount++
	}

	for _, m := range matches {
		if isMismatch(primer, m, c) {
			primerCount--
		}

		if primerCount < 0 {
			return true, m, nil
		}
	}

	return false, match{}, nil
}

// runs blast on the query file against another subject file (rather than blastdb)
func (b *blastExec) runAgainst() (err error) {
	// create the blast command
	// https://www.ncbi.nlm.nih.gov/books/NBK279682/
	blastCmd := exec.Command(
		"blastn",
		"-task", "blastn",
		"-query", b.in.Name(),
		"-subject", b.subject,
		"-out", b.out.Name(),
		"-outfmt", "7 sseqid qstart qend sstart send sseq mismatch gaps stitle",
	)

	// execute BLAST and wait on it to finish
	if output, err := blastCmd.CombinedOutput(); err != nil {
		return fmt.Errorf("failed to execute blastn against %s: %v: %s", b.subject, err, string(output))
	}

	return
}

// isMismatch returns whether the match constitutes a mismatch
// between it and the would be primer sequence
//
// estimate the ntthal and check against the max offtarget tm
// from the settings
func isMismatch(primer string, m match, c *config.Config) bool {
	// we want the reverse complement of one to the other
	ectopic := m.seq
	if m.forward {
		ectopic = reverseComplement(ectopic)
	}

	ntthalCmd := exec.Command(
		"ntthal",
		"-a", "END1", // end of primer sequence
		"-r", // temperature only
		"-s1", primer,
		"-s2", ectopic,
		"-path", config.Primer3Config,
	)

	ntthalOut, err := ntthalCmd.CombinedOutput()
	if err != nil {
		stderr.Fatal(err)
	}

	ntthalOutString := string(ntthalOut)
	temp, err := strconv.ParseFloat(strings.TrimSpace(ntthalOutString), 64)
	if err != nil {
		stderr.Fatalln(err)
	}

	return temp > c.PCRMaxOfftargetTm
}

func init() {
	var err error

	blastnDir, err = ioutil.TempDir("", "blastn")
	if err != nil {
		stderr.Fatal(err)
	}

	blastdbcmdDir, err = ioutil.TempDir("", "blastdbcmd")
	if err != nil {
		stderr.Fatal(err)
	}
}
