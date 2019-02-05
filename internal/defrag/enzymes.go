package defrag

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"os"
	"regexp"
	"sort"
	"strings"
	"text/tabwriter"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
)

// enzyme is a single enzyme that can be used to linearize a backbone before
// inserting a sequence
type enzyme struct {
	recog   string
	cutInd  int
	hangInd int
}

// parses a recognition sequence into a hangInd, cutInd for overhang calculation
func newEnzyme(recogSeq string) enzyme {
	cutIndex := strings.Index(recogSeq, "^")
	hangIndex := strings.Index(recogSeq, "_")

	if cutIndex < hangIndex {
		hangIndex--
	} else {
		cutIndex--
	}

	recogSeq = strings.Replace(recogSeq, "^", "", -1)
	recogSeq = strings.Replace(recogSeq, "_", "", -1)

	return enzyme{
		recog:   recogSeq,
		cutInd:  cutIndex,
		hangInd: hangIndex,
	}
}

// digest a Frag (backbone) with an enzyme's first recogition site
//
// remove the 5' end of the fragment post-cleaving. it will be degraded.
// keep exposed 3' ends. good visual explanation:
// https://warwick.ac.uk/study/csde/gsp/eportfolio/directory/pg/lsujcw/gibsonguide/
func digest(frag Frag, enz enzyme) (digested Frag, err error) {
	wrappedBp := 38 // largest current recognition site in the list of enzymes
	if len(frag.Seq) < wrappedBp {
		return Frag{}, fmt.Errorf("%s is too short for digestion", frag.ID)
	}

	if frag.fragType == circular {
		// double check that it's circular, first 20 bp should show up at least twice
		firstBpsRegex := regexp.MustCompile(frag.Seq[:20])
		firstBpsCount := firstBpsRegex.FindAllStringIndex(frag.Seq, -1)

		if len(firstBpsCount) > 1 {
			frag.Seq = frag.Seq[:len(frag.Seq)/2] // undo the doubling of sequence for circular parts
		}
	}

	// turn recognition site (with ambigous bps) into a recognition seq
	reg := regexp.MustCompile(recogRegex(enz.recog))
	seq := frag.Seq + frag.Seq[0:wrappedBp]
	revCompSeq := revComp(frag.Seq) + revComp(frag.Seq[0:wrappedBp])

	// positive if template strand has overhang
	// negative if rev comp strand has overhang
	templateOverhangLength := enz.cutInd - enz.hangInd
	cutIndex := -1
	digestedSeq := ""
	if reg.MatchString(seq) {
		// template
		cutIndex = reg.FindStringIndex(seq)[0] // first int is the start of match
	}
	if reg.MatchString(revCompSeq) {
		// reverse complement
		revCutIndex := reg.FindStringIndex(revCompSeq)[0]
		revCutIndex = len(frag.Seq) - revCutIndex - len(enz.recog) // flip it to account for being on rev comp
		revCutIndex = (revCutIndex + len(frag.Seq)) % len(frag.Seq)
		if revCutIndex >= 0 && (revCutIndex < cutIndex || cutIndex < 0) {
			cutIndex = revCutIndex // take whichever occurs sooner in the sequence
		}
	}
	if cutIndex == -1 {
		// no valid cutsites in the sequence
		return Frag{}, fmt.Errorf("no %s cutsites found in %s", enz.recog, frag.ID)
	}

	if templateOverhangLength >= 0 {
		cutIndex = (cutIndex + enz.cutInd) % len(frag.Seq)
		digestedSeq = frag.Seq[cutIndex:] + frag.Seq[:cutIndex]
	} else {
		bottomIndex := (cutIndex + enz.cutInd) % len(frag.Seq)
		topIndex := (cutIndex + enz.hangInd) % len(frag.Seq)
		digestedSeq = frag.Seq[topIndex:] + frag.Seq[:bottomIndex]
	}

	return Frag{
		ID:  frag.ID,
		Seq: digestedSeq,
	}, nil
}

// recogRegex turns a recognition sequence into a regex sequence for searching
// sequence for searching the template sequence for digestion sites
func recogRegex(recog string) (decoded string) {
	regexDecode := map[rune]string{
		'A': "A",
		'C': "C",
		'G': "G",
		'T': "T",
		'M': "(A|C)",
		'R': "(A|G)",
		'W': "(A|T)",
		'Y': "(C|T)",
		'S': "(C|G)",
		'K': "(G|T)",
		'H': "(A|C|T)",
		'D': "(A|G|T)",
		'V': "(A|C|G)",
		'B': "(C|G|T)",
		'N': "(A|C|G|T)",
		'X': "(A|C|G|T)",
	}

	var regexDecoder strings.Builder
	for _, c := range recog {
		regexDecoder.WriteString(regexDecode[c])
	}

	return regexDecoder.String()
}

// EnzymeDB is a struct for accessing defrags enzymes db
type EnzymeDB struct {
	// enzymes is a map between a enzymes name and its sequence
	enzymes map[string]string
}

// NewEnzymeDB returns a new copy of the enzymes db
func NewEnzymeDB() *EnzymeDB {
	enzymes := make(map[string]string)

	enzymeFile, err := os.Open(config.EnzymeDB)
	if err != nil {
		stderr.Fatal(err)
	}

	// https://golang.org/pkg/bufio/#example_Scanner_lines
	scanner := bufio.NewScanner(enzymeFile)
	for scanner.Scan() {
		columns := strings.Split(scanner.Text(), "	")
		enzymes[columns[0]] = columns[1] // enzyme name = enzyme seq
	}

	if err := enzymeFile.Close(); err != nil {
		stderr.Fatal(err)
	}

	return &EnzymeDB{enzymes: enzymes}
}

// ReadCmd returns enzymes that are similar in name to the enzyme name requested.
// if multiple enzyme names include the enzyme name, they are all returned.
// otherwise a list of enzyme names are returned (those beneath a levenshtein distance cutoff)
func (f *EnzymeDB) ReadCmd(cmd *cobra.Command, args []string) {
	// from https://golang.org/pkg/text/tabwriter/
	w := tabwriter.NewWriter(os.Stdout, 0, 0, 3, ' ', tabwriter.TabIndent)

	if len(args) < 1 {
		enzymeNames := make([]string, len(f.enzymes), len(f.enzymes))
		i := 0
		for name := range f.enzymes {
			enzymeNames[i] = name
			i++
		}
		sort.Strings(enzymeNames)

		for _, name := range enzymeNames {
			fmt.Fprintf(w, "%s\t%s\n", name, f.enzymes[name])
		}
		w.Flush()
		return
	}

	name := args[0]

	// if there's an exact match, just long that one
	if seq, exists := f.enzymes[name]; exists {
		fmt.Printf("%s	%s\n", name, seq)
		return
	}

	ldCutoff := 2
	containing := []string{}
	lowDistance := []string{}

	for fName, fSeq := range f.enzymes {
		if strings.Contains(fName, name) {
			containing = append(containing, fName+"\t"+fSeq)
		} else if len(fName) > ldCutoff && ld(name, fName, true) <= ldCutoff {
			lowDistance = append(lowDistance, fName+"\t"+fSeq)
		}
	}

	if len(containing) < 3 {
		lowDistance = append(lowDistance, containing...)
		containing = []string{} // clear
	}
	if len(containing) > 0 {
		fmt.Fprintf(w, strings.Join(containing, "\n"))
	} else if len(lowDistance) > 0 {
		fmt.Fprintf(w, strings.Join(lowDistance, "\n"))
	} else {
		fmt.Fprintf(w, fmt.Sprintf("failed to find any enzymes for %s", name))
	}
	w.Write([]byte("\n"))
	w.Flush()
}

// SetCmd the enzyme's seq in the database (or create if it isn't in the enzyme db)
func (f *EnzymeDB) SetCmd(cmd *cobra.Command, args []string) {
	if len(args) < 2 {
		stderr.Fatalf("expecting two args: a name and recognition sequence. see 'defrag find enzyme --help'\n")
	}

	name := args[0]
	seq := args[1]
	if len(args) > 2 {
		name = strings.Join(args[:len(args)-1], " ")
		seq = args[len(args)-1]
	}
	seq = strings.ToUpper(seq)

	invalidChars := regexp.MustCompile("[^ATGCMRWYSKHDVBNX_\\^]")
	seq = invalidChars.ReplaceAllString(seq, "")

	if strings.Count(seq, "^") != 1 || strings.Count(seq, "_") != 1 {
		stderr.Fatalf("%s is not a valid enzyme recognition sequence. see 'defrag find enzyme --help'\n", seq)
	}

	enzymeFile, err := os.Open(config.EnzymeDB)
	if err != nil {
		stderr.Fatal(err)
	}

	// https://golang.org/pkg/bufio/#example_Scanner_lines
	var output strings.Builder
	updated := false
	scanner := bufio.NewScanner(enzymeFile)
	for scanner.Scan() {
		columns := strings.Split(scanner.Text(), "	")
		if columns[0] == name {
			output.WriteString(fmt.Sprintf("%s	%s\n", name, seq))
			updated = true
		} else {
			output.WriteString(scanner.Text())
		}
	}

	// create from nothing
	if !updated {
		output.WriteString(fmt.Sprintf("%s	%s\n", name, seq))
	}

	if err := enzymeFile.Close(); err != nil {
		stderr.Fatal(err)
	}

	if err := ioutil.WriteFile(config.EnzymeDB, []byte(output.String()), 0644); err != nil {
		stderr.Fatal(err)
	}

	if updated {
		fmt.Printf("updated %s in the enzymes database\n", name)
	}

	// update in memory
	f.enzymes[name] = seq
}

// DeleteCmd the enzyme from the database
func (f *EnzymeDB) DeleteCmd(cmd *cobra.Command, args []string) {
	if len(args) < 1 {
		stderr.Fatalf("expecting one arg: a enzymes name. %d passed\n", len(args))
	}

	name := args[0]
	if len(args) > 1 {
		name = strings.Join(args, " ")
	}

	if _, contained := f.enzymes[name]; !contained {
		fmt.Printf("failed to find %s in the enzymes database\n", name)
	}

	enzymeFile, err := os.Open(config.EnzymeDB)
	if err != nil {
		stderr.Fatal(err)
	}

	// https://golang.org/pkg/bufio/#example_Scanner_lines
	var output strings.Builder
	deleted := false
	scanner := bufio.NewScanner(enzymeFile)
	for scanner.Scan() {
		columns := strings.Split(scanner.Text(), "	")
		if columns[0] != name {
			output.WriteString(scanner.Text())
		} else {
			deleted = true
		}
	}

	if err := enzymeFile.Close(); err != nil {
		stderr.Fatal(err)
	}

	if err := ioutil.WriteFile(config.EnzymeDB, []byte(output.String()), 0644); err != nil {
		stderr.Fatal(err)
	}

	// delete from memory
	delete(f.enzymes, name)

	if deleted {
		fmt.Printf("deleted %s from the enzymes database\n", name)
	} else {
		fmt.Printf("failed to find %s in the enzymes database\n", name)
	}
}
