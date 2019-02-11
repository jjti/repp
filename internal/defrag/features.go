package defrag

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"os"
	"sort"
	"strconv"
	"strings"
	"text/tabwriter"
	"time"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
)

// FeatureDB is a struct for accessing defrags features db
type FeatureDB struct {
	// features is a map between a features name and its sequence
	features map[string]string
}

type featureMatch struct {
	featureIndex int
	match        match
}

// FeaturesCmd accepts a cobra commands and assembles a vector containing all the features
func FeaturesCmd(cmd *cobra.Command, args []string) {
	features(parseCmdFlags(cmd, args))
}

// features assembles a vector with all the features requested with the 'defrag features [feature ...]' command
// defrag assemble features p10 promoter, mEGFP, T7 terminator
func features(flags *Flags, conf *config.Config) {
	start := time.Now()

	// turn feature names into sequences
	insertFeatures, backboneFeature := queryFeatures(flags)
	targetFeatures := append(insertFeatures, backboneFeature)

	// find matches in the databases
	fragToFeatureMatches := blastFeatures(flags, targetFeatures)

	// build assemblies containing the matched fragments
	syntheticVector, solutions := buildFeatureSolutions(targetFeatures, fragToFeatureMatches, flags.dbs, conf)

	// write the output file
	insertLength := 0
	for _, f := range targetFeatures {
		insertLength += len(f[1])
	}
	write(flags.out, flags.in, syntheticVector, solutions, insertLength, conf, time.Since(start).Seconds())
}

// queryFeatures takes the list of feature names and finds them in the available databases
func queryFeatures(flags *Flags) ([][]string, []string) {
	insertFeatures := [][]string{} // slice of tuples [feature name, feature sequence]

	if readFeatures, err := read(flags.in, true); err == nil {
		// see if the features are in a file (multi-FASTA or features in a Genbank)
		seenFeatures := make(map[string]string) // map feature name to sequence
		for _, f := range readFeatures {
			if seq := seenFeatures[f.ID]; seq != f.Seq {
				stderr.Fatalf("failed to parse features, %s has two different sequences:\n\t%s\n\t%s\n", f.ID, f.Seq, seq)
			}
			insertFeatures = append(insertFeatures, []string{f.ID, f.Seq})
		}
	} else {
		// if the features weren't in a file, try and find each in the features database
		// or one of the databases passed as a source of building fragments
		featureNames := []string{}
		if strings.Contains(flags.in, ",") {
			featureNames = strings.Split(flags.in, ",") // comma separated
		} else {
			featureNames = strings.Fields(flags.in) // spaces
		}

		if len(featureNames) < 1 {
			stderr.Fatal("no features chosen. see 'defrag features --help'")
		}

		featureDB := NewFeatureDB()
		for _, f := range featureNames {
			if seq, contained := featureDB.features[f]; contained {
				insertFeatures = append(insertFeatures, []string{f, seq})
			} else if dbFrag, err := queryDatabases(f, flags.dbs); err == nil {
				insertFeatures = append(insertFeatures, []string{f, dbFrag.Seq})
			} else {
				stderr.Fatalf(
					"failed to find '%s' in the features database (%s) or any of:"+
						"%s\ncheck features database with 'defrag features find [feature name]'",
					f,
					config.FeatureDB,
					"\n  "+strings.Join(flags.dbs, "\n  "),
				)
			}
		}
	}

	// add in the backbone as a "feature"
	backboneFeature := []string{}
	if flags.backbone.ID != "" {
		backboneFeature = []string{flags.backbone.ID, flags.backbone.Seq}
	}

	return insertFeatures, backboneFeature
}

// blastFeatures returns matches between the target features and entries in the databases with those features
func blastFeatures(flags *Flags, targetFeatures [][]string) map[string][]featureMatch {
	fragToFeatureMatches := make(map[string][]featureMatch) // a list from each entry (by id) to its list of matched features
	tw := blastWriter()

	for i, target := range targetFeatures {
		matches, err := blast(target[0], target[1], false, flags.dbs, flags.filters, flags.identity, tw)
		if err != nil {
			stderr.Fatalln(err)
		}

		for _, m := range matches {
			m.queryStart = i
			m.queryEnd = i
			m.uniqueID = m.entry + strconv.Itoa(m.subjectStart)

			if _, exists := fragToFeatureMatches[m.entry]; !exists {
				fragToFeatureMatches[m.entry] = []featureMatch{featureMatch{featureIndex: i, match: m}}
			} else {
				fragToFeatureMatches[m.entry] = append(fragToFeatureMatches[m.entry], featureMatch{featureIndex: i, match: m})
			}
		}
	}

	tw.Flush()

	return fragToFeatureMatches
}

// buildFeatureSolultions creates and fills the assemblies using the matched fragments
func buildFeatureSolutions(
	targetFeatures [][]string,
	fragToFeatureMatches map[string][]featureMatch,
	dbs []string,
	conf *config.Config) (string, [][]*Frag) {
	accMatches := []match{}

	// turn matches into frags but with query start and end referring to feature index
	for _, matches := range fragToFeatureMatches {
		sort.Slice(matches, func(i, j int) bool {
			return matches[i].featureIndex < matches[j].featureIndex
		})

		// expand the range the matches range based on its continuous feature stretches
		m := matches[0].match
		stretchStart := matches[0].featureIndex

		// features doubled on self to account for feature runs that cross the zero index
		for mIndex, featureMatch := range append(matches, matches...) {
			if featureMatch.featureIndex == (stretchStart+1)%len(targetFeatures) && mIndex < len(matches)*2-1 {
				continue // continue the feature match stretch
			}

			// cannot extend this stretch, create a new fragment
			m.queryStart = stretchStart
			m.queryEnd = featureMatch.featureIndex
			m.subjectEnd = matches[mIndex%len(matches)].match.subjectEnd

			accMatches = append(accMatches, featureMatch.match)

			// start on the next stretch
			m = matches[(mIndex+1)%len(matches)].match
			stretchStart = matches[(mIndex+1)%len(matches)].featureIndex
		}
	}

	// filter out matches that are completely contained in others
	accMatches = filter(accMatches, len(targetFeatures), conf.PCRMinLength)

	fmt.Printf("%d matches after filtering\n", len(accMatches))

	// get the full vector length if just synthesizing each feature next to one another
	var fullSynthSeqBuilder strings.Builder
	featureToStart := make(map[int]int) // map from feature index to start index on fullSynthSeq
	for i, feat := range targetFeatures {
		featureToStart[i] = fullSynthSeqBuilder.Len()
		fullSynthSeqBuilder.WriteString(feat[1])
	}
	syntheticVector := fullSynthSeqBuilder.String()

	// get the matches back out of the databases (the full parts)
	frags := []*Frag{}
	for _, m := range accMatches {
		frag, err := queryDatabases(m.entry, dbs)
		if err != nil {
			stderr.Fatalln(err)
		}

		start := m.subjectStart
		end := m.subjectEnd
		if end < start {
			end += len(frag.Seq)
		}

		frag.ID = m.entry
		frag.uniqueID = m.uniqueID
		frag.Seq = (frag.Seq + frag.Seq)[start : end+1]
		frag.conf = conf

		frag.start = featureToStart[m.queryStart]
		frag.end = featureToStart[(m.queryEnd+1)%len(targetFeatures)]

		if frag.end <= frag.start {
			frag.end += len(frag.Seq)
		}

		if m.queryStart >= len(targetFeatures) {
			frag.start += len(syntheticVector)
			frag.end += len(syntheticVector)
		}

		frags = append(frags, &frag)
	}

	// traverse the fragments, accumulate assemblies that span all the features
	assemblies := createAssemblies(frags, len(syntheticVector), conf)

	// build up a map from fragment count to a sorted list of assemblies with that number
	assemblyCounts, countToAssemblies := groupAssembliesByCount(assemblies)

	// fill each assembly and accumulate the pareto optimal solutions
	solutions := fillSolutions(syntheticVector, assemblyCounts, countToAssemblies, conf)

	return syntheticVector, solutions
}

// NewFeatureDB returns a new copy of the features db
func NewFeatureDB() *FeatureDB {
	features := make(map[string]string)

	featureFile, err := os.Open(config.FeatureDB)
	if err != nil {
		stderr.Fatal(err)
	}

	// https://golang.org/pkg/bufio/#example_Scanner_lines
	scanner := bufio.NewScanner(featureFile)
	for scanner.Scan() {
		columns := strings.Split(scanner.Text(), "	")
		features[columns[0]] = columns[1] // feature name = feature seq
	}

	if err := featureFile.Close(); err != nil {
		stderr.Fatal(err)
	}

	return &FeatureDB{features: features}
}

// ReadCmd returns features that are similar in name to the feature name requested.
// if multiple feature names include the feature name, they are all returned.
// otherwise a list of feature names are returned (those beneath a levenshtein distance cutoff)
func (f *FeatureDB) ReadCmd(cmd *cobra.Command, args []string) {
	w := tabwriter.NewWriter(os.Stdout, 0, 0, 3, ' ', tabwriter.TabIndent)

	if len(args) < 1 {
		// no feature name passed, log all of them
		featNames := []string{}
		for feat := range f.features {
			featNames = append(featNames, feat)
		}
		sort.Slice(
			featNames,
			func(i, j int) bool {
				return strings.ToLower(featNames[i]) < strings.ToLower(featNames[j])
			},
		)

		// print all their names to the console and the first few bp
		for _, feat := range featNames {
			seq := f.features[feat]
			if len(seq) > 20 {
				seq = seq[:20] + "..."
			}
			fmt.Fprintf(w, "%s\t%s\n", feat, seq)
		}

		w.Flush()
		return
	}

	name := args[0]
	if len(args) > 1 {
		name = strings.Join(args, " ")
	}

	ldCutoff := len(name) / 3
	if 2 > ldCutoff {
		ldCutoff = 2
	}
	containing := []string{}
	lowDistance := []string{}

	for fName, fSeq := range f.features {
		if strings.Contains(fName, name) {
			containing = append(containing, fName+"\t"+fSeq)
		} else if len(fName) > ldCutoff && ld(name, fName, true) <= ldCutoff {
			lowDistance = append(lowDistance, fName+"\t"+fSeq)
		}
	}

	// check for an exact match
	matchedFeature, exactMatch := f.features[name]
	if exactMatch && len(containing) < 2 {
		fmt.Fprintf(w, name+"\t"+matchedFeature)
		w.Write([]byte("\n"))
		w.Flush()
		return
	}

	// from https://golang.org/pkg/text/tabwriter/
	if len(containing) < 3 {
		lowDistance = append(lowDistance, containing...)
		containing = []string{} // clear
	}
	if len(containing) > 0 {
		fmt.Fprintf(w, strings.Join(containing, "\n"))
	} else if len(lowDistance) > 0 {
		fmt.Fprintf(w, strings.Join(lowDistance, "\n"))
	} else {
		fmt.Fprintf(w, fmt.Sprintf("failed to find any features for %s", name))
	}
	w.Write([]byte("\n"))
	w.Flush()
}

// SetCmd the feature's seq in the database (or create if it isn't in the feature db)
func (f *FeatureDB) SetCmd(cmd *cobra.Command, args []string) {
	if len(args) < 2 {
		stderr.Fatalf("expecting two args: a features name and sequence. %d passed\n", len(args))
	}

	name := args[0]
	seq := args[1]

	if len(args) > 2 {
		name = strings.Join(args[:len(args)-1], " ")
		seq = args[len(args)-1]
	}

	featureFile, err := os.Open(config.FeatureDB)
	if err != nil {
		stderr.Fatal(err)
	}

	// https://golang.org/pkg/bufio/#example_Scanner_lines
	var output strings.Builder
	updated := false
	scanner := bufio.NewScanner(featureFile)
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

	if err := featureFile.Close(); err != nil {
		stderr.Fatal(err)
	}

	if err := ioutil.WriteFile(config.FeatureDB, []byte(output.String()), 0644); err != nil {
		stderr.Fatal(err)
	}

	if updated {
		fmt.Printf("updated %s in the features database\n", name)
	}

	// update in memory
	f.features[name] = seq
}

// DeleteCmd the feature from the database
func (f *FeatureDB) DeleteCmd(cmd *cobra.Command, args []string) {
	if len(args) < 1 {
		stderr.Fatalf("expecting one arg: a features name. %d passed\n", len(args))
	}

	name := args[0]
	if len(args) > 1 {
		name = strings.Join(args, " ")
	}

	if _, contained := f.features[name]; !contained {
		fmt.Printf("failed to find %s in the features database\n", name)
	}

	featureFile, err := os.Open(config.FeatureDB)
	if err != nil {
		stderr.Fatal(err)
	}

	// https://golang.org/pkg/bufio/#example_Scanner_lines
	var output strings.Builder
	deleted := false
	scanner := bufio.NewScanner(featureFile)
	for scanner.Scan() {
		columns := strings.Split(scanner.Text(), "	")
		if columns[0] != name {
			output.WriteString(scanner.Text())
		} else {
			deleted = true
		}
	}

	if err := featureFile.Close(); err != nil {
		stderr.Fatal(err)
	}

	if err := ioutil.WriteFile(config.FeatureDB, []byte(output.String()), 0644); err != nil {
		stderr.Fatal(err)
	}

	// delete from memory
	delete(f.features, name)

	if deleted {
		fmt.Printf("deleted %s from the features database\n", name)
	} else {
		fmt.Printf("failed to find %s in the features database\n", name)
	}
}

// ld compares two strings and returns the levenshtein distance between them.
// This was copied verbatim from https://github.com/spf13/cobra
func ld(s, t string, ignoreCase bool) int {
	if ignoreCase {
		s = strings.ToUpper(s)
		t = strings.ToUpper(t)
	}
	d := make([][]int, len(s)+1)
	for i := range d {
		d[i] = make([]int, len(t)+1)
	}
	for i := range d {
		d[i][0] = i
	}
	for j := range d[0] {
		d[0][j] = j
	}
	for j := 1; j <= len(t); j++ {
		for i := 1; i <= len(s); i++ {
			if s[i-1] == t[j-1] {
				d[i][j] = d[i-1][j-1]
			} else {
				min := d[i-1][j]
				if d[i][j-1] < min {
					min = d[i][j-1]
				}
				if d[i-1][j-1] < min {
					min = d[i-1][j-1]
				}
				d[i][j] = min + 1
			}
		}
	}
	return d[len(s)][len(t)]
}
