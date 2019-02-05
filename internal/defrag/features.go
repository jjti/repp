package defrag

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"os"
	"strings"
	"text/tabwriter"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
)

// FeatureDB is a struct for accessing defrags features db
type FeatureDB struct {
	// features is a map between a features name and its sequence
	features map[string]string
}

// FeaturesCmd accepts a cobra commands and assembles a vector containing all the features
func FeaturesCmd(cmd *cobra.Command, args []string) {
	features(parseCmdFlags(cmd, args))
}

// features assembles a vector with all the features requested with the 'defrag features [feature ...]' command
func features(flags *Flags, conf *config.Config) {
	targetFeatures := [][]string{} // slice of tuples [feature name, feature sequence]

	if parsedFeatures, err := read(flags.in, true); err == nil {
		seenFeatures := make(map[string]string) // map feature name to sequence
		// see if the features are in a file (multi-FASTA or features in a Genbank)
		for _, f := range parsedFeatures {
			if seq := seenFeatures[f.ID]; seq != f.Seq {
				stderr.Fatalf("failed to parse features, %s has two different sequences:\n\t%s\n\t%s\n", f.ID, f.Seq, seq)
			}

			targetFeatures = append(targetFeatures, []string{f.ID, f.Seq})
		}
	} else {
		// if the features weren't in a file, try and find each in the features database
		// or one of the databases passed as a source of building fragments
		featureNames := strings.Fields(flags.in)
		if len(featureNames) < 1 {
			stderr.Fatal("no features chosen. see 'defrag features --help'")
		}

		featureDB := NewFeatureDB()
		for _, f := range featureNames {
			if seq, contained := featureDB.features[f]; contained {
				targetFeatures = append(targetFeatures, []string{f, seq})
			} else if dbFrag, err := queryDatabases(f, flags.dbs); err == nil {
				targetFeatures = append(targetFeatures, []string{f, dbFrag.Seq})
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
	if flags.backbone.ID != "" {
		targetFeatures = append(targetFeatures, []string{flags.backbone.ID, flags.backbone.Seq})
	}

	fragsToMatches := make(map[string][]match) // map from building fragments to feature matches
	fragsToFeats := make(map[string][]int)     // list of features (by index) that the frag/vector with unique id contains
	featsToFrags := make(map[int][]string)     // map from feature index to building fragments
	for i, target := range targetFeatures {
		// TODO: don't re-BLAST here. if we already know a features' matches, just return those
		matches, err := blast(target[0], target[1], false, flags.dbs, flags.filters, flags.identity)
		if err != nil {
			stderr.Fatalln(err)
		}

		if _, exists := fragsToFeats[target[0]]; !exists {
			fragsToFeats[target[0]] = []int{i}
		} else {
			fragsToFeats[target[0]] = append(fragsToFeats[target[0]], i)
		}

		featsToFrags[i] = []string{}
		for _, m := range matches {
			if _, exists := fragsToMatches[m.entry]; !exists {
				fragsToMatches[m.entry] = []match{m}
			} else {
				fragsToMatches[m.entry] = append(fragsToMatches[m.entry], m)
			}
			featsToFrags[i] = append(featsToFrags[i], m.entry)
		}
	}

	// step through each feature index and create an assembly

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
	if len(args) < 1 {
		stderr.Fatalf("expecting one arg: a features name. %d passed\n", len(args))
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

	// from https://golang.org/pkg/text/tabwriter/
	w := tabwriter.NewWriter(os.Stdout, 0, 0, 3, ' ', tabwriter.TabIndent)
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

// containNext returns fragments that have the next feature after the last feature
// contained in the one being tested. Returns two slices. One with frags that
// already overlap the frag enough for annealing without preparation, and another
// of fragments that contain the next feature but will need to be annealed via PCR (more expensive)
func containNext(last *Frag, lastF int, featsToFrags map[int][]string, fragsToFeats map[string][]int) (overlap []*Frag, haveNext []*Frag) {
	return nil, nil
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
