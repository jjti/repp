package defrag

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
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

// NewFeatureDB returns a new copy of the features db
func NewFeatureDB() *FeatureDB {
	features := make(map[string]string)

	featureFile, err := os.Open(config.FeatureDB)
	if err != nil {
		log.Fatal(err)
	}

	// https://golang.org/pkg/bufio/#example_Scanner_lines
	scanner := bufio.NewScanner(featureFile)
	for scanner.Scan() {
		columns := strings.Split(scanner.Text(), "	")
		features[columns[0]] = columns[1] // feature name = feature seq
	}

	if err := featureFile.Close(); err != nil {
		log.Fatal(err)
	}

	return &FeatureDB{features: features}
}

// Create adds an additional feature to the db (if it's not in it already)
func (f *FeatureDB) Create(cmd *cobra.Command, args []string) {
	if len(args) < 2 {
		log.Fatalf("expecting two args: a features name and sequence. %d passed\n", len(args))
	}

	name := args[0]
	seq := args[1]

	argNames := []string{"create", "add", "read", "find", "update", "change", "delete", "remove"}
	for _, argName := range argNames {
		if name == argName {
			log.Fatalf("cannot create a feature named %s. invalid feature names: %s\n", name, strings.Join(argNames, ", "))
		}
	}

	if len(args) > 2 {
		name = strings.Join(args[:len(args)-1], " ")
		seq = args[len(args)-1]
	}

	if _, contained := f.features[name]; contained {
		log.Fatalf("cannot create feature %s, already in %s. try update\n", name, config.FeatureDB)
	}

	// https://golang.org/pkg/os/#example_OpenFile_append
	featureFile, err := os.OpenFile(config.FeatureDB, os.O_APPEND|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatal(err)
	}

	if _, err := featureFile.Write([]byte(fmt.Sprintf("%s	%s\n", name, seq))); err != nil {
		log.Fatal(err)
	}

	if err := featureFile.Close(); err != nil {
		log.Fatal(err)
	}

	fmt.Printf("created %s in the features database\n", name)

	// add to features in the map
	f.features[name] = seq
}

// Read returns features that are similar in name to the feature name requested.
// if multiple feature names include the feature name, they are all returned.
// otherwise a list of feature names are returned (those beneath a levenshtein distance cutoff)
func (f *FeatureDB) Read(cmd *cobra.Command, args []string) {
	if len(args) < 1 {
		log.Fatalf("expecting one arg: a features name. %d passed\n", len(args))
	}

	name := args[0]
	if len(args) > 1 {
		name = strings.Join(args, " ")
	}

	ldCutoff := 5
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
	w := tabwriter.NewWriter(os.Stdout, 0, 0, 1, ' ', tabwriter.TabIndent)
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

// Update the feature's seq in the database (or create if it isn't in the feature db)
func (f *FeatureDB) Update(cmd *cobra.Command, args []string) {
	if len(args) < 2 {
		log.Fatalf("expecting two args: a features name and sequence. %d passed\n", len(args))
	}

	name := args[0]
	seq := args[1]

	if len(args) > 2 {
		name = strings.Join(args[:len(args)-1], " ")
		seq = args[len(args)-1]
	}

	if _, contained := f.features[name]; !contained {
		f.Create(nil, []string{name, seq}) // it doesn't exist, create
	}

	featureFile, err := os.Open(config.FeatureDB)
	if err != nil {
		log.Fatal(err)
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

	if err := featureFile.Close(); err != nil {
		log.Fatal(err)
	}

	if err := ioutil.WriteFile(config.FeatureDB, []byte(output.String()), 0644); err != nil {
		log.Fatal(err)
	}

	if updated {
		fmt.Printf("updated %s in the features database\n", name)
	}

	// update in memory
	f.features[name] = seq
}

// Delete the feature from the database
func (f *FeatureDB) Delete(cmd *cobra.Command, args []string) {
	if len(args) < 1 {
		log.Fatalf("expecting one arg: a features name. %d passed\n", len(args))
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
		log.Fatal(err)
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
		log.Fatal(err)
	}

	if err := ioutil.WriteFile(config.FeatureDB, []byte(output.String()), 0644); err != nil {
		log.Fatal(err)
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
