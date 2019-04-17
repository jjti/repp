package rvec

import (
	"fmt"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
	"text/tabwriter"

	"github.com/spf13/cobra"
)

// Annotate is for annotating a vector sequence given the features in the feature database.
// If an output path is provided, the annotated vector is writen to that file. Otherwise,
// the feature matches are written to stdout.
func Annotate(cmd *cobra.Command, args []string) {
	output, _ := cmd.Flags().GetString("out")

	identity, err := cmd.Flags().GetInt("identity")
	if err != nil {
<<<<<<< HEAD
		identity = 100 // might be something other than `rvec vector`
=======
		identity = 100 // might be something other than `defrag vector`
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
	}

	p := inputParser{}
	filters, err := cmd.Flags().GetString("exclude")
	excludeFilters := []string{}
	if err == nil {
		excludeFilters = p.getFilters(filters)
	}

	name := ""
	query := "" // the vector sequence that we're querying
	if len(args) > 0 {
		query = args[0]
	} else {
		in, err := cmd.Flags().GetString("in")
		if in == "" || err != nil {
			cmd.Help()
			stderr.Fatalln("must pass a file with a vector sequence or the vector sequence as an argument.")
		}

		frags, err := read(in, false)
		if err != nil {
			stderr.Fatalln(err)
		}
		name = frags[0].ID
		query = frags[0].Seq
	}

	enclosed, _ := cmd.Flags().GetBool("enclosed")

	annotate(name, query, output, identity, excludeFilters, enclosed)
}

// annotate is for executing blast against the query sequence.
func annotate(name, seq, output string, identity int, filters []string, enclosed bool) {
	handleErr := func(err error) {
		if err != nil {
			stderr.Fatalln(err)
		}
	}

	in, err := ioutil.TempFile(blastnDir, name+"in-*")
	handleErr(err)
	defer os.Remove(in.Name())

	out, err := ioutil.TempFile(blastnDir, name+"out-*")
	handleErr(err)
	defer os.Remove(out.Name())

	// create a subject file with all the blast features
	fDB := NewFeatureDB()
	featIndex := 0
	var featureSubjects strings.Builder
	indexToFeature := make(map[int]string)
	for feat, seq := range fDB.features {
		indexToFeature[featIndex] = feat
		featureSubjects.WriteString(fmt.Sprintf(">%d\n%s\n", featIndex, seq))
		featIndex++
	}
	subjectFile, err := ioutil.TempFile(blastnDir, name+"features-*")
	handleErr(err)
	defer os.Remove(subjectFile.Name())

	_, err = subjectFile.WriteString(featureSubjects.String())
	handleErr(err)

	b := &blastExec{
		in:       in,
		out:      out,
		name:     name,
		subject:  subjectFile.Name(),
		seq:      seq,
		identity: identity,
		circular: true,
	}

	handleErr(b.input())
	handleErr(b.runAgainst())
	features, err := b.parse(filters)
	handleErr(err)

	// get rid of features that start past the zero index, wrap that those that go around it
	// get rid of features matches that aren't 100% of the feature in the feature database
	var cleanedFeatures []match
	for _, f := range features {
		if f.queryStart >= len(seq) {
			continue
		}

		featureIndex, _ := strconv.Atoi(f.entry)
		f.entry = indexToFeature[featureIndex]
		if len(f.seq) < len(fDB.features[f.entry]) {
			continue
		}

		f.queryEnd %= len(seq)
		if f.queryEnd == 0 {
			f.queryEnd = len(seq)
		}
		cleanedFeatures = append(cleanedFeatures, f)
	}
	features = cleanedFeatures

	sortMatches(features)
	if !enclosed {
		features = properize(features)
	}

	if output != "" {
		writeGenbank(output, name, seq, []*Frag{}, features)
	} else {
		tw := tabwriter.NewWriter(os.Stdout, 0, 4, 3, ' ', 0)
		fmt.Fprintf(tw, "\nfeatures (%d)\tstart\tend\t\n", len(features))
		for _, feat := range features {
			fmt.Fprintf(tw, "%s\t%d\t%d\t\n", feat.entry, feat.queryStart+1, feat.queryEnd+1)
		}
		tw.Flush()
	}
}
