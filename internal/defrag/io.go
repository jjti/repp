package defrag

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/jjtimmons/defrag/config"
)

// Payload is the JSON shape of an input payload from the API endpoint
type Payload struct {
	// Target sequence that we want to build
	Target string `json:"target"`

	// Fragments that should be assembled
	Fragments []Frag `json:"fragments"`

	// Addgene is whether to use the Addgene repository
	Addgene bool `json:"addgene"`

	// IGEM is whether to use the IGEM repository
	IGEM bool `json:"igem"`

	// Backbone to insert the inserts into
	Backbone string `json:"backbone"`

	// Enzyme to cut the backbone with (name)
	Enzyme string `json:"enzyme"`
}

// Solution is a single solution to build up the target vector
type Solution struct {
	// Count is the number of fragments in this solution
	Count int `json:"count"`

	// Cost estimated from the primer and sequence lengths
	Cost float64 `json:"cost"`

	// Fragments used to build this solution
	Fragments []*Frag `json:"fragments"`
}

// Output is a struct containing design results for the assembly
type Output struct {
	// Target's name. In >example_CDS FASTA its "example_CDS"
	Target string `json:"target"`

	// Target's sequence
	TargetSeq string `json:"seq"`

	// Time, ex:
	// "2018-01-01 20:41:00"
	Time string `json:"time"`

	// FullSynthesisCost cost of a full synthesis
	FullSynthesisCost float64 `json:"fullSynthesisCost"`

	// Solutions builds
	Solutions []Solution `json:"solutions"`
}

// read a FASTA file (by its path on local FS) to a slice of Fragments
func read(path string) (fragments []Frag, err error) {
	if !filepath.IsAbs(path) {
		path, err = filepath.Abs(path)
		if err != nil {
			return nil, fmt.Errorf("failed to create path to FASTA file: %s", err)
		}
	}

	dat, err := ioutil.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("failed to read input FASTA path: %s", err)
	}
	file := string(dat)

	// split by newlines
	lines := strings.Split(file, "\n")

	// read in the fragments
	var headerIndices []int
	var ids []string
	for i, line := range lines {
		if strings.HasPrefix(line, ">") {
			headerIndices = append(headerIndices, i)
			ids = append(ids, line[1:])
		}
	}

	// create a regex for cleaning the sequence
	var unwantedChars = regexp.MustCompile(`(?im)[^atgc]|\W`)

	// accumulate the sequences from between the headers
	var seqs []string
	for i, headerIndex := range headerIndices {
		nextLine := len(lines)
		if i < len(headerIndices)-1 {
			nextLine = headerIndices[i+1]
		}
		seqLines := lines[headerIndex+1 : nextLine]
		seqJoined := strings.Join(seqLines, "")
		seq := unwantedChars.ReplaceAllString(seqJoined, "")
		seqs = append(seqs, seq)
	}

	// build and return the new fragments
	for i, id := range ids {
		fragments = append(fragments, Frag{
			ID:  id,
			Seq: seqs[i],
		})
	}

	// opened and parsed file but found nothing
	if len(fragments) < 1 {
		return fragments, fmt.Errorf("failed to parse fragment(s) from %s", file)
	}

	return
}

// write a slice of assemblies to the fs at the output path
//
// filename is the output file to write to
// target is the vector we tried to assemble
// assemblies are the solutions that can build up the target vector
func write(filename string, target Frag, assemblies [][]*Frag, insertSeqLength int, conf *config.Config) (output []byte, err error) {
	// store save time, using same format as log.Println https://golang.org/pkg/log/#Println
	t := time.Now() // https://gobyexample.com/time-formatting-parsing
	time := fmt.Sprintf("%d/%02d/%02d %02d:%02d:%02d", t.Year(), t.Month(), t.Day(), t.Hour(), t.Minute(), t.Second())

	// calculate final cost of the assembly and fragment count
	solutions := []Solution{}
	for _, assembly := range assemblies {
		assemblyCost := 0.0
		for _, f := range assembly {
			// freeze fragment type
			f.Type = f.fragType.String()

			// round to two decimal places
			f.Cost, err = roundCost(f.cost())
			if err != nil {
				return nil, err
			}

			// add to actual assembly cost
			assemblyCost += f.Cost
		}

		solutionCost, err := roundCost(assemblyCost)
		if err != nil {
			return nil, err
		}

		solutions = append(solutions, Solution{
			Count:     len(assembly),
			Cost:      solutionCost,
			Fragments: assembly,
		})
	}

	// sort solutions in increasing fragment count order
	sort.Slice(solutions, func(i, j int) bool {
		return solutions[i].Count < solutions[j].Count
	})

	// get the cost of full synthesis
	fullSynthCost, err := roundCost(conf.SynthGeneCost(insertSeqLength))
	if err != nil {
		return nil, err
	}

	out := Output{
		Time:              time,
		Target:            target.ID,
		TargetSeq:         strings.ToUpper(target.Seq),
		Solutions:         solutions,
		FullSynthesisCost: fullSynthCost,
	}

	output, err = json.MarshalIndent(out, "", "  ")
	if err != nil {
		return output, fmt.Errorf("failed to serialize output: %v", err)
	}

	if err = ioutil.WriteFile(filename, output, 0666); err != nil {
		return output, fmt.Errorf("failed to write the output: %v", err)
	}
	return output, nil
}

// roundCost returns a float for cost to 2 decimal places
func roundCost(cost float64) (float64, error) {
	roundedString := fmt.Sprintf("%.2f", cost)
	return strconv.ParseFloat(roundedString, 64)
}
