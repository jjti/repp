package defrag

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"path/filepath"
	"regexp"
	"sort"
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
	Cost float64 `json:"dollars"`

	// Fragments used to build this solution
	Fragments []*Frag `json:"fragments"`
}

// Output is a struct containing design results for the assembly
type Output struct {
	// Time, ex:
	// "2006-01-02 15:04:05.999999999 -0700 MST"
	// https://golang.org/pkg/time/#Time.String
	Time string `json:"time"`

	// Target sequence
	Target string `json:"target"`

	// Solutions builds
	Solutions []Solution `json:"solutions"`

	// SynthInsertCost cost of a full synthesis
	SynthInsertCost float64 `json:"synthInsertCost"`
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
	for i, ID := range ids {
		fragments = append(fragments, Frag{
			ID:  ID,
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
	// calculate final cost of the assembly and fragment count
	solutions := []Solution{}
	for _, assembly := range assemblies {
		assemblyCost := 0.0
		for _, f := range assembly {
			f.Dollars = f.cost()
			assemblyCost += f.Dollars
		}

		solutions = append(solutions, Solution{
			Count:     len(assembly),
			Cost:      roundCost(assemblyCost),
			Fragments: assembly,
		})
	}
	// sort solutions in increasing fragment count order
	sort.Slice(solutions, func(i, j int) bool {
		return solutions[i].Count < solutions[j].Count
	})
	out := Output{
		Time:            time.Now().String(),
		Target:          strings.ToUpper(target.Seq),
		Solutions:       solutions,
		SynthInsertCost: conf.SynthCost(insertSeqLength),
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
func roundCost(cost float64) float64 {
	// roundedString := fmt.Sprintf("%.2f", cost)
	// result, _ := strconv.ParseFloat(roundedString, 64)
	return cost
}
