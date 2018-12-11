package defrag

import (
	"encoding/json"
	"io/ioutil"
	"log"
	"sort"
	"time"

	"github.com/jjtimmons/defrag/config"
)

var (
	conf = config.New()
)

// Solution is a single solution to build up the target vector
type Solution struct {
	// count is the number of fragments in this solution
	Count int `json:"count"`

	// cost estimated from the primer and sequence lengths
	Cost float32 `json:"costDollars"`

	// Fragments used to build this solution
	Fragments []Fragment `json:"fragments"`
}

// Out is the result output from this assembly
type Out struct {
	// local time, ex:
	// "2006-01-02 15:04:05.999999999 -0700 MST"
	// https://golang.org/pkg/time/#Time.String
	Time string `json:"time"`

	// target sequence
	Target string `json:"target"`

	// solution builds
	Solutions []Solution `json:"solutions"`
}

// Write a slice of pareto optimal assemblies to the fs at the output path
func Write(filename string, target Fragment, assemblies [][]Fragment) {
	// calculate final cost of the assembly and fragment count
	solutions := []Solution{}
	for _, assembly := range assemblies {
		solutions = append(solutions, Solution{
			Count:     len(assembly),
			Cost:      solutionCost(assembly, conf.PCR.BPCost, conf.SynthCost),
			Fragments: assembly,
		})
	}
	// sort solutions in increasing fragment count order
	sort.Slice(solutions, func(i, j int) bool {
		return solutions[i].Count < solutions[j].Count
	})
	out := Out{
		Time:      time.Now().String(),
		Target:    target.Seq,
		Solutions: solutions,
	}

	output, err := json.MarshalIndent(out, "", "  ")
	if err != nil {
		log.Fatalf("Failed to serialize the output data: %v", err)
	}

	err = ioutil.WriteFile(filename, output, 0666)
	if err != nil {
		log.Fatalf("Failed to write the results to the file system: %v", err)
	}
}

// solutionCost returns an estimated cost of an assembly (a list of fragments)
// be estimating the total cost of the primers and any synthesized basepairs
func solutionCost(frags []Fragment, primerBP float32, synthCost func(int) float32) (cost float32) {
	for _, frag := range frags {
		if frag.Type == PCR {
			cost += primerBP * float32(len(frag.Primers[0].Seq))
			cost += primerBP * float32(len(frag.Primers[1].Seq))
		} else if frag.Type == Synthetic {
			cost += synthCost(len(frag.Seq))
		}
	}

	return
}