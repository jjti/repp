package io

import (
	"encoding/json"
	"io/ioutil"
	"log"
	"sort"
	"time"

	"github.com/jjtimmons/defrag/config"
	"github.com/jjtimmons/defrag/internal/defrag"
)

// Solution is a single solution to build up the target vector
type Solution struct {
	// count is the number of fragments in this solution
	Count int `json:"count"`

	// cost estimated from the primer and sequence lengths
	Cost int `json:"costDollars"`

	// Fragments used to build this solution
	Fragments []defrag.Fragment `json:"fragments"`
}

// Out is the result output from this assembly
type Out struct {
	// unix
	Time int64 `json:"time"`

	// result
	Result []Solution `json:"result"`
}

// Write a slice of possible assemblies to the fs at the output path
func Write(filename string, assemblies [][]defrag.Fragment) {
	c := config.New()

	// calculate final cost of the assembly and fragment count
	results := []Solution{}
	for _, assembly := range assemblies {
		var cost int
		for _, frag := range assembly {
			if frag.Type == defrag.PCR {
				cost += int(c.PCR.BPCost * float32(len(frag.Primers[0].Seq)))
				cost += int(c.PCR.BPCost * float32(len(frag.Primers[1].Seq)))
			} else if frag.Type == defrag.Synthetic {
				cost += int(c.Synthesis.BPCost * float32(len(frag.Seq)))
			}
		}

		results = append(results, Solution{
			Count:     len(assembly),
			Cost:      cost,
			Fragments: assembly,
		})
	}
	// sort results in increasing fragment count order
	sort.Slice(results, func(i, j int) bool {
		return results[i].Count < results[j].Count
	})
	out := Out{
		Time:   time.Now().Unix(),
		Result: results,
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
