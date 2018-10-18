package io

import (
	"encoding/json"
	"io/ioutil"
	"log"
	"time"

	"github.com/jjtimmons/decvec/config"
	"github.com/jjtimmons/decvec/internal/dvec"
)

// Solution is a single solution to build up the target vector
type Solution struct {
	// count is the number of fragments in this solution
	Count int `json:"count"`

	// cost estimated from the primer and sequence lengths
	Cost float32 `json:"costDollars"`

	// Fragments used to build this solution
	Fragments []dvec.Fragment `json:"fragments"`
}

// Out is the result output from this assembly
type Out struct {
	// unix
	Time int64 `json:"time"`

	// result
	Result []Solution `json:"result"`
}

// Write a slice of possible assemblies to the fs at the output path
func Write(filename string, assemblies [][]dvec.Fragment) {
	c := config.New()

	results := []Solution{}
	for _, assembly := range assemblies {
		var cost float32
		for _, frag := range assembly {
			if frag.Type == dvec.PCR {
				cost += c.PCR.BPCost * float32(len(frag.Primers[0].Seq))
				cost += c.PCR.BPCost * float32(len(frag.Primers[1].Seq))
			} else if frag.Type == dvec.Synthetic {
				cost += c.Synthesis.BPCost * float32(len(frag.Seq))
			}
		}

		results = append(results, Solution{
			Count:     len(assembly),
			Cost:      cost,
			Fragments: assembly,
		})
	}
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
