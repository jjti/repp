package repp

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/jjtimmons/repp/config"
)

// Solution is a single solution to build up the target plasmid.
type Solution struct {
	// Count is the number of fragments in this solution
	Count int `json:"count"`

	// Cost estimated from the primer and sequence lengths
	Cost float64 `json:"cost"`

	// Fragments used to build this solution
	Fragments []*Frag `json:"fragments"`
}

// Output is a struct containing design results for the assembly.
type Output struct {
	// Target's name. In >example_CDS FASTA its "example_CDS"
	Target string `json:"target"`

	// Target's sequence
	TargetSeq string `json:"seq"`

	// Time, ex: "2018-01-01 20:41:00"
	Time string `json:"time"`

	// Execution is the number of seconds it took to execute the command
	Execution float64 `json:"execution"`

	// PlasmidSynthesisCost is the cost of a full gene synthesis within a plasmid
	// PlasmidSynthesisCost float64 `json:"plasmidSynthesisCost"`

	// InsertSynthesisCost is the cost of just synthesizing the insert with homology for a linearized backbone
	// InsertSynthesisCost float64 `json:"insertSynthesisCost"`

	// Solutions builds
	Solutions []Solution `json:"solutions"`

	// Backbone is the user linearized a backbone fragment
	Backbone *Backbone `json:"backbone,omitempty"`
}

// writeJSON turns a list of solutions into a Solution object and writes to the filename requested.
func writeJSON(
	filename,
	targetName,
	targetSeq string,
	assemblies [][]*Frag,
	insertSeqLength int,
	seconds float64,
	backbone *Backbone,
	conf *config.Config,
) (output []byte, err error) {
	// store save time, using same format as log.Println https://golang.org/pkg/log/#Println
	t := time.Now() // https://gobyexample.com/time-formatting-parsing
	time := fmt.Sprintf(
		"%d/%02d/%02d %02d:%02d:%02d",
		t.Year(), t.Month(), t.Day(), t.Hour(), t.Minute(), t.Second(),
	)
	roundCost := func(cost float64) (float64, error) {
		return strconv.ParseFloat(fmt.Sprintf("%.2f", cost), 64)
	}

	// calculate final cost of the assembly and fragment count
	solutions := []Solution{}
	for _, assembly := range assemblies {
		assemblyCost := 0.0
		assemblyFragmentIDs := make(map[string]bool)
		gibson := false // whether it will be assembled via Gibson assembly
		hasPCR := false // whether there will be a batch PCR

		for _, f := range assembly {
			if f.fragType != linear && f.fragType != circular {
				gibson = true
			}

			if f.fragType == pcr {
				hasPCR = true
			}

			f.Type = f.fragType.String() // freeze fragment type

			if f.URL == "" && f.fragType != synthetic {
				f.URL = parseURL(f.ID, f.db)
			}

			if f.URL != "" {
				f.ID = "" // just log one or the other
			}

			// round to two decimal places
			if f.Cost, err = roundCost(f.cost(true)); err != nil {
				return nil, err
			}

			// if it's already in the assembly, don't count cost twice
			if _, contained := assemblyFragmentIDs[f.ID]; f.ID != "" && contained {
				if f.Cost, err = roundCost(f.cost(false)); err != nil {
					return nil, err // ignore repo procurement costs
				}
			} else {
				assemblyFragmentIDs[f.ID] = true
			}

			// accumulate assembly cost
			assemblyCost += f.Cost
		}

		if gibson {
			assemblyCost += conf.CostGibson + conf.CostTimeGibson
		}

		if hasPCR {
			assemblyCost += conf.CostTimePCR
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

	// get the cost of full synthesis (comes in plasmid)
	// fullSynthCost, err := roundCost(conf.SynthPlasmidCost(insertSeqLength))
	// if err != nil {
	// 	return nil, err
	// }

	// estimate the cost of making the insert, with overhang for the backbone,
	// from the synthesis provider and then the Gibson assembly cost
	insertSynthCost, err := roundCost(conf.SynthFragmentCost(insertSeqLength + conf.FragmentsMinHomology*2))
	if err != nil {
		return nil, err
	}
	insertSynthCost += conf.CostGibson

	if backbone.Seq == "" {
		backbone = nil
	}

	out := Output{
		Time:      time,
		Target:    targetName,
		TargetSeq: strings.ToUpper(targetSeq),
		Execution: seconds,
		Solutions: solutions,
		Backbone:  backbone,
		// PlasmidSynthesisCost: fullSynthCost,
		// InsertSynthesisCost: insertSynthCost,
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

// writeGenbank writes a slice of fragments/features to a genbank output file.
func writeGenbank(filename, name, seq string, frags []*Frag, feats []match) {
	// header row
	d := time.Now().Local()
	h1 := fmt.Sprintf("LOCUS       %s", name)
	h2 := fmt.Sprintf("%d bp DNA      circular      %s\n", len(seq), strings.ToUpper(d.Format("02-Jan-2006")))
	space := strings.Repeat(" ", 81-len(h1+h2))
	header := h1 + space + h2

	// feature rows
	var fsb strings.Builder
	fsb.WriteString("DEFINITION  .\nACCESSION   .\nFEATURES             Location/Qualifiers\n")
	for _, m := range feats {
		cS := ""
		cE := ""
		if !m.forward {
			cS = "complement("
			cE = ")"
		}

		s := (m.queryStart + 1) % len(seq)
		e := (m.queryEnd + 1) % len(seq)

		if s == 0 {
			s = len(seq)
		}
		if e == 0 {
			e = len(seq)
		}

		fsb.WriteString(
			fmt.Sprintf("     misc_feature    %s%d..%d%s\n", cS, s, e, cE) +
				fmt.Sprintf("                     /label=\"%s\"\n", m.entry),
		)
	}

	// origin row
	var ori strings.Builder
	ori.WriteString("ORIGIN\n")
	for i := 0; i < len(seq); i += 60 {
		n := strconv.Itoa(i + 1)
		ori.WriteString(strings.Repeat(" ", 9-len(n)) + n)
		for s := i; s < i+60 && s < len(seq); s += 10 {
			e := s + 10
			if e > len(seq) {
				e = len(seq)
			}
			ori.WriteString(fmt.Sprintf(" %s", seq[s:e]))
		}
		ori.WriteString("\n")
	}
	ori.WriteString("//\n")

	gb := strings.Join([]string{header, fsb.String(), ori.String()}, "")
	err := ioutil.WriteFile(filename, []byte(gb), 0644)
	if err != nil {
		stderr.Fatalln(err)
	}
}
