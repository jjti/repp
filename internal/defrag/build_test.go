package defrag

import (
	"reflect"
	"sort"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

func Test_build(t *testing.T) {
	c := config.New()

	c.Synthesis.MaxLength = 5
	c.Fragments.MinHomology = 1
	c.Fragments.MaxCount = 3

	// ignore cost for now
	c.PCR.BPCost = 0
	c.Synthesis.Cost = map[int]config.SynthCost{
		100000: {
			Fixed:   true,
			Dollars: float32(0),
		},
	}

	// input
	n21 := node{
		id:       "1",
		uniqueID: "1",
		start:    0,
		end:      10,
		conf:     &c,
	}
	n22 := node{
		id:       "2",
		uniqueID: "2",
		start:    5,
		end:      15,
		conf:     &c,
	}
	n23 := node{
		id:       "3",
		uniqueID: "1",
		start:    14,
		end:      30,
		conf:     &c,
	}
	n24 := node{
		id:       "4",
		uniqueID: "2",
		start:    28,
		end:      35,
		conf:     &c,
	}

	// output
	a1 := assembly{
		// wrap w/o synthesis
		nodes:  []node{n21, n22, n23},
		synths: 0,
	}
	a2 := assembly{
		// wrap w/ synthesis
		nodes:  []node{n21, n23},
		synths: 1,
	}
	a3 := assembly{
		// wrap w/o synthsis
		nodes:  []node{n22, n23, n24},
		synths: 0,
	}

	type args struct {
		nodes []node
	}
	tests := []struct {
		name           string
		args           args
		wantAssemblies []assembly
	}{
		{
			"test building of assemblies",
			args{
				nodes: []node{n21, n22, n23, n24},
			},
			[]assembly{a1, a2, a3},
		},
	}

	// turn list list of assemblies into a list of node ids
	getNodeSet := func(as []assembly) (assemblyIDs []string) {
		for _, a := range as {
			thisID := ""
			for _, n := range a.nodes {
				thisID += n.id
			}
			assemblyIDs = append(assemblyIDs, thisID)
		}

		sort.Strings(assemblyIDs)

		return
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assemblies := build(tt.args.nodes, c.Fragments.MaxCount)

			// concatenate node ids together
			actualIds := getNodeSet(assemblies)

			// wanted id sets
			wantedIds := getNodeSet(tt.wantAssemblies)

			if !reflect.DeepEqual(actualIds, wantedIds) {
				t.Errorf("build() = %v, want %v", actualIds, wantedIds)
			}
		})
	}
}
