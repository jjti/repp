package assemble

import (
	"reflect"
	"sort"
	"testing"
)

func Test_build(t *testing.T) {
	conf.Synthesis.MaxLength = 5
	conf.Fragments.MinHomology = 1
	conf.Fragments.MaxCount = 3

	// ignore cost for now
	conf.PCR.BPCost = 0
	conf.Synthesis.BPCost = 0

	// input
	n21 := node{
		id:       "1",
		uniqueID: "1",
		start:    0,
		end:      10,
	}
	n22 := node{
		id:       "2",
		uniqueID: "2",
		start:    5,
		end:      15,
	}
	n23 := node{
		id:       "3",
		uniqueID: "1",
		start:    14,
		end:      30,
	}
	n24 := node{
		id:       "4",
		uniqueID: "2",
		start:    28,
		end:      35,
	}

	// output
	a1 := assembly{
		// wrap w/o synthesis
		nodes:  []node{n21, n22},
		synths: 0,
	}
	a2 := assembly{
		// wrap w/ synthesis
		nodes:  []node{n21},
		synths: 1,
	}
	a3 := assembly{
		// wrap w/o synthsis
		nodes:  []node{n22, n23},
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
			assemblies := build(tt.args.nodes)

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
