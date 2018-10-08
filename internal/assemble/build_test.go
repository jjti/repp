package assemble

import (
	"reflect"
	"testing"
)

func Test_build(t *testing.T) {
	// too small for a node to synth all the way back to itself
	conf.Synthesis.MaxLength = 5
	conf.Fragments.MinHomology = 1

	// input
	n21 := node{
		uniqueID: "1",
		start:    0,
		end:      10,
	}
	n22 := node{
		uniqueID: "2",
		start:    5,
		end:      15,
	}
	n23 := node{
		uniqueID: "3",
		start:    17,
		end:      20,
	}
	n24 := node{
		uniqueID: "1",
		start:    19,
		end:      30,
	}
	n25 := node{
		uniqueID: "2",
		start:    28,
		end:      35,
	}

	// output
	a1 := assembly{
		nodes:  []node{n21, n22, n24},
		synths: 0,
	}
	a2 := assembly{
		nodes: []node{n22, n24, n25},
	}
	a3 := assembly{
		nodes: []node{n21, n22, n23},
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
				nodes: []node{n21, n22, n23},
			},
			[]assembly{a1, a2, a3},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotAssemblies := build(tt.args.nodes); !reflect.DeepEqual(gotAssemblies, tt.wantAssemblies) {
				t.Errorf("build() = %v, want %v", gotAssemblies, tt.wantAssemblies)
			}
		})
	}
}
