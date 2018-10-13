// Package assemble turns blast Matches into building Fragments
package assemble

import (
	"reflect"
	"testing"
)

func Test_filterNode(t *testing.T) {
	// blacklisted node
	bNode := node{
		uniqueID: "b1",
	}

	// other node 1
	oNode1 := node{
		uniqueID: "o1",
	}

	// other node 2
	oNode2 := node{
		uniqueID: "o2",
	}

	// other node 3
	oNode3 := node{
		uniqueID: "o3",
	}

	// assembly w/ the black listed node
	a1 := assembly{
		nodes: []node{bNode, oNode1, oNode2},
	}

	// assembly without the blast listed node
	a2 := assembly{
		nodes: []node{oNode1, oNode2, oNode3},
	}

	type args struct {
		black node
		old   []assembly
	}
	tests := []struct {
		name    string
		args    args
		wantNew []assembly
	}{
		{
			"remove assemblies with blacklisted node",
			args{
				black: bNode,
				old:   []assembly{a1, a2},
			},
			[]assembly{a2},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotNew := filterNode(tt.args.black, tt.args.old); !reflect.DeepEqual(gotNew, tt.wantNew) {
				t.Errorf("filterNode() = %v, want %v", gotNew, tt.wantNew)
			}
		})
	}
}
