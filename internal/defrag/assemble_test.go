// package defrag turns blast Matches into building Fragments
package defrag

import (
	"reflect"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

func Test_countMaps(t *testing.T) {
	a1 := assembly{
		nodes: []*node{
			n1, n1,
		},
		cost: 11.0,
	}
	a2 := assembly{
		nodes: []*node{
			n1, n2, n1,
		},
		cost: 12.5,
	}
	a3 := assembly{
		nodes: []*node{
			n2, n3, n2,
		},
		cost: 12.0,
	}
	a4 := assembly{
		nodes: []*node{
			n1, n2, n3, n1,
		},
		cost: 10.0,
	}
	a5 := assembly{
		nodes: []*node{
			n2, n3, n1, n2,
		},
		cost: 10.5,
	}

	type args struct {
		assemblies []assembly
	}
	tests := []struct {
		name          string
		args          args
		wantParetoSet map[int][]assembly
	}{
		{
			"gen pSet up to 3",
			args{
				assemblies: []assembly{a1, a2, a3, a4, a5},
			},
			map[int][]assembly{
				1: []assembly{a1},
				2: []assembly{a3, a2},
				3: []assembly{a4, a5},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotParetoSet := countMap(tt.args.assemblies); !reflect.DeepEqual(gotParetoSet, tt.wantParetoSet) {
				t.Errorf("pareto() = %v, want %v", gotParetoSet, tt.wantParetoSet)
			}
		})
	}
}

func Test_assembleFWD(t *testing.T) {
	c := config.New()
	c.Fragments.MinHomology = 8
	c.Fragments.MaxHomology = 20

	type args struct {
		inputFragments []Fragment
		conf           *config.Config
	}
	tests := []struct {
		name             string
		args             args
		wantTargetVector Fragment
		wantFragments    []Fragment
	}{
		{
			"fragments with enough overlap",
			args{
				[]Fragment{
					Fragment{
						Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
					},
					Fragment{
						Seq: "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
					},
					Fragment{
						Seq: "TCGACTAGCTAGAACTGATCTAGACGTGCTAGCTACA",
					},
				},
				&c,
			},
			Fragment{
				Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCGACTGATCACTAGCATCGACTAGCTAGAACTGATCTAG",
			},
			[]Fragment{
				Fragment{
					Seq:  "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
					Type: existing,
				},
				Fragment{
					Seq:  "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
					Type: existing,
				},
				Fragment{
					Seq:  "TCGACTAGCTAGAACTGATCTAGACGTGCTAGCTACA",
					Type: existing,
				},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotTargetVector, gotFragments := assembleFWD(tt.args.inputFragments, tt.args.conf)
			if !reflect.DeepEqual(gotTargetVector, tt.wantTargetVector) {
				t.Errorf("assembleFWD() gotTargetVector = %v, want %v", gotTargetVector, tt.wantTargetVector)
			}
			if !reflect.DeepEqual(gotFragments, tt.wantFragments) {
				t.Errorf("assembleFWD() gotFragments = %v, want %v", gotFragments, tt.wantFragments)
			}
		})
	}
}
