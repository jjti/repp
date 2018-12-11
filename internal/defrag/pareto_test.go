package defrag

import (
	"reflect"
	"testing"
)

func Test_pareto(t *testing.T) {
	a1 := assembly{
		nodes: []node{
			n1,
		},
		cost: 11.0,
	}
	a2 := assembly{
		nodes: []node{
			n1, n2,
		},
		cost: 12.5,
	}
	a3 := assembly{
		nodes: []node{
			n2, n3,
		},
		cost: 12.0,
	}
	a4 := assembly{
		nodes: []node{
			n1, n2, n3,
		},
		cost: 10.0,
	}
	a5 := assembly{
		nodes: []node{
			n2, n3, n1,
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
			if gotParetoSet := pareto(tt.args.assemblies); !reflect.DeepEqual(gotParetoSet, tt.wantParetoSet) {
				t.Errorf("pareto() = %v, want %v", gotParetoSet, tt.wantParetoSet)
			}
		})
	}
}
