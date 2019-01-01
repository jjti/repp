package defrag

import (
	"path"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
)

func Test_vector(t *testing.T) {
	in, _ := filepath.Abs(path.Join("..", "..", "test", "109049.addgene.fa"))
	out, _ := filepath.Abs(path.Join("..", "..", "bin", "109049.addgene.json"))
	dbs, _ := filepath.Abs(path.Join("..", "..", "test", "db", "db"))

	// https://stackoverflow.com/a/50880663
	cmd := &cobra.Command{}
	cmd.Flags().String("in", in, "")
	cmd.Flags().String("out", out, "")
	cmd.Flags().String("dbs", dbs, "")
	cmd.Flags().Bool("addgene", true, "")

	type args struct {
		cmd  *cobra.Command
		args []string
	}
	tests := []struct {
		name string
		args args
	}{
		{
			"end to end design of a test vector",
			args{
				cmd: cmd,
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			Vector(tt.args.cmd, tt.args.args)
		})
	}
}

// if an input fragment being built is exactly the same as one in a DB, it should be used
// as is and without PCR or any preparation
func Test_vector_single_vector(t *testing.T) {
	c := config.New()
	in := path.Join("..", "..", "test", "109049.addgene.fa")
	out := path.Join("..", "..", "bin", "109049.output.json")

	assemblies := vector(&flags{
		in:  in,
		out: out,
	}, c) // use addgene database

	if len(assemblies) != 1 {
		t.Fatal("failed to return the pareto optimal solution: 109049 alone")
	}

	if len(assemblies[0]) != 1 || !strings.Contains(assemblies[0][0].ID, "109049") {
		t.Fatal("failed to use 109049 to build the vector")
	}

	if assemblies[0][0].Type != circular {
		t.Fatalf("failed to recognize 109049 as a Type.Vector, was %d", assemblies[0][0].Type)
	}
}

func Test_countMaps(t *testing.T) {
	a1 := assembly{
		nodes: []*Frag{
			n1, n1,
		},
		cost: 11.0,
	}
	a2 := assembly{
		nodes: []*Frag{
			n1, n2, n1,
		},
		cost: 12.5,
	}
	a3 := assembly{
		nodes: []*Frag{
			n2, n3, n2,
		},
		cost: 12.0,
	}
	a4 := assembly{
		nodes: []*Frag{
			n1, n2, n3, n1,
		},
		cost: 10.0,
	}
	a5 := assembly{
		nodes: []*Frag{
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
