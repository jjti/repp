package assemble

import (
	"reflect"
	"testing"

	"github.com/jjtimmons/decvec/internal/dvec"
)

// create t
var (
	n1 = node{
		uniqueID: "1",
		start:    0,
		end:      50,
	}
	n2 = node{
		uniqueID: "2",
		start:    20,
		end:      80,
	}
	n3 = node{
		uniqueID: "3",
		start:    60,
		end:      100,
	}
)

func Test_assembly_add(t *testing.T) {
	conf.Fragments.MaxCount = 5
	conf.Synthesis.MaxLength = 100
	conf.Synthesis.BPCost = 0.08

	// create the nodes for testing
	type fields struct {
		nodes  []node
		cost   float32
		synths int
	}
	type args struct {
		n node
	}
	tests := []struct {
		name            string
		fields          fields
		args            args
		wantNewAssembly assembly
		wantCreated     bool
		wantComplete    bool
	}{
		{
			"add with overlap",
			fields{
				nodes:  []node{n1},
				cost:   0,
				synths: 0,
			},
			args{
				n: n2,
			},
			assembly{
				nodes:  []node{n1, n2},
				cost:   n1.costTo(n2),
				synths: 0,
			},
			true,
			false,
		},
		{
			"add with synthesis",
			fields{
				nodes:  []node{n1},
				cost:   10.0,
				synths: 0,
			},
			args{
				n: n3,
			},
			assembly{
				nodes:  []node{n1, n3},
				cost:   10.0 + n1.costTo(n3),
				synths: 1,
			},
			true,
			false,
		},
		{
			"add with completion",
			fields{
				nodes:  []node{n1, n2, n3},
				cost:   10.0,
				synths: 0,
			},
			args{
				n: n1,
			},
			assembly{
				nodes:  []node{n1, n2, n3, n1},
				cost:   10.0,
				synths: 0,
			},
			true,
			true,
		},
		{
			"add with completion requiring synthesis",
			fields{
				nodes:  []node{n1, n2, n3},
				cost:   10.0,
				synths: 0,
			},
			args{
				// a node that's too far away for straightforward annealing
				n: node{
					uniqueID: n1.uniqueID,
					start:    n3.start + conf.Synthesis.MaxLength,
					end:      n3.end + conf.Synthesis.MaxLength,
				},
			},
			assembly{
				nodes: []node{n1, n2, n3, node{
					uniqueID: n1.uniqueID,
					start:    n3.start + conf.Synthesis.MaxLength,
					end:      n3.end + conf.Synthesis.MaxLength,
				}},
				cost:   16.4,
				synths: 1,
			},
			true,
			true,
		},
		{
			"don't exceed fragment limit",
			fields{
				nodes:  []node{n1, n2, n3, n2, n3},
				cost:   10.0,
				synths: 0,
			},
			args{
				n: n2,
			},
			assembly{},
			false,
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := &assembly{
				nodes:  tt.fields.nodes,
				cost:   tt.fields.cost,
				synths: tt.fields.synths,
			}
			gotNewAssembly, gotCreated, gotComplete := a.add(tt.args.n)
			if !reflect.DeepEqual(gotNewAssembly, tt.wantNewAssembly) {
				t.Errorf("assembly.add() gotNewAssembly = %v, want %v", gotNewAssembly, tt.wantNewAssembly)
			}
			if gotCreated != tt.wantCreated {
				t.Errorf("assembly.add() gotCreated = %v, want %v", gotCreated, tt.wantCreated)
			}
			if gotComplete != tt.wantComplete {
				t.Errorf("assembly.add() gotComplete = %v, want %v", gotComplete, tt.wantComplete)
			}
		})
	}
}

func Test_assembly_contains(t *testing.T) {
	type fields struct {
		nodes  []node
		cost   float32
		synths int
	}
	type args struct {
		n node
	}
	tests := []struct {
		name            string
		fields          fields
		args            args
		wantIsContained bool
	}{
		{
			"contains node",
			fields{
				[]node{n1, n2},
				0.0,
				0,
			},
			args{
				node{
					uniqueID: "1",
				},
			},
			true,
		},
		{
			"doesn't contain node",
			fields{
				[]node{n1, n2},
				0.0,
				0,
			},
			args{
				node{
					uniqueID: "3",
				},
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := &assembly{
				nodes:  tt.fields.nodes,
				cost:   tt.fields.cost,
				synths: tt.fields.synths,
			}
			if gotIsContained := a.contains(tt.args.n); gotIsContained != tt.wantIsContained {
				t.Errorf("assembly.contains() = %v, want %v", gotIsContained, tt.wantIsContained)
			}
		})
	}
}

func Test_assembly_len(t *testing.T) {
	type fields struct {
		nodes  []node
		cost   float32
		synths int
	}
	tests := []struct {
		name   string
		fields fields
		want   int
	}{
		{
			"length without synths",
			fields{
				nodes:  []node{n1, n2},
				synths: 0,
			},
			2,
		},
		{
			"length with synths",
			fields{
				nodes:  []node{n1, n2},
				synths: 2,
			},
			4,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := &assembly{
				nodes:  tt.fields.nodes,
				cost:   tt.fields.cost,
				synths: tt.fields.synths,
			}
			if got := a.len(); got != tt.want {
				t.Errorf("assembly.len() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_assembly_fill(t *testing.T) {
	type fields struct {
		nodes  []node
		cost   float32
		synths int
	}
	type args struct {
		seq string
	}
	tests := []struct {
		name          string
		fields        fields
		args          args
		wantBlacklist node
		wantFrags     []dvec.Fragment
	}{
		{
			"fills in a basic assembly",
			fields{
				nodes: []node{},
			},
			args{},
			node{},
			nil,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := &assembly{
				nodes:  tt.fields.nodes,
				cost:   tt.fields.cost,
				synths: tt.fields.synths,
			}
			gotBlacklist, gotFrags := a.fill(tt.args.seq)
			if !reflect.DeepEqual(gotBlacklist, tt.wantBlacklist) {
				t.Errorf("assembly.fill() gotBlacklist = %v, want %v", gotBlacklist, tt.wantBlacklist)
			}
			if !reflect.DeepEqual(gotFrags, tt.wantFrags) {
				t.Errorf("assembly.fill() gotFrags = %v, want %v", gotFrags, tt.wantFrags)
			}
		})
	}
}
