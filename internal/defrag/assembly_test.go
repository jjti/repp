package defrag

import (
	"reflect"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

func Test_assembly_add(t *testing.T) {
	c := config.New()

	c.FragmentsMaxCount = 5
	c.PCRMaxEmbedLength = 0
	c.SynthesisMaxLength = 100
	c.SynthesisFragmentCost = map[int]config.SynthCost{
		100000: {
			Fixed: true,
			Cost:  0.0,
		},
	}

	sl := 100

	n1 := &Frag{
		uniqueID: "1",
		start:    0,
		end:      50,
		conf:     c,
	}
	n2 := &Frag{
		uniqueID: "2",
		start:    20,
		end:      80,
		conf:     c,
	}
	n3 := &Frag{
		uniqueID: "3",
		start:    60,
		end:      100,
		conf:     c,
	}
	n4 := &Frag{
		uniqueID: "1",
		start:    100,
		end:      150,
		conf:     c,
	}

	// create the nodes for testing
	type fields struct {
		nodes    []*Frag
		cost     float64
		synths   int
		maxCount int
	}
	type args struct {
		n *Frag
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
				nodes:  []*Frag{n1},
				cost:   0,
				synths: 0,
			},
			args{
				n: n2,
			},
			assembly{
				frags:  []*Frag{n1, n2},
				cost:   n1.costTo(n2),
				synths: 0,
			},
			true,
			false,
		},
		{
			"add with synthesis",
			fields{
				nodes:  []*Frag{n1},
				cost:   10.0,
				synths: 0,
			},
			args{
				n: n3,
			},
			assembly{
				frags:  []*Frag{n1, n3},
				cost:   10.0 + n1.costTo(n3),
				synths: 1,
			},
			true,
			false,
		},
		{
			"add with completion/circularization",
			fields{
				nodes:  []*Frag{n1, n2, n3},
				cost:   10.0,
				synths: 0,
			},
			args{
				n: n4,
			},
			assembly{
				frags:  []*Frag{n1, n2, n3},
				cost:   10.0,
				synths: 0,
			},
			true,
			true,
		},
		{
			"add with completion requiring synthesis",
			fields{
				nodes:  []*Frag{n1, n2},
				cost:   16.4,
				synths: 0,
			},
			args{
				// a Frag that's too far away for straightforward annealing
				n: &Frag{
					uniqueID: n1.uniqueID,
					start:    n3.start + c.SynthesisMaxLength,
					end:      n3.end + c.SynthesisMaxLength,
				},
			},
			assembly{
				frags:  []*Frag{n1, n2},
				cost:   16.4,
				synths: 1,
			},
			true,
			true,
		},
		{
			"don't exceed fragment limit",
			fields{
				nodes:  []*Frag{n1, n2, n3, n2, n3},
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
				frags:  tt.fields.nodes,
				cost:   tt.fields.cost,
				synths: tt.fields.synths,
			}
			gotNewAssembly, gotCreated, gotComplete := a.add(tt.args.n, 5, sl)
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
func Test_assembly_len(t *testing.T) {
	c := config.New()

	n1 := &Frag{
		uniqueID: "1",
		start:    0,
		end:      50,
		conf:     c,
	}
	n2 := &Frag{
		uniqueID: "2",
		start:    20,
		end:      80,
		conf:     c,
	}

	type fields struct {
		nodes    []*Frag
		cost     float64
		synths   int
		maxCount int
	}
	tests := []struct {
		name   string
		fields fields
		want   int
	}{
		{
			"length without synths",
			fields{
				nodes:  []*Frag{n1, n2},
				synths: 0,
			},
			2,
		},
		{
			"length with synths",
			fields{
				nodes:  []*Frag{n1, n2},
				synths: 2,
			},
			4,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := &assembly{
				frags:  tt.fields.nodes,
				cost:   tt.fields.cost,
				synths: tt.fields.synths,
			}
			if got := a.len(); got != tt.want {
				t.Errorf("assembly.len() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_countMaps(t *testing.T) {
	c := config.New()
	n1 := &Frag{
		uniqueID: "1",
		start:    0,
		end:      50,
		conf:     c,
	}
	n2 := &Frag{
		uniqueID: "2",
		start:    20,
		end:      80,
		conf:     c,
	}
	n3 := &Frag{
		uniqueID: "3",
		start:    60,
		end:      100,
		conf:     c,
	}

	a1 := assembly{
		frags: []*Frag{
			n1, n1,
		},
		cost: 11.0,
	}
	a2 := assembly{
		frags: []*Frag{
			n1, n2, n1,
		},
		cost: 12.5,
	}
	a3 := assembly{
		frags: []*Frag{
			n2, n3, n2,
		},
		cost: 12.0,
	}
	a4 := assembly{
		frags: []*Frag{
			n1, n2, n3, n1,
		},
		cost: 10.0,
	}
	a5 := assembly{
		frags: []*Frag{
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
			if _, gotParetoSet := groupAssembliesByCount(tt.args.assemblies); !reflect.DeepEqual(gotParetoSet, tt.wantParetoSet) {
				t.Errorf("pareto() = %v, want %v", gotParetoSet, tt.wantParetoSet)
			}
		})
	}
}

func Test_assembly_duplicates(t *testing.T) {
	type fields struct {
		nodes  []*Frag
		cost   float64
		synths int
	}
	type args struct {
		nodes       []*Frag
		minHomology int
		maxHomology int
	}
	tests := []struct {
		name         string
		fields       fields
		args         args
		want         bool
		wantJunction string
	}{
		{
			"no false positive",
			fields{},
			args{
				nodes: []*Frag{
					&Frag{
						Seq: "ATACCTACTATGGATGACGTAGCAAC",
					},
					&Frag{
						Seq: "AGCAACTCGTTGATATCCACGTA",
					},
					&Frag{
						Seq: "CCACGTAGGTGCATGATGAGATGA",
					},
					&Frag{
						Seq: "TGAGATGATCTACTGTATACCTACT",
					},
				},
				minHomology: 5,
				maxHomology: 10,
			},
			false,
			"",
		},
		{
			"assembly with a self-annealing Frag",
			fields{},
			args{
				nodes: []*Frag{
					&Frag{
						Seq: "CAGATGACGATGGCAACTGAGATGAGACCAGATGACGATG", // <- Frag (if much larger) has the chance to circularize
					},
					&Frag{
						Seq: "CAGATGACGATGTCGTTGATATACCTACTGGAGAGCACAG",
					},
					&Frag{
						Seq: "TGGAGAGCACAGATGGATGACGTAATGATGATGACCGCAAC",
					},
					&Frag{
						Seq: "ACCGCAACTCGTTGATATACCTACTCAGATGACGAT",
					},
				},
				minHomology: 5,
				maxHomology: 20,
			},
			true,
			"CAGATGACGATG",
		},
		{
			"assembly with a duplicate junction",
			fields{},
			args{
				nodes: []*Frag{
					&Frag{
						Seq: "ATGATGCCACGTGCAACTGAGATGAGACCAGATGACGATG", // <- same junction
					},
					&Frag{
						Seq: "CAGATGACGATGTCGTTGATATACCTACTGGAGAGCACAG",
					},
					&Frag{
						Seq: "TGGAGAGCACAGATGGATGACGTAATGACAGATGACGATG", // <- same junction
					},
					&Frag{
						Seq: "CAGATGACGATGACCGCAACTCGTTGATGATGCCAC",
					},
				},
				minHomology: 5,
				maxHomology: 20,
			},
			true,
			"CAGATGACGATG",
		},
		{
			"another false positive to avoid",
			fields{},
			args{
				nodes: []*Frag{
					&Frag{
						Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG", // this shouldn't be flagged as anything
					},
					&Frag{
						Seq: "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
					},
					&Frag{
						Seq: "TCGACTAGCTAGAACTGATGCTAGACGTGCTAGCTACA",
					},
				},
				minHomology: 8,
				maxHomology: 20,
			},
			false,
			"",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := &assembly{
				frags:  tt.fields.nodes,
				cost:   tt.fields.cost,
				synths: tt.fields.synths,
			}
			isDuplicate, _, _, duplicateSeq := a.duplicates(tt.args.nodes, tt.args.minHomology, tt.args.maxHomology)

			if isDuplicate != tt.want {
				t.Errorf("assembly.duplicates() = %v, want %v", isDuplicate, tt.want)
			}

			if duplicateSeq != tt.wantJunction {
				t.Errorf("assembly.duplicates() = %v, want %v", duplicateSeq, tt.wantJunction)
			}
		})
	}
}
