package defrag

import (
	"reflect"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

// create t
var (
	c = config.New()

	n1 = &node{
		uniqueID: "1",
		start:    0,
		end:      50,
		conf:     &c,
	}
	n2 = &node{
		uniqueID: "2",
		start:    20,
		end:      80,
		conf:     &c,
	}
	n3 = &node{
		uniqueID: "3",
		start:    60,
		end:      100,
		conf:     &c,
	}
)

func Test_assembly_add(t *testing.T) {
	c.Fragments.MaxCount = 5
	c.Synthesis.MaxLength = 100
	c.Synthesis.Cost = map[int]config.SynthCost{
		100000: {
			Fixed:   true,
			Dollars: 0.0,
		},
	}

	// create the nodes for testing
	type fields struct {
		nodes    []*node
		cost     float64
		synths   int
		maxCount int
	}
	type args struct {
		n *node
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
				nodes:  []*node{n1},
				cost:   0,
				synths: 0,
			},
			args{
				n: n2,
			},
			assembly{
				nodes:  []*node{n1, n2},
				cost:   n1.costTo(n2),
				synths: 0,
			},
			true,
			false,
		},
		{
			"add with synthesis",
			fields{
				nodes:  []*node{n1},
				cost:   10.0,
				synths: 0,
			},
			args{
				n: n3,
			},
			assembly{
				nodes:  []*node{n1, n3},
				cost:   10.0 + n1.costTo(n3),
				synths: 1,
			},
			true,
			false,
		},
		{
			"add with completion",
			fields{
				nodes:  []*node{n1, n2, n3},
				cost:   10.0,
				synths: 0,
			},
			args{
				n: n1,
			},
			assembly{
				nodes:  []*node{n1, n2, n3, n1},
				cost:   10.0,
				synths: 0,
			},
			true,
			true,
		},
		{
			"add with completion requiring synthesis",
			fields{
				nodes:  []*node{n1, n2, n3},
				cost:   16.4,
				synths: 0,
			},
			args{
				// a node that's too far away for straightforward annealing
				n: &node{
					uniqueID: n1.uniqueID,
					start:    n3.start + c.Synthesis.MaxLength,
					end:      n3.end + c.Synthesis.MaxLength,
				},
			},
			assembly{
				nodes: []*node{n1, n2, n3, &node{
					uniqueID: n1.uniqueID,
					start:    n3.start + c.Synthesis.MaxLength,
					end:      n3.end + c.Synthesis.MaxLength,
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
				nodes:  []*node{n1, n2, n3, n2, n3},
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
			gotNewAssembly, gotCreated, gotComplete := a.add(tt.args.n, 5)
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
		nodes    []*node
		cost     float64
		synths   int
		maxCount int
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
				[]*node{n1, n2},
				0.0,
				0,
				5,
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
				[]*node{n1, n2},
				0.0,
				0,
				5,
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
		nodes    []*node
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
				nodes:  []*node{n1, n2},
				synths: 0,
			},
			2,
		},
		{
			"length with synths",
			fields{
				nodes:  []*node{n1, n2},
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
	c := config.New()
	c.Fragments.MinHomology = 5
	c.Fragments.MaxHomology = 50
	c.Synthesis.MaxLength = 1000
	c.PCR.P3MaxPenalty = 40.0

	// All of these ids correspond to entires in the test BLAST db
	f1 := &node{
		id:       "gnl|addgene|113726(circular)",
		uniqueID: "01",
		start:    5,
		end:      119 + 5,
		seq:      "GTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTA",
		conf:     &c,
	}

	f2 := &node{
		id:       "gnl|addgene|85039.1",
		uniqueID: "02",
		start:    102,
		end:      102 + 141,
		seq:      "CTCGTAGGGCTTGCCTTCGCCCTCGGATGTGCACTTGAAGTGGTGGTTGTTCACGGTGCCCTCCATGTACAGCTTCATGTGCATGTTCTCCTTGATCAGCTCGCTCATAGGTCCAGGGTTCTCCTCCACGTCTCCAGCCTG",
		conf:     &c,
	}

	f3 := &node{
		id:       "gnl|addgene|39412.1",
		uniqueID: "03",
		start:    102 + 121,
		end:      102 + 121 + 224,
		seq:      "CACGTCTCCAGCCTGCTTCAGCAGGCTGAAGTTAGTAGCTCCGCTTCCGGATCCCCCGGGGAGCATGTCAAGGTCAAAATCGTCAAGAGCGTCAGCAGGCAGCATATCAAGGTCAAAGTCGTCAAGGGCATCGGCTGGGAgCATGTCTAAgTCAAAATCGTCAAGGGCGTCGGCCGGCCCGCCGCTTTcgcacGCCCTGGCAATCGAGATGCTGGACAGGCATC",
		conf:     &c,
	}

	// starts 150bp past the end of f3, will require synthesis
	f4 := &node{
		id:       "gnl|addgene|107006(circular)",
		uniqueID: "01",
		start:    102 + 121 + 224 + 150,
		end:      102 + 121 + 224 + 150 + 57,
		seq:      "CGTCGGCCGGCCCGCCGCTTTcgcacGCCCTGGCAATCGAGATGCTGGACAGGCATC",
		conf:     &c,
	}

	type fields struct {
		nodes    []*node
		cost     float64
		synths   int
		maxCount int
	}
	type args struct {
		seq string
	}
	tests := []struct {
		name   string
		fields fields
		args   args
	}{
		{
			"fill in an assembly",
			fields{
				nodes: []*node{f1, f2, f3, f4},
			},
			args{
				seq: "GGCCGCAATAAAATATCTTTATTTTCATTACATCTGTGTGTTGGTTTTTTGTGTGAATCGATAGTACTAACATGACCACCTTGATCTTCATGGTCTGGGTGCCCTCGTAGGGCTTGCCTTCGCCCTCGGATGTGCACTTGAAGTGGTGGTTGTTCACGGTGCCCTCCATGTACAGCTTCATGTGCATGTTCTCCTTGATCAGCTCGCTCATAGGTCCAGGGTTCTCCTCCACGTCTCCAGCCTGCTTCAGCAGGCTGAAGTTAGTAGCTCCGCTTCCGGATCCCCCGGGGAGCATGTCAAGGTCAAAATCGTCAAGAGCGTCAGCAGGCAGCATATCAAGGTCAAAGTCGTCAAGGGCATCGGCTGGGAgCATGTCTAAgTCAAAATCGTCAAGGGCGTCGGCCGGCCCGCCGCTTTcgcacGCCCTGGCAATCGAGATGCTGGACAGGCATCATACCCACTTCTGCCCCCTGGAAGGCGAGTCATGGCAAGACTTTCTGCGGAACAACGCCAAGTCATTCCGCTGTGCTCTCCTCTCACATCGCGACGGGGCTAAAGTGCATCTCGGCACCCGCCCAACAGAGAAACAGTACGAAACCCTGGAAAATCAGCTCGCGTTCCTGTGTCAGCAAGGCTTCTCCCTGGAGAACGCACTGTACGCTCTGTCCGCCGTGGGCCACTTTACACTGGGCTGCGTATTGGAGGATCAGGAGCATCAAGTAGCAAAAGAGGAAAGAGAGACACCTACCACCGATTCTATGCCTGACTGTGGCGGGTGAGCTTAGGGGGCCTCCGCTCCAGCTCGACACCGGGCAGCTGCTGAAGATCGCGAAGAGAGGGGGAGTAACAGCGGTAGAGGCAGTGCACGCCTGGCGCAATGCGCTCACCGGGGCCCCCTTGAACCTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAGCAATGGGGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTTTAACAAAATATTAACGCTTACAATTTAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAA",
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := &assembly{
				nodes:  tt.fields.nodes,
				cost:   tt.fields.cost,
				synths: tt.fields.synths,
			}
			frags, err := a.fill(tt.args.seq, &c)

			if err != nil {
				t.Error(err)
			}

			pcrCount := 0
			synthCount := 0

			for _, f := range frags {
				if f.Type == pcr {
					pcrCount++
				} else if f.Type == synthetic {
					synthCount++
				}
			}

			if pcrCount != 3 {
				t.Errorf("expected 3 pcr fragments, got %d", pcrCount)
			}

			if synthCount != 1 {
				t.Errorf("expected 1 synthetic fragment, got %d", synthCount)
			}
		})
	}
}

func Test_assembly_duplicates(t *testing.T) {
	type fields struct {
		nodes  []*node
		cost   float64
		synths int
	}
	type args struct {
		nodes       []*node
		minHomology int
		maxHomology int
	}
	tests := []struct {
		name   string
		fields fields
		args   args
		want   bool
	}{
		{
			"no false positive",
			fields{},
			args{
				nodes: []*node{
					&node{
						seq: "ATACCTACTATGGATGACGTAGCAAC",
					},
					&node{
						seq: "AGCAACTCGTTGATATCCACGTA",
					},
					&node{
						seq: "CCACGTAGGTGCATGATGAGATGA",
					},
					&node{
						seq: "TGAGATGATCTACTGTATACCTACT",
					},
				},
				minHomology: 5,
				maxHomology: 10,
			},
			false,
		},
		{
			"assembly with a self-annealing node",
			fields{},
			args{
				nodes: []*node{
					&node{
						seq: "CAGATGACGATGGCAACTGAGATGAGACCAGATGACGATG", // <- node (if much larger) has the chance to circularize
					},
					&node{
						seq: "CAGATGACGATGTCGTTGATATACCTACTGGAGAGCACAG",
					},
					&node{
						seq: "TGGAGAGCACAGATGGATGACGTAATGATGATGACCGCAAC",
					},
					&node{
						seq: "ACCGCAACTCGTTGATATACCTACTCAGATGACGAT",
					},
				},
				minHomology: 5,
				maxHomology: 20,
			},
			true,
		},
		{
			"assembly with a duplicate junction",
			fields{},
			args{
				nodes: []*node{
					&node{
						seq: "ATGATGCCACGTGCAACTGAGATGAGACCAGATGACGATG", // <- same junction
					},
					&node{
						seq: "CAGATGACGATGTCGTTGATATACCTACTGGAGAGCACAG",
					},
					&node{
						seq: "TGGAGAGCACAGATGGATGACGTAATGACAGATGACGATG", // <- same junction
					},
					&node{
						seq: "CAGATGACGATGACCGCAACTCGTTGATGATGCCAC",
					},
				},
				minHomology: 5,
				maxHomology: 20,
			},
			true,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := &assembly{
				nodes:  tt.fields.nodes,
				cost:   tt.fields.cost,
				synths: tt.fields.synths,
			}
			if got := a.duplicates(tt.args.nodes, tt.args.minHomology, tt.args.maxHomology); got != tt.want {
				t.Errorf("assembly.duplicates() = %v, want %v", got, tt.want)
			}
		})
	}
}
