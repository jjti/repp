package defrag

import (
	"reflect"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

func Test_node_distTo(t *testing.T) {
	c := config.New()

	type fields struct {
		id         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		other *node
	}
	tests := []struct {
		name       string
		fields     fields
		args       args
		wantBpDist int
	}{
		{
			"negative distance with overlap",
			fields{
				id:         "1",
				uniqueID:   "1",
				start:      0,
				end:        40,
				assemblies: []assembly{},
			},
			args{
				other: &node{
					id:         "2",
					uniqueID:   "2",
					start:      20,
					end:        60,
					assemblies: []assembly{},
					conf:       c,
				},
			},
			-20,
		},
		{
			"positive distance without overlap",
			fields{
				start: 0,
				end:   40,
			},
			args{
				other: &node{
					start: 60,
					end:   100,
				},
			},
			20,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &node{
				id:         tt.fields.id,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
				conf:       c,
			}
			if gotBpDist := n.distTo(tt.args.other); gotBpDist != tt.wantBpDist {
				t.Errorf("node.distTo() = %v, want %v", gotBpDist, tt.wantBpDist)
			}
		})
	}
}

func Test_node_synthDist(t *testing.T) {
	c := config.New()
	c.Synthesis.MaxLength = 100

	type fields struct {
		id         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		other *node
	}
	tests := []struct {
		name           string
		fields         fields
		args           args
		wantSynthCount int
	}{
		{
			"synth dist is 0 with overlap",
			fields{
				start: 0,
				end:   40,
			},
			args{
				other: &node{
					start: 20,
					end:   60,
					conf:  c,
				},
			},
			0,
		},
		{
			"synth dist is 1 without overlap",
			fields{
				start: 0,
				end:   40,
			},
			args{
				other: &node{
					start: 60,
					end:   80,
					conf:  c,
				},
			},
			1,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &node{
				id:         tt.fields.id,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
				conf:       c,
			}
			if gotSynthCount := n.synthDist(tt.args.other); gotSynthCount != tt.wantSynthCount {
				t.Errorf("node.synthDist() = %v, want %v", gotSynthCount, tt.wantSynthCount)
			}
		})
	}
}

func Test_node_costTo(t *testing.T) {
	c := config.New()
	c.Fragments.MinHomology = 20
	c.PCR.BPCost = 0.03
	c.Synthesis.Cost = map[int]config.SynthCost{
		100000: {
			Fixed:   false,
			Dollars: 0.05,
		},
	}

	type fields struct {
		id         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		other *node
	}
	tests := []struct {
		name     string
		fields   fields
		args     args
		wantCost float64
	}{
		{
			"just cost of PCR of new node if they overlap",
			fields{
				start: 0,
				end:   50,
			},
			args{
				other: &node{
					start: 20,
					end:   100,
					conf:  c,
				},
			},
			1.38,
		},
		{
			"cost of synthesis if they don't overlap",
			fields{
				start: 0,
				end:   50,
			},
			args{
				other: &node{
					start: 80,
					end:   120,
					conf:  c,
				},
			},
			(20.0 + 30.0) * 0.05,
		},
		{
			"cost to self should just be for PCR",
			fields{
				start: n1.start,
				end:   n1.end,
			},
			args{
				other: n1,
			},
			1.38,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &node{
				id:         tt.fields.id,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
				conf:       c,
			}
			if gotCost := n.costTo(tt.args.other); gotCost != tt.wantCost {
				t.Errorf("node.costTo() = %v, want %v", gotCost, tt.wantCost)
			}
		})
	}
}

func Test_node_reach(t *testing.T) {
	c := config.New()

	c.Fragments.MinHomology = 2

	n11 := &node{
		start: 0,
		end:   10,
		conf:  c,
	}
	n12 := &node{
		start: 5,
		end:   15,
		conf:  c,
	}
	n13 := &node{
		start: 6,
		end:   16,
		conf:  c,
	}
	n14 := &node{
		start: 7,
		end:   17,
		conf:  c,
	}
	n15 := &node{
		start: 15,
		end:   20,
		conf:  c,
	}
	n16 := &node{
		start: 16,
		end:   21,
		conf:  c,
	}
	n17 := &node{
		start: 17,
		end:   22,
		conf:  c,
	}

	type fields struct {
		id         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		nodes      []*node
		i          int
		synthCount int
	}
	tests := []struct {
		name          string
		fields        fields
		args          args
		wantReachable []int
	}{
		{
			"gather all reachable nodes",
			fields{
				start: n11.start,
				end:   n11.end,
			},
			args{
				[]*node{n11, n12, n13, n14, n15, n16, n17},
				0,
				2, // limit to 3 "synthable" nodes
			},
			// get all the over-lappable nodes plus two more that
			// can be synthesized to
			[]int{1, 2, 3, 4, 5},
		},
		{
			"return nothing at end",
			fields{
				start: n11.start,
				end:   n11.end,
			},
			args{
				[]*node{n11, n12, n13, n14, n15, n16, n17},
				6,
				2, // limit to 3 "synthable" nodes
			},
			// nothing is reachable
			[]int{},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &node{
				id:         tt.fields.id,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
				conf:       c,
			}
			if gotReachable := n.reach(tt.args.nodes, tt.args.i, tt.args.synthCount); !reflect.DeepEqual(gotReachable, tt.wantReachable) {
				t.Errorf("node.reach() = %v, want %v", gotReachable, tt.wantReachable)
			}
		})
	}
}

func Test_node_synthTo(t *testing.T) {
	c := config.New()
	c.Fragments.MinHomology = 2
	c.Synthesis.MinLength = 4
	c.Synthesis.MaxLength = 100

	type fields struct {
		id         string
		seq        string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		next *node
		seq  string
	}
	tests := []struct {
		name             string
		fields           fields
		args             args
		wantSynthedFrags []Fragment
	}{
		{
			"return nothing when there's enough overlap",
			fields{
				start: 10,
				end:   16,
			},
			args{
				next: &node{
					start: 13,
					end:   20,
					conf:  c,
				},
			},
			nil,
		},
		{
			"return synthetic fragments when there's no overlap",
			fields{
				id:    "first",
				start: 10,
				end:   16,
			},
			args{
				next: &node{
					id:    "second",
					start: 25,
					end:   30,
					conf:  c,
				},
				seq: "TGCTGACTGTGGCGGGTGAGCTTAGGGGGCCTCCGCTCCAGCTCGACACCGGGCAGCTGC",
			},
			[]Fragment{
				Fragment{
					ID:   "first-synthetic-1",
					Seq:  "GGTGAGCTTAGGGGG",
					Type: synthetic,
					Cost: 160,
				},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &node{
				id:         tt.fields.id,
				seq:        tt.fields.seq,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
				conf:       c,
			}
			if gotSynthedFrags := n.synthTo(tt.args.next, tt.args.seq); !reflect.DeepEqual(gotSynthedFrags, tt.wantSynthedFrags) {
				t.Errorf("node.synthTo() = %v, want %v", gotSynthedFrags, tt.wantSynthedFrags)
			}
		})
	}
}

func Test_node_fragment(t *testing.T) {
	c := config.New()

	type fields struct {
		id         string
		seq        string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	tests := []struct {
		name   string
		fields fields
		want   Fragment
	}{
		{
			"convert to fragment from node",
			fields{
				id:  "frag1",
				seq: "atgctgac",
			},
			Fragment{
				ID:    "frag1",
				Seq:   "ATGCTGAC",
				Entry: "frag1",
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &node{
				id:         tt.fields.id,
				seq:        tt.fields.seq,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
				conf:       c,
			}
			if got := n.fragment(); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("node.fragment() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_new(t *testing.T) {
	c := config.New()

	type args struct {
		m    match
		seqL int
	}
	tests := []struct {
		name string
		args args
		want *node
	}{
		{
			"create a node from a match",
			args{
				m: match{
					entry: "testMatch",
					seq:   "atgctagctagtg",
					start: 0,
					end:   12,
				},
				seqL: 50,
			},
			&node{
				id:         "testMatch",
				seq:        "ATGCTAGCTAGTG",
				uniqueID:   "0testMatch",
				start:      0,
				end:        12,
				assemblies: nil,
				conf:       c,
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := newNode(tt.args.m, tt.args.seqL, c); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("new() = %+v, want %+v", got, tt.want)
			}
		})
	}
}

func Test_node_junction(t *testing.T) {
	type fields struct {
		id         string
		uniqueID   string
		seq        string
		start      int
		end        int
		db         string
		assemblies []assembly
		cost       float64
		url        string
		primers    []Primer
		conf       *config.Config
	}
	type args struct {
		other       *node
		minHomology int
		maxHomology int
	}
	tests := []struct {
		name         string
		fields       fields
		args         args
		wantJunction string
	}{
		{
			"find a junction",
			fields{
				seq: "ATGACACGATACGTTATCCACACAGATAGTAGAGATGACACAGATACGAGCGCCTTGAATAACGTACTCATCTCTA",
			},
			args{
				other: &node{
					seq: "GAGCGCCTTGAATAACGTACTCATCTCTATACATTCTCGTGCGCATCACTCTGAATGTACAAGCAACCCAAGAGGGCTGAGCCTGGACTCAGCTGGTTCCTGGGTGAGCTCGAGACTCGGGGTGACAGCTCTTCA",
				},
				minHomology: 5,
				maxHomology: 40,
			},
			"GAGCGCCTTGAATAACGTACTCATCTCTA",
		},
		{
			"fails to find a junction with mismatch",
			fields{
				seq: "ATGACACGATACGTTATCCACACAGATAGTAGAGATGACACAGATACGAGCGCCTTGAATAACGTACTCATCTCTAg", // <- extra g at the end that prevents this from being a junction
			},
			args{
				other: &node{
					seq: "GAGCGCCTTGAATAACGTACTCATCTCTATACATTCTCGTGCGCATCACTCTGAATGTACAAGCAACCCAAGAGGGCTGAGCCTGGACTCAGCTGGTTCCTGGGTGAGCTCGAGACTCGGGGTGACAGCTCTTCA",
				},
				minHomology: 5,
				maxHomology: 40,
			},
			"",
		},
		{
			"finds a small junction",
			fields{
				seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
			},
			args{
				other: &node{
					seq: "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
				},
				minHomology: 5,
				maxHomology: 40,
			},
			"AGCTAGCATCG",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &node{
				id:         tt.fields.id,
				uniqueID:   tt.fields.uniqueID,
				seq:        tt.fields.seq,
				start:      tt.fields.start,
				end:        tt.fields.end,
				db:         tt.fields.db,
				assemblies: tt.fields.assemblies,
				cost:       tt.fields.cost,
				url:        tt.fields.url,
				primers:    tt.fields.primers,
				conf:       tt.fields.conf,
			}
			if gotJunction := n.junction(tt.args.other, tt.args.minHomology, tt.args.maxHomology); gotJunction != tt.wantJunction {
				t.Errorf("node.junction() = %v, want %v", gotJunction, tt.wantJunction)
			}
		})
	}
}
