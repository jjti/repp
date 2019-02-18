package defrag

import (
	"path"
	"reflect"
	"strings"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

func Test_Frag_distTo(t *testing.T) {
	c := config.New()

	type fields struct {
		ID         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		other *Frag
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
				ID:         "1",
				uniqueID:   "1",
				start:      0,
				end:        40,
				assemblies: []assembly{},
			},
			args{
				other: &Frag{
					ID:         "2",
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
				other: &Frag{
					start: 60,
					end:   100,
				},
			},
			20,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &Frag{
				ID:         tt.fields.ID,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
				conf:       c,
			}
			if gotBpDist := n.distTo(tt.args.other); gotBpDist != tt.wantBpDist {
				t.Errorf("Frag.distTo() = %v, want %v", gotBpDist, tt.wantBpDist)
			}
		})
	}
}

func Test_Frag_synthDist(t *testing.T) {
	c := config.New()

	c.PCRMaxEmbedLength = 0
	c.SynthesisMaxLength = 100

	type fields struct {
		ID         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		other *Frag
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
				other: &Frag{
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
				other: &Frag{
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
			n := &Frag{
				ID:         tt.fields.ID,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
				conf:       c,
			}
			if gotSynthCount := n.synthDist(tt.args.other); gotSynthCount != tt.wantSynthCount {
				t.Errorf("Frag.synthDist() = %v, want %v", gotSynthCount, tt.wantSynthCount)
			}
		})
	}
}

func Test_Frag_costTo(t *testing.T) {
	c := config.New()
	c.FragmentsMinHomology = 20
	c.PCRBPCost = 0.03
	c.SynthesisFragmentCost = map[int]config.SynthCost{
		100000: {
			Fixed: false,
			Cost:  0.05,
		},
	}

	n1 := &Frag{
		uniqueID: "1",
		start:    0,
		end:      50,
		conf:     c,
	}

	type fields struct {
		ID         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		other *Frag
	}
	tests := []struct {
		name     string
		fields   fields
		args     args
		wantCost float64
	}{
		{
			"just cost of PCR of new Frag if they overlap",
			fields{
				start: 0,
				end:   50,
			},
			args{
				other: &Frag{
					start: 20,
					end:   100,
					conf:  c,
				},
			},
			1.5,
		},
		{
			"cost of synthesis if they don't overlap",
			fields{
				start: 0,
				end:   50,
			},
			args{
				other: &Frag{
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
			1.5,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &Frag{
				ID:         tt.fields.ID,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
				conf:       c,
			}
			if gotCost := n.costTo(tt.args.other); gotCost != tt.wantCost {
				t.Errorf("Frag.costTo() = %v, want %v", gotCost, tt.wantCost)
			}
		})
	}
}

func Test_Frag_reach(t *testing.T) {
	c := config.New()

	c.FragmentsMinHomology = 2
	c.PCRMaxEmbedLength = 0

	n11 := &Frag{
		start: 0,
		end:   10,
		conf:  c,
	}
	n12 := &Frag{
		start: 5,
		end:   15,
		conf:  c,
	}
	n13 := &Frag{
		start: 6,
		end:   16,
		conf:  c,
	}
	n14 := &Frag{
		start: 7,
		end:   17,
		conf:  c,
	}
	n15 := &Frag{
		start: 15,
		end:   20,
		conf:  c,
	}
	n16 := &Frag{
		start: 16,
		end:   21,
		conf:  c,
	}
	n17 := &Frag{
		start: 17,
		end:   22,
		conf:  c,
	}

	type fields struct {
		ID         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		nodes      []*Frag
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
				[]*Frag{n11, n12, n13, n14, n15, n16, n17},
				0,
				2, // limit to 2 "synthable" nodes
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
				[]*Frag{n11, n12, n13, n14, n15, n16, n17},
				6,
				2, // limit to 3 "synthable" nodes
			},
			// nothing is reachable
			[]int{},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &Frag{
				ID:         tt.fields.ID,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
				conf:       c,
			}
			if gotReachable := n.reach(tt.args.nodes, tt.args.i, tt.args.synthCount); !reflect.DeepEqual(gotReachable, tt.wantReachable) {
				t.Errorf("Frag.reach() = %v, want %v", gotReachable, tt.wantReachable)
			}
		})
	}
}

func Test_new(t *testing.T) {
	c := config.New()

	type args struct {
		m match
	}
	tests := []struct {
		name string
		args args
		want *Frag
	}{
		{
			"create a Frag from a match",
			args{
				m: match{
					entry:      "testMatch",
					uniqueID:   "0testMatch",
					seq:        "atgctagctagtg",
					queryStart: 0,
					queryEnd:   12,
				},
			},
			&Frag{
				ID:         "testMatch",
				fragType:   existing,
				Seq:        "ATGCTAGCTAGTG",
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
			if got := newFrag(tt.args.m, c); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("new() = %+v, want %+v", got, tt.want)
			}
		})
	}
}

func Test_Frag_junction(t *testing.T) {
	type fields struct {
		ID         string
		uniqueID   string
		Seq        string
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
		other       *Frag
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
				Seq: "ATGACACGATACGTTATCCACACAGATAGTAGAGATGACACAGATACGAGCGCCTTGAATAACGTACTCATCTCTA",
			},
			args{
				other: &Frag{
					Seq: "GAGCGCCTTGAATAACGTACTCATCTCTATACATTCTCGTGCGCATCACTCTGAATGTACAAGCAACCCAAGAGGGCTGAGCCTGGACTCAGCTGGTTCCTGGGTGAGCTCGAGACTCGGGGTGACAGCTCTTCA",
				},
				minHomology: 5,
				maxHomology: 40,
			},
			"GAGCGCCTTGAATAACGTACTCATCTCTA",
		},
		{
			"fails to find a junction with mismatch",
			fields{
				Seq: "ATGACACGATACGTTATCCACACAGATAGTAGAGATGACACAGATACGAGCGCCTTGAATAACGTACTCATCTCTAg", // <- extra g at the end that prevents this from being a junction
			},
			args{
				other: &Frag{
					Seq: "GAGCGCCTTGAATAACGTACTCATCTCTATACATTCTCGTGCGCATCACTCTGAATGTACAAGCAACCCAAGAGGGCTGAGCCTGGACTCAGCTGGTTCCTGGGTGAGCTCGAGACTCGGGGTGACAGCTCTTCA",
				},
				minHomology: 5,
				maxHomology: 40,
			},
			"",
		},
		{
			"finds a small junction",
			fields{
				Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
			},
			args{
				other: &Frag{
					Seq: "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
				},
				minHomology: 5,
				maxHomology: 40,
			},
			"AGCTAGCATCG",
		},
		{
			"find an off-target junction",
			fields{
				Seq: "ATGATGCCACGTGCAACTGAGATGAGACCAGATGACGATG",
			},
			args{
				other: &Frag{
					Seq: "CAGATGACGATGACCGCAACTCGTTGATGATGCCAC",
				},
				minHomology: 5,
				maxHomology: 20,
			},
			"CAGATGACGATG",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &Frag{
				ID:         tt.fields.ID,
				uniqueID:   tt.fields.uniqueID,
				Seq:        tt.fields.Seq,
				start:      tt.fields.start,
				end:        tt.fields.end,
				db:         tt.fields.db,
				assemblies: tt.fields.assemblies,
				URL:        tt.fields.url,
				Primers:    tt.fields.primers,
				conf:       tt.fields.conf,
			}
			if gotJunction := n.junction(tt.args.other, tt.args.minHomology, tt.args.maxHomology); gotJunction != tt.wantJunction {
				t.Errorf("Frag.junction() = %v, want %v", gotJunction, tt.wantJunction)
			}
		})
	}
}

// this is little more than a deprecation test right now
func Test_setPrimers(t *testing.T) {
	c := config.New()

	c.FragmentsMinHomology = 20
	c.FragmentsMaxHomology = 80
	c.PCRP3MaxPenalty = 50.0
	c.PCRMaxEmbedLength = 10
	db := path.Join("..", "..", "test", "db", "db")

	type args struct {
		last *Frag
		next *Frag
		Seq  string
	}
	tests := []struct {
		name        string
		n           Frag
		args        args
		wantPrimers []string // two primers seqs, first is fwd, second is rev
		wantErr     bool
		wantStart   int
		wantEnd     int
	}{
		{
			"create without added homology",
			Frag{
				ID:    "gnl|addgene|85039.2",
				start: 0,
				end:   1050,
				conf:  c,
				db:    db,
			},
			args{
				last: &Frag{ // close enough that no homology should be added
					ID:    "last|fragment",
					start: -50,
					end:   20,
					conf:  c,
				},
				next: &Frag{
					ID:    "next|fragment",
					start: 1030,
					end:   1090,
					conf:  c,
				},
				Seq: "CAGTCAATCTTTCACAAATTTTGTaATCCAGAGGTTGATTATCGATAAGCTTGATATCGAATTCATAACttCgTATAGCATACATTATACGAAGTTATTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTAAGGCAAATAACTTCGTATAGGATACTTTATACGAAGTTATGCTAGCGCCACCATGTCTAGACTGGACAAGAGCAAAGTCATAAACTCTGCTCTGGAATTACTCAATGAAGTCGGTATCGAAGGCCTGACGACAAGGAAACTCGCTCAAAAGCTGGGAGTTGAGCAGCCTACCCTGTACTGGCACGTGAAGAACAAGCGGGCCCTGCTCGATGCCCTGGCAATCGAGATGCTGGACAGGCATCATACCCACTTCTGCCCCCTGGAAGGCGAGTCATGGCAAGACTTTCTGCGGAACAACGCCAAGTCATTCCGCTGTGCTCTCCTCTCACATCGCGACGGGGCTAAAGTGCATCTCGGCACCCGCCCAACAGAGAAACAGTACGAAACCCTGGAAAATCAGCTCGCGTTCCTGTGTCAGCAAGGCTTCTCCCTGGAGAACGCACTGTACGCTCTGTCCGCCGTGGGCCACTTTACACTGGGCTGCGTATTGGAGGATCAGGAGCATCAAGTAGCAAAAGAGGAAAGAGAGACACCTACCACCGATTCTATGCCCCCACTTCTGAGACAAGCAATTGAGCTGTTCGACCATCAGGGAGCCGAACCTGCCTTCCTTTTCGGCCTGGAACTAATCATATGTGGCCTGGAGAAACAGCTAAAGTGCGAAAGCGGCGGGCCGGCCGACGCCCTTGACGATTTTGACTTAGACATGCTCCCAGCCGATGCCCTTGACGACTTTGACCTTGATATGCTGCCTGCTGACGCTCTTGACGATTTtGACCTTGACATGCTCCCCGGGGGATCCGGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGAgACGTGGAGGAGAACCCTGGACCTATGAGCGAGCTGATCAAGGAGAACATGCACATGAAGCTGTACATGGAgggc",
			},
			[]string{
				"CAGTCAATCTTTCACAAATTTTGT",
				"ACAGCTTCATGTGCATGTTCTC", // rev-comp: ACATGAAGCTGTACATGGAGGG
			},
			false,
			0,
			1049,
		},
		{
			"create with added homology",
			Frag{
				ID:    "add_homology",
				start: 500,
				end:   800,
				conf:  c,
				db:    db,
			},
			args{
				last: &Frag{
					ID:    "last",
					start: 200,
					end:   500,
					conf:  c,
				},
				next: &Frag{
					ID:    "next",
					start: 795,
					end:   900,
					conf:  c,
				},
				Seq: "TGCTGACTGTGGCGGGTGAGCTTAGGGGGCCTCCGCTCCAGCTCGACACCGGGCAGCTGCTGAAGATCGCGAAGAGAGGGGGAGTAACAGCGGTAGAGGCAGTGCACGCCTGGCGCAATGCGCTCACCGGGGCCCCCTTGAACCTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAGCAATGGGGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCCACGGGCTGACTCCCGATCAAGTTGTAGCGATTGCGTCGCATGACGGAGGGAAACAAGCATTGGAGACTGTCCAACGGCTCCTTCCCGTGTTGTGTCAAGCCCACGGTTTGACGCCTGCACAAGTGGTCGCCATCGCCAGCCATGATGGCGGTAAGCAGGCGCTGGAAACAGTACAGCGCCTGCTGCCTGTACTGTGCCAGGATCATGGACTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAATAATAACGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCCACGGGCTGACTCCCGATCAAGTTGTAGCGATTGCGTCGCATGACGGAGGGAAACAAGCATTGGAGACTGTCCAACGGCTCCTTCCCGTGTTGTGTCAAGCCCACGGTTTGACGCCTGCACAAGTGGTCGCCATCGCCAACAACAACGGCGGTAAGCAGGCGCTGGAAACAGTACAGCGCCTGCTGCCTGTACTGTGCCAGGATCATGGACTGACCCCAGACCAGGTAGTCGCAATCGCGTCGAACATTGGGGGAAAGCAAGCCCTGGAAACCG",
			},
			[]string{
				"CGGTAAGCAGGCGCTGGAAACAGTACAG",
				"TGTTTCCCTCCGTCATGCGACGCAATCG",
			},
			false,
			490,
			804,
		},
		{
			"embed additional sequence between fragments",
			Frag{
				ID:    "embedded_primer_seq",
				start: 50,
				end:   350,
				conf:  c,
				db:    db,
			},
			args{
				last: &Frag{
					ID:    "last",
					start: 0,
					end:   45,
					conf:  c,
				},
				next: &Frag{
					ID:    "next",
					start: 355,
					end:   400,
					conf:  c,
				},
				Seq: "GTAAATCCTGGGATCATTCAGTAGTAACCACAAACTTACGCTGGGGCTTCTTTGGCGGATTTTTACAGATACTAACCAGGTGATTTGAAGTAAATTAGTTGAGGATTTAGCCGCGCTATCCGGTAATCTCCAAATTAAAACATACCGTTCCATGAGGGCTAGAATTACTTACCGGCCTTCACCATGCCTGCGCTATACGCGCCCACTCTCCCGTTTATCCGTCCAAGCGGATGCAATGCGATCCTCCGCTAAGATATTCTTACGTGTAACGTAGCTATGTATTTTACAGAGCTGGCGTACGCGTTGAACACTTCACAGATGATAGGGATTCGGGTAAAGAGCGTGTTATTGGGGACTTACACAGGCGTAGACTACAATGGGCCCAACTCAATCACAGCTC",
			},
			[]string{
				"TTACGCTGGGGCTTCTTTGGCGGATTTTTACAGATACT",
				"CCTGTGTAAGTCCCCAATAACACGCTCTTTACCCGA", // rev comp is TCGGGTAAAGAGCGTGTTATTGGGGGACTTACACAGGC
			},
			false,
			50,
			349,
		},
		{
			"optimize when synthesizing neighbors",
			Frag{
				ID:    "optimize_synth",
				start: 125,
				end:   700,
				conf:  c,
				db:    db,
			},
			args{
				last: &Frag{
					ID:    "last",
					start: 0,
					end:   85,
					conf:  c,
				},
				next: &Frag{
					ID:    "next",
					start: 820,
					end:   900,
					conf:  c,
				},
				Seq: "CaGTCAaTCTTTCaCAAaTTTTGTaATCCAGAGGTTGATTATCGATAAGCTTGATATCGAATTCATAACttCgTATAGCATACATTATACGAAGTTATTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTAAGGCAAATAACTTCGTATAGGATACTTTATACGAAGTTATGCTAGCGCCACCATGTCTAGACTGGACAAGAGCAAAGTCATAAACTCTGCTCTGGAATTACTCAATGAAGTCGGTATCGAAGGCCTGACGACAAGGAAACTCGCTCAAAAGCTGGGAGTTGAGCAGCCTACCCTGTACTGGCACGTGAAGAACAAGCGGGCCCTGCTCGATGCCCTGGCAATCGAGATGCTGGACAGGCATCATACCCACTTCTGCCCCCTGGAAGGCGAGTCATGGCAAGACTTTCTGCGGAACAACGCCAAGTCATTCCGCTGTGCTCTCCTCTCACATCGCGACGGGGCTAAAGTGCATCTCGGCACCCGCCCAACAGAGAAACAGTACGAAACCCTGGAAAATCAGCTCGCGTTCCTGTGTCAGCAAGGCTTCTCCCTGGAGAACGCACTGTACGCTCTGTCCGCCGTGGGCCACTTTACACTGGGCTGCGTATTGGAGGATCAGGAGCATCAAGTAGCAAAAGAGGAAAGAGAGACACCTACCACCGATTCTATGCCCCCACTTCTGAGACAAGCAATTGAGCTGTTCGACCATCAGGGAGCCGAACCTGCCTTCCTTTTCGGCCTGGAACTAATCATATGTGGCCTGGAGAAACAGCTAAAGTGCGAAAGCGGCGGGCCGGCCGACGCCCTTGACGATTTTGACTTAGACATGCTCCCAGCCGATGCCCTTGACGACTTTGACCTTGATATGCTGCCTGCTGACGCTCTTGACGATTTtGACCTTGACATGCTCCCCGGGGGATCCGGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGAgACGTGGAGGAGAACCCTGGACCTATGAGCGAGCTGATCAAGGAGAACATGCACATGAAGCTGTACATGGAgggc",
			},
			[]string{
				"ACGAAGTTATGCTAGCGCCA", // rev comp is TGATCCTCCAATACGCAGCC
				"TGATCCTCCAATACGCAGCC", // rev comp is GGCTGCGTATTGGAGGATCA
			},
			false,
			172, // deprecation test only, got this from the output and confirmed it made sense
			639,
		},
		{
			"optimize when there are large overlaps with neighbors",
			Frag{
				ID:    "optimize_overlaps",
				start: 50,
				end:   350,
				conf:  c,
				db:    db,
			},
			args{
				last: &Frag{
					ID:    "last",
					start: 0,
					end:   150,
					conf:  c,
				},
				next: &Frag{
					ID:    "next",
					start: 250,
					end:   400,
					conf:  c,
				},
				Seq: "GTAAATCCTGGGATCATTCAGTAGTAACCACAAACTTACGCTGGGGCTTCTTTGGCGGATTTTTACAGATACTAACCAGGTGATTTGAAGTAAATTAGTTGAGGATTTAGCCGCGCTATCCGGTAATCTCCAAATTAAAACATACCGTTCCATGAGGGCTAGAATTACTTACCGGCCTTCACCATGCCTGCGCTATACGCGCCCACTCTCCCGTTTATCCGTCCAAGCGGATGCAATGCGATCCTCCGCTAAGATATTCTTACGTGTAACGTAGCTATGTATTTTACAGAGCTGGCGTACGCGTTGAACACTTCACAGATGATAGGGATTCGGGTAAAGAGCGTGTTATTGGGGACTTACACAGGCGTAGACTACAATGGGCCCAACTCAATCACAGCTC",
			},
			[]string{
				"AGTTGAGGATTTAGCCGCGC",
				"ACGCTCTTTACCCGAATCCC", // rev comp is TCGGGTAAAGAGCGTGTTATTGGGGGACTTACACAGGC
			},
			false,
			96, // deprecation test only, got this from the output and confirmed it made sense
			343,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			err := tt.n.setPrimers(tt.args.last, tt.args.next, tt.args.Seq, c)
			if (err != nil) != tt.wantErr {
				t.Errorf("setPrimers() error = %v, wantErr %v", err, tt.wantErr)
				return
			}

			if !tt.wantErr && len(tt.n.Primers) < 2 {
				t.Errorf("setPrimers() failed to set primers on %s", tt.n.ID)
			}

			if len(tt.wantPrimers) != len(tt.n.Primers) {
				t.Errorf("setPrimers() on %s, got %d primers wanted %d", tt.n.ID, len(tt.n.Primers), len(tt.wantPrimers))
			}

			for _, primer := range tt.n.Primers {
				if primer.Strand && !strings.Contains(strings.ToUpper(tt.args.Seq), primer.Seq) {
					t.Errorf("setPrimers() FWD primer not contained in the args's Seq: %s", primer.Seq)
				}

				if !primer.Strand && !strings.Contains(strings.ToUpper(tt.args.Seq), revComp(primer.Seq)) {
					t.Errorf("setPrimers() REV primer not contained in the args's Seq: %s", revComp(primer.Seq))
				}

				if primer.Strand && primer.Seq != tt.wantPrimers[0] {
					t.Errorf("setPrimers() FWD primer = %s, want %s", primer.Seq, tt.wantPrimers[0])
				}

				if !primer.Strand && primer.Seq != tt.wantPrimers[1] {
					t.Errorf("setPrimers() REV primer = %s, want %s", primer.Seq, tt.wantPrimers[1])
				}
			}

			if !tt.wantErr && (tt.wantStart != tt.n.start || tt.wantEnd != tt.n.end) {
				t.Errorf("wrong range for setPrimers() on %s, want %d-%d, got %d-%d", tt.n.ID, tt.wantStart, tt.wantEnd, tt.n.start, tt.n.end)
			}
		})
	}
}
