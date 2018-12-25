package defrag

import (
	"path"
	"reflect"
	"strings"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

// this is little more than a deprecation test right now
func Test_setPrimers(t *testing.T) {
	c := config.New()

	c.Fragments.MinHomology = 20
	c.PCR.P3MaxPenalty = 50.0
	c.PCR.MaxEmbedLength = 10
	db := path.Join("..", "..", "test", "db", "db")

	type args struct {
		last *node
		next *node
		seq  string
	}
	tests := []struct {
		name        string
		n           node
		args        args
		wantPrimers []string // two primers seqs, first is fwd, second is rev
		wantErr     bool
		wantStart   int
		wantEnd     int
	}{
		{
			"create primers without added homology",
			node{
				id:    "gnl|addgene|85039.2",
				start: 0,
				end:   1060,
				conf:  &c,
				db:    db,
			},
			args{
				last: &node{ // close enough that no homology should be added
					id:    "last|fragment",
					start: -50,
					end:   25,
					conf:  &c,
				},
				next: &node{
					id:    "next|fragment",
					start: 1035,
					end:   1090,
					conf:  &c,
				},
				seq: "CaGTCAaTCTTTCaCAAaTTTTGTaATCCAGAGGTTGATTATCGATAAGCTTGATATCGAATTCATAACttCgTATAGCATACATTATACGAAGTTATTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTAAGGCAAATAACTTCGTATAGGATACTTTATACGAAGTTATGCTAGCGCCACCATGTCTAGACTGGACAAGAGCAAAGTCATAAACTCTGCTCTGGAATTACTCAATGAAGTCGGTATCGAAGGCCTGACGACAAGGAAACTCGCTCAAAAGCTGGGAGTTGAGCAGCCTACCCTGTACTGGCACGTGAAGAACAAGCGGGCCCTGCTCGATGCCCTGGCAATCGAGATGCTGGACAGGCATCATACCCACTTCTGCCCCCTGGAAGGCGAGTCATGGCAAGACTTTCTGCGGAACAACGCCAAGTCATTCCGCTGTGCTCTCCTCTCACATCGCGACGGGGCTAAAGTGCATCTCGGCACCCGCCCAACAGAGAAACAGTACGAAACCCTGGAAAATCAGCTCGCGTTCCTGTGTCAGCAAGGCTTCTCCCTGGAGAACGCACTGTACGCTCTGTCCGCCGTGGGCCACTTTACACTGGGCTGCGTATTGGAGGATCAGGAGCATCAAGTAGCAAAAGAGGAAAGAGAGACACCTACCACCGATTCTATGCCCCCACTTCTGAGACAAGCAATTGAGCTGTTCGACCATCAGGGAGCCGAACCTGCCTTCCTTTTCGGCCTGGAACTAATCATATGTGGCCTGGAGAAACAGCTAAAGTGCGAAAGCGGCGGGCCGGCCGACGCCCTTGACGATTTTGACTTAGACATGCTCCCAGCCGATGCCCTTGACGACTTTGACCTTGATATGCTGCCTGCTGACGCTCTTGACGATTTtGACCTTGACATGCTCCCCGGGGGATCCGGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGAgACGTGGAGGAGAACCCTGGACCTATGAGCGAGCTGATCAAGGAGAACATGCACATGAAGCTGTACATGGAgggc",
			},
			[]string{
				"CAGTCAATCTTTCACAAATTTTG",
				"GCCCTCCATGTACAGCTTCA", // rev-comp: TGAAGCTGTACATGGAGGGC
			},
			false,
			0,
			1060,
		},
		{
			"create primers with homology",
			node{
				id:    "add_homology",
				start: 500,
				end:   800,
				conf:  &c,
				db:    db,
			},
			args{
				last: &node{
					id:    "last",
					start: 200,
					end:   500,
					conf:  &c,
				},
				next: &node{
					id:    "next",
					start: 795,
					end:   900,
					conf:  &c,
				},
				seq: "TGCTGACTGTGGCGGGTGAGCTTAGGGGGCCTCCGCTCCAGCTCGACACCGGGCAGCTGCTGAAGATCGCGAAGAGAGGGGGAGTAACAGCGGTAGAGGCAGTGCACGCCTGGCGCAATGCGCTCACCGGGGCCCCCTTGAACCTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAGCAATGGGGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCCACGGGCTGACTCCCGATCAAGTTGTAGCGATTGCGTCGCATGACGGAGGGAAACAAGCATTGGAGACTGTCCAACGGCTCCTTCCCGTGTTGTGTCAAGCCCACGGTTTGACGCCTGCACAAGTGGTCGCCATCGCCAGCCATGATGGCGGTAAGCAGGCGCTGGAAACAGTACAGCGCCTGCTGCCTGTACTGTGCCAGGATCATGGACTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAATAATAACGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCCACGGGCTGACTCCCGATCAAGTTGTAGCGATTGCGTCGCATGACGGAGGGAAACAAGCATTGGAGACTGTCCAACGGCTCCTTCCCGTGTTGTGTCAAGCCCACGGTTTGACGCCTGCACAAGTGGTCGCCATCGCCAACAACAACGGCGGTAAGCAGGCGCTGGAAACAGTACAGCGCCTGCTGCCTGTACTGTGCCAGGATCATGGACTGACCCCAGACCAGGTAGTCGCAATCGCGTCGAACATTGGGGGAAAGCAAGCCCTGGAAACCG",
			},
			[]string{
				"CGGTAAGCAGGCGCTGGAAACAGTACAG",
				"TTGTTTCCCTCCGTCATGCGACGCAATC",
			},
			false,
			490,
			805,
		},
		{
			"fail to create primers with off-targets",
			node{
				id:    "gnl|addgene|39412.1",
				start: 0,
				end:   1022,
				conf:  &c,
				db:    db,
			},
			args{
				last: &node{
					id:    "last",
					start: -50,
					end:   25,
					conf:  &c,
				},
				next: &node{
					id:    "next",
					start: 1000,
					end:   1050,
					conf:  &c,
				},
				seq: "TGCTGACTGTGGCGGGTGAGCTTAGGGGGCCTCCGCTCCAGCTCGACACCGGGCAGCTGCTGAAGATCGCGAAGAGAGGGGGAGTAACAGCGGTAGAGGCAGTGCACGCCTGGCGCAATGCGCTCACCGGGGCCCCCTTGAACCTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAGCAATGGGGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCCACGGGCTGACTCCCGATCAAGTTGTAGCGATTGCGTCGCATGACGGAGGGAAACAAGCATTGGAGACTGTCCAACGGCTCCTTCCCGTGTTGTGTCAAGCCCACGGTTTGACGCCTGCACAAGTGGTCGCCATCGCCAGCCATGATGGCGGTAAGCAGGCGCTGGAAACAGTACAGCGCCTGCTGCCTGTACTGTGCCAGGATCATGGACTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAATAATAACGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCCACGGGCTGACTCCCGATCAAGTTGTAGCGATTGCGTCGCATGACGGAGGGAAACAAGCATTGGAGACTGTCCAACGGCTCCTTCCCGTGTTGTGTCAAGCCCACGGTTTGACGCCTGCACAAGTGGTCGCCATCGCCAACAACAACGGCGGTAAGCAGGCGCTGGAAACAGTACAGCGCCTGCTGCCTGTACTGTGCCAGGATCATGGACTGACCCCAGACCAGGTAGTCGCAATCGCGTCGAACATTGGGGGAAAGCAAGCCCTGGAAACCG",
			},
			nil,
			true,
			0,
			1022,
		},
		{
			"primers with extra space to avoid bad scores",
			node{
				id:    "add_extra_space",
				start: 125,
				end:   700,
				conf:  &c,
				db:    db,
			},
			args{
				last: &node{
					id:    "last",
					start: 0,
					end:   85,
					conf:  &c,
				},
				next: &node{
					id:    "next",
					start: 820,
					end:   900,
					conf:  &c,
				},
				seq: "CaGTCAaTCTTTCaCAAaTTTTGTaATCCAGAGGTTGATTATCGATAAGCTTGATATCGAATTCATAACttCgTATAGCATACATTATACGAAGTTATTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTAAGGCAAATAACTTCGTATAGGATACTTTATACGAAGTTATGCTAGCGCCACCATGTCTAGACTGGACAAGAGCAAAGTCATAAACTCTGCTCTGGAATTACTCAATGAAGTCGGTATCGAAGGCCTGACGACAAGGAAACTCGCTCAAAAGCTGGGAGTTGAGCAGCCTACCCTGTACTGGCACGTGAAGAACAAGCGGGCCCTGCTCGATGCCCTGGCAATCGAGATGCTGGACAGGCATCATACCCACTTCTGCCCCCTGGAAGGCGAGTCATGGCAAGACTTTCTGCGGAACAACGCCAAGTCATTCCGCTGTGCTCTCCTCTCACATCGCGACGGGGCTAAAGTGCATCTCGGCACCCGCCCAACAGAGAAACAGTACGAAACCCTGGAAAATCAGCTCGCGTTCCTGTGTCAGCAAGGCTTCTCCCTGGAGAACGCACTGTACGCTCTGTCCGCCGTGGGCCACTTTACACTGGGCTGCGTATTGGAGGATCAGGAGCATCAAGTAGCAAAAGAGGAAAGAGAGACACCTACCACCGATTCTATGCCCCCACTTCTGAGACAAGCAATTGAGCTGTTCGACCATCAGGGAGCCGAACCTGCCTTCCTTTTCGGCCTGGAACTAATCATATGTGGCCTGGAGAAACAGCTAAAGTGCGAAAGCGGCGGGCCGGCCGACGCCCTTGACGATTTTGACTTAGACATGCTCCCAGCCGATGCCCTTGACGACTTTGACCTTGATATGCTGCCTGCTGACGCTCTTGACGATTTtGACCTTGACATGCTCCCCGGGGGATCCGGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGAgACGTGGAGGAGAACCCTGGACCTATGAGCGAGCTGATCAAGGAGAACATGCACATGAAGCTGTACATGGAgggc",
			},
			[]string{
				"ACGAAGTTATGCTAGCGCCA", // rev comp is TGATCCTCCAATACGCAGCC
				"TGATCCTCCAATACGCAGCC", // rev comp is GGCTGCGTATTGGAGGATCA
			},
			false,
			172, // deprecation test only, got this from the output
			639,
		},
		{
			"primers that embed additional sequence between fragments",
			node{
				id:    "embedded_primer_seq",
				start: 50,
				end:   350,
				conf:  &c,
				db:    db,
			},
			args{
				last: &node{
					id:    "last",
					start: 0,
					end:   45,
					conf:  &c,
				},
				next: &node{
					id:    "next",
					start: 355,
					end:   400,
					conf:  &c,
				},
				seq: "GTAAATCCTGGGATCATTCAGTAGTAACCACAAACTTACGCTGGGGCTTCTTTGGCGGATTTTTACAGATACTAACCAGGTGATTTGAAGTAAATTAGTTGAGGATTTAGCCGCGCTATCCGGTAATCTCCAAATTAAAACATACCGTTCCATGAGGGCTAGAATTACTTACCGGCCTTCACCATGCCTGCGCTATACGCGCCCACTCTCCCGTTTATCCGTCCAAGCGGATGCAATGCGATCCTCCGCTAAGATATTCTTACGTGTAACGTAGCTATGTATTTTACAGAGCTGGCGTACGCGTTGAACACTTCACAGATGATAGGGATTCGGGTAAAGAGCGTGTTATTGGGGACTTACACAGGCGTAGACTACAATGGGCCCAACTCAATCACAGCTC",
			},
			[]string{
				"TTACGCTGGGGCTTCTTTGGCGGATTTTTACAGATACT",
				"GCCTGTGTAAGTCCCCAATAACACGCTCTTTACCCGA", // rev comp is TCGGGTAAAGAGCGTGTTATTGGGGGACTTACACAGGC
			},
			false,
			35, // deprecation test only, got this from the output
			365,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			err := tt.n.setPrimers(tt.args.last, tt.args.next, tt.args.seq, &c)
			if (err != nil) != tt.wantErr {
				t.Errorf("setPrimers() error = %v, wantErr %v", err, tt.wantErr)
				return
			}

			if !tt.wantErr && len(tt.n.primers) < 2 {
				t.Errorf("setPrimers() failed to set primers on %s", tt.n.id)
			}

			if len(tt.wantPrimers) != len(tt.n.primers) {
				t.Errorf("setPrimers() on %s, got %d primers wanted %d", tt.n.id, len(tt.n.primers), len(tt.wantPrimers))
			}

			for _, primer := range tt.n.primers {
				if primer.Strand && !strings.Contains(strings.ToUpper(tt.n.seq), primer.Seq) {
					t.Errorf("setPrimers() FWD primer not contained in the node's seq: %s", primer.Seq)
				}

				if !primer.Strand && !strings.Contains(strings.ToUpper(tt.n.seq), revComp(primer.Seq)) {
					t.Errorf("setPrimers() REV primer not contained in the node's seq: %s", revComp(primer.Seq))
				}

				if primer.Strand && primer.Seq != tt.wantPrimers[0] {
					t.Errorf("setPrimers() FWD primer = %s, want %s", primer.Seq, tt.wantPrimers[0])
				}

				if !primer.Strand && primer.Seq != tt.wantPrimers[1] {
					t.Errorf("setPrimers() REV primer = %s, want %s", primer.Seq, tt.wantPrimers[1])
				}
			}

			if !tt.wantErr && (tt.wantStart != tt.n.start || tt.wantEnd != tt.n.end) {
				t.Errorf("wrong range for setPrimers() on %s, want %d-%d, got %d-%d", tt.n.id, tt.wantStart, tt.wantEnd, tt.n.start, tt.n.end)
			}
		})
	}
}

func Test_p3Exec_input(t *testing.T) {
	type fields struct {
		n      *node
		last   *node
		next   *node
		seq    string
		in     string
		out    string
		p3Path string
		p3Conf string
		p3Dir  string
	}
	type args struct {
		minHomology    int
		maxEmbedLength int
	}
	tests := []struct {
		name           string
		fields         fields
		args           args
		wantBpAddLeft  int
		wantBpAddRight int
		wantErr        bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			p := &p3Exec{
				n:      tt.fields.n,
				last:   tt.fields.last,
				next:   tt.fields.next,
				seq:    tt.fields.seq,
				in:     tt.fields.in,
				out:    tt.fields.out,
				p3Path: tt.fields.p3Path,
				p3Conf: tt.fields.p3Conf,
				p3Dir:  tt.fields.p3Dir,
			}
			gotBpAddLeft, gotBpAddRight, err := p.input(tt.args.minHomology, tt.args.maxEmbedLength)
			if (err != nil) != tt.wantErr {
				t.Errorf("p3Exec.input() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if gotBpAddLeft != tt.wantBpAddLeft {
				t.Errorf("p3Exec.input() gotBpAddLeft = %v, want %v", gotBpAddLeft, tt.wantBpAddLeft)
			}
			if gotBpAddRight != tt.wantBpAddRight {
				t.Errorf("p3Exec.input() gotBpAddRight = %v, want %v", gotBpAddRight, tt.wantBpAddRight)
			}
		})
	}
}

func Test_bpToShare(t *testing.T) {
	c := config.New()

	type args struct {
		left        *node
		right       *node
		minHomology int
	}
	tests := []struct {
		name        string
		args        args
		wantBpToAdd int
	}{
		{
			"no added homology is needed",
			args{
				left: &node{
					start: 0,
					end:   20,
					conf:  &c,
				},
				right: &node{
					start: 10,
					end:   30,
					conf:  &c,
				},
				minHomology: 10,
			},
			0,
		},
		{
			"no added homology is needed - lots of overlap",
			args{
				left: &node{
					start: 0,
					end:   50,
					conf:  &c,
				},
				right: &node{
					start: 10,
					end:   30,
					conf:  &c,
				},
				minHomology: 10,
			},
			0,
		},
		{
			"add homology to each node",
			args{
				left: &node{
					start: 0,
					end:   10,
					conf:  &c,
				},
				right: &node{
					start: 16,
					end:   30,
					conf:  &c,
				},
				minHomology: 3,
			},
			7,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotBpToAdd := bpToShare(tt.args.left, tt.args.right, tt.args.minHomology); gotBpToAdd != tt.wantBpToAdd {
				t.Errorf("bpToShare() = %v, want %v", gotBpToAdd, tt.wantBpToAdd)
			}
		})
	}
}

func Test_mutateNodePrimers(t *testing.T) {
	type args struct {
		n        *node
		seq      string
		addLeft  int
		addRight int
	}
	tests := []struct {
		name        string
		args        args
		wantMutated *node
	}{
		{
			"add homology to both sides of a node",
			args{
				n: &node{
					seq:   "TGACCTCGGCTCCCCATTGCTACTACGGCG",
					start: 10,
					end:   29,
					primers: []Primer{
						Primer{
							Seq: "TGACCTCGGC",
							Range: ranged{
								start: 10,
								end:   20,
							},
							Strand: true,
						},
						Primer{
							Seq: "CGCCGTAGTA", // rev comp TACTACGGCG
							Range: ranged{
								start: 30,
								end:   39,
							},
							Strand: false,
						},
					},
				},
				seq:      "GATCACTCGATGACCTCGGCTCCCCATTGCTACTACGGCGATTCTTGGAG",
				addLeft:  5,
				addRight: 6,
			},
			&node{
				seq:   "CTCGATGACCTCGGCTCCCCATTGCTACTACGGCGATTCTT",
				start: 5,
				end:   45,
				primers: []Primer{
					Primer{
						Seq: "CTCGATGACCTCGGC",
						Range: ranged{
							start: 5,
							end:   20,
						},
						Strand: true,
					},
					Primer{
						Seq: "AAGAATCGCCGTAGTA", //  rev comp TACTACGGCGATTCTT
						Range: ranged{
							start: 30,
							end:   45,
						},
						Strand: false,
					},
				},
			},
		},
		{
			"add nothing",
			args{
				n: &node{
					seq:   "TGACCTCGGCTCCCCATTGCTACTACGGCG",
					start: 10,
					end:   39,
					primers: []Primer{
						Primer{
							Seq: "TGACCTCGGC",
							Range: ranged{
								start: 10,
								end:   20,
							},
							Strand: true,
						},
						Primer{
							Seq: "CGCCGTAGTA", // rev comp TACTACGGCG
							Range: ranged{
								start: 30,
								end:   39,
							},
							Strand: false,
						},
					},
				},
				seq:      "GATCACTCGATGACCTCGGCTCCCCATTGCTACTACGGCGATTCTTGGAG",
				addLeft:  0,
				addRight: 0,
			},
			&node{
				seq:   "TGACCTCGGCTCCCCATTGCTACTACGGCG",
				start: 10,
				end:   39,
				primers: []Primer{
					Primer{
						Seq: "TGACCTCGGC",
						Range: ranged{
							start: 10,
							end:   20,
						},
						Strand: true,
					},
					Primer{
						Seq: "CGCCGTAGTA", // rev comp TACTACGGCG
						Range: ranged{
							start: 30,
							end:   39,
						},
						Strand: false,
					},
				},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotMutated := mutateNodePrimers(tt.args.n, tt.args.seq, tt.args.addLeft, tt.args.addRight); !reflect.DeepEqual(gotMutated, tt.wantMutated) {
				t.Errorf("mutateNodePrimers() = %v, want %v", gotMutated, tt.wantMutated)
			}
		})
	}
}

func Test_revComp(t *testing.T) {
	type args struct {
		seq string
	}
	tests := []struct {
		name string
		args args
		want string
	}{
		{
			"generates reverse complment",
			args{
				seq: "ATGtgca",
			},
			"TGCACAT",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := revComp(tt.args.seq); got != tt.want {
				t.Errorf("revComp() = %v, want %v", got, tt.want)
			}
		})
	}
}
