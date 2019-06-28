package repp

import (
	"math"
	"os"
	"reflect"
	"testing"

	"github.com/jjtimmons/repp/config"
)

func Test_primer3_shrink(t *testing.T) {
	type fields struct {
		n      *Frag
		last   *Frag
		next   *Frag
		seq    string
		in     *os.File
		out    *os.File
		p3Path string
		p3Conf string
		p3Dir  string
	}
	type args struct {
		last        *Frag
		n           *Frag
		next        *Frag
		maxHomology int
		minLength   int
	}
	tests := []struct {
		name   string
		fields fields
		args   args
		want   *Frag
	}{
		{
			"shrink Frag with an excessive amount of homology",
			fields{}, // not relevant, nothing used from primer3
			args{
				last: &Frag{
					start: 0,
					end:   100,
				},
				n: &Frag{
					Seq:   "GGGGGAACGCTGAAGATCTCTTCTTCTCATGACTGAACTCGCGAGGGTCGTGATGTCGGTTCCTTCAAAGGTTAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTAGCGGCTGGCGTCTTGAATCCTCGGTCCCCCTTGTCTTTCCAGATTAATCCATTTCCCTCATTCACGAGCTTACCAAGTCAACATTGGTATATGAATGCGACCTTGAAGAGGCCGCTTAAAAATGGCAGTGGTTGAT",
					start: 90,
					end:   300,
				},
				next: &Frag{
					start: 250,
					end:   500,
				},
				maxHomology: 10, // much less than normal,
				minLength:   20,
			},
			&Frag{
				Seq:   "GGGGGAACGCTGAAGATCTCTTCTTCTCATGACTGAACTCGCGAGGGTCGTGATGTCGGTTCCTTCAAAGGTTAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTAGCGGCTGGCGTCTTGAATCCTCGGTCCCCCTTGTCTTTCCAGATTAATCCATTTCCCTCATTCACGAGCTTACCAAGTCAACATTGGTATATGAAT",
				start: 90,
				end:   260,
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			p := &primer3{
				f:              tt.fields.n,
				last:           tt.fields.last,
				next:           tt.fields.next,
				seq:            tt.fields.seq,
				in:             tt.fields.in,
				out:            tt.fields.out,
				primer3Path:    tt.fields.p3Path,
				primer3ConfDir: tt.fields.p3Conf,
			}
			if got := p.shrink(tt.args.last, tt.args.n, tt.args.next, tt.args.maxHomology, tt.args.minLength); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("primer3.shrink() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_bpToAdd(t *testing.T) {
	c := config.New()
	c.PCRMaxEmbedLength = 20
	c.FragmentsMinHomology = 10

	p := primer3{}

	type args struct {
		left  *Frag
		right *Frag
	}
	tests := []struct {
		name        string
		args        args
		wantBpToAdd int
	}{
		{
			"no added homology is needed",
			args{
				left: &Frag{
					start: 0,
					end:   20,
					conf:  c,
				},
				right: &Frag{
					start: 10,
					end:   30,
					conf:  c,
				},
			},
			0,
		},
		{
			"no added homology is needed - lots of overlap",
			args{
				left: &Frag{
					start: 0,
					end:   50,
					conf:  c,
				},
				right: &Frag{
					start: 10,
					end:   30,
					conf:  c,
				},
			},
			0,
		},
		{
			"add homology to each Frag",
			args{
				left: &Frag{
					start: 0,
					end:   10,
					conf:  c,
				},
				right: &Frag{
					start: 16,
					end:   30,
					conf:  c,
				},
			},
			12,
		},
		{
			"correct bp to share when negative",
			args{
				left: &Frag{
					start: -50,
					end:   20,
					conf:  c,
				},
				right: &Frag{
					start: 0,
					end:   1050,
					conf:  c,
				},
			},
			0,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotBpToAdd := p.bpToAdd(tt.args.left, tt.args.right); gotBpToAdd != tt.wantBpToAdd {
				t.Errorf("bpToAdd() = %v, want %v", gotBpToAdd, tt.wantBpToAdd)
			}
		})
	}
}

func Test_mutatePrimers(t *testing.T) {
	type args struct {
		n        *Frag
		Seq      string
		addLeft  int
		addRight int
	}
	tests := []struct {
		name        string
		args        args
		wantMutated *Frag
	}{
		{
			"add homology to both sides of a Frag",
			args{
				n: &Frag{
					Seq:   "TGACCTCGGCTCCCCATTGCTACTACGGCG",
					start: 10,
					end:   39,
					Primers: []Primer{
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
				Seq:      "GATCACTCGATGACCTCGGCTCCCCATTGCTACTACGGCGATTCTTGGAG",
				addLeft:  5,
				addRight: 6,
			},
			&Frag{
				PCRSeq: "CTCGATGACCTCGGCTCCCCATTGCTACTACGGCGATTCTT",
				Seq:    "TGACCTCGGCTCCCCATTGCTACTACGGCG",
				start:  10,
				end:    39,
				Primers: []Primer{
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
				n: &Frag{
					Seq:   "TGACCTCGGCTCCCCATTGCTACTACGGCG",
					start: 10,
					end:   39,
					Primers: []Primer{
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
				Seq:      "GATCACTCGATGACCTCGGCTCCCCATTGCTACTACGGCGATTCTTGGAG",
				addLeft:  0,
				addRight: 0,
			},
			&Frag{
				Seq:    "TGACCTCGGCTCCCCATTGCTACTACGGCG",
				PCRSeq: "TGACCTCGGCTCCCCATTGCTACTACGGCG",
				start:  10,
				end:    39,
				Primers: []Primer{
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
			if gotMutated := mutatePrimers(tt.args.n, tt.args.Seq, tt.args.addLeft, tt.args.addRight); !reflect.DeepEqual(gotMutated, tt.wantMutated) {
				t.Errorf("mutateNodePrimers() = %v, want %v", gotMutated, tt.wantMutated)
			}
		})
	}
}

func Test_reverseComplement(t *testing.T) {
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
		{
			"correct reverse complement of enzyme recog seqs",
			args{
				seq: "ATG^_CAT",
			},
			"ATG^_CAT",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := reverseComplement(tt.args.seq); got != tt.want {
				t.Errorf("revComp() = %v, want %v", got, tt.want)
			}
		})
	}
}

// these estimated hairpin tms jump around when the primer3 version changes
func Test_hairpin(t *testing.T) {
	c := config.New()

	type args struct {
		seq  string
		conf *config.Config
	}
	tests := []struct {
		name     string
		args     args
		wantMelt float64
	}{
		{
			"find hairpin of ~75 degrees",
			args{
				"TGTGCACTCATCATCATCATCGGGGGGGGGGGGTGAACACTATCCCCCCCCCCCCCCA",
				c,
			},
			75.0,
		},
		{
			"return 0 when no hairpin found",
			args{
				"TGTGcactcatcatcCCCA",
				c,
			},
			0.0,
		},
		{
			"return the right-most hairpin when >60bp",
			args{
				"TGTGcactcatcatcaacacaactacgtcgatcagctacgatcgatcgatgctgatcgatatttatatcgagctagctacggatcatcGGGGGGGGGGGGTGAACACTATCCCCCCCCCCCCCCA",
				c,
			},
			75.0,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotMelt := hairpin(tt.args.seq, tt.args.conf); math.Abs(gotMelt-tt.wantMelt) > 1 {
				t.Errorf("hairpin() = %v, want %v", gotMelt, tt.wantMelt)
			}
		})
	}
}
