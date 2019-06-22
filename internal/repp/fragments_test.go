package repp

import (
	"testing"
)

func Test_annealFragments(t *testing.T) {
	type args struct {
		min   int
		max   int
		frags []*Frag
	}
	tests := []struct {
		name      string
		args      args
		wantFrags []*Frag
		wantVec   string
	}{
		{
			"don't change two fragments without overlap",
			args{
				min: 5,
				max: 10,
				frags: []*Frag{
					&Frag{
						Seq: "GGCTAATATAGCGAATTGCCGAGAACCCGGCCCCACGCAATGGAACGTCTTTAGCTCCGGCAGGCAATTAAGGACAACGTAAGTATAGCGCATATAAACA",
					},
					&Frag{
						Seq: "GAGAAATGGGCGAATGAACCTATTCGTACCGTATCGAAGAATAGCCTCGCGGAGGCATGTGCCATGCTAGCGTGCGGGGCACTCTAGTTATGCATATGGT",
					},
				},
			},
			[]*Frag{
				&Frag{
					start: 0,
					end:   99,
				},
				&Frag{
					start: 100,
					end:   199,
				},
			},
			"GGCTAATATAGCGAATTGCCGAGAACCCGGCCCCACGCAATGGAACGTCTTTAGCTCCGGCAGGCAATTAAGGACAACGTAAGTATAGCGCATATAAACAGAGAAATGGGCGAATGAACCTATTCGTACCGTATCGAAGAATAGCCTCGCGGAGGCATGTGCCATGCTAGCGTGCGGGGCACTCTAGTTATGCATATGGT",
		},
		{
			"change the range of two fragments with overlap on both ends",
			args{
				min: 5,
				max: 10,
				frags: []*Frag{
					&Frag{
						Seq: "TGCATATGGTGCGAATTGCCGAGAACCCGGCCCCACGCAATGGAACGTCTTTAGCTCCGGCAGGCAATTAAGGACAACGTAAGTATAGCGCATATAAACA",
					},
					&Frag{
						Seq: "CATATAAACACGAATGAACCTATTCGTACCGTATCGAAGAATAGCCTCGCGGAGGCATGTGCCATGCTAGCGTGCGGGGCACTCTAGTTATGCATATGGT",
					},
				},
			},
			[]*Frag{
				&Frag{
					start: 0,
					end:   99,
				},
				&Frag{
					start: 90,
					end:   189,
				},
			},
			"TGCATATGGTGCGAATTGCCGAGAACCCGGCCCCACGCAATGGAACGTCTTTAGCTCCGGCAGGCAATTAAGGACAACGTAAGTATAGCGCATATAAACACGAATGAACCTATTCGTACCGTATCGAAGAATAGCCTCGCGGAGGCATGTGCCATGCTAGCGTGCGGGGCACTCTAGTTA",
		},
		{
			"change three fragments all annealing",
			args{
				min: 5,
				max: 15,
				frags: []*Frag{
					&Frag{
						Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
					},
					&Frag{
						Seq: "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
					},
					&Frag{
						Seq: "TCGACTAGCTAGAACTGATCTAGACGTGCTAGCTACA",
					},
				},
			},
			[]*Frag{
				&Frag{
					start: 0,
					end:   36,
				},
				&Frag{
					start: 26,
					end:   62,
				},
				&Frag{
					start: 51,
					end:   87,
				},
			},
			"ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCGACTGATCACTAGCATCGACTAGCTAGAACTGATCTAG",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotVec := annealFragments(tt.args.min, tt.args.max, tt.args.frags)

			if gotVec != tt.wantVec {
				t.Errorf("annealFragments() = %v, want %v", gotVec, tt.wantVec)
			}

			for i, f := range tt.wantFrags {
				if f.start != tt.args.frags[i].start {
					t.Errorf("annealFragments() = %v, want %v", tt.args.frags[i].start, f.start)
				}

				if f.end != tt.args.frags[i].end {
					t.Errorf("annealFragments() = %v, want %v", tt.args.frags[i].end, f.end)
				}
			}

		})
	}
}
