package defrag

import (
	"testing"
)

func Test_solutionCost(t *testing.T) {

	synthCostMock := func(length int) float64 {
		return 0.1 * float64(length)
	}

	type args struct {
		frags         []Fragment
		primerBP      float64
		synthCostFunc func(int) float64
	}
	tests := []struct {
		name     string
		args     args
		wantCost float64
	}{
		{
			"pcr only",
			args{
				[]Fragment{
					Fragment{
						Type: PCR,
						Seq:  "agtgcatgcatgcatgctagctagctagctagctacg",
						Primers: []Primer{
							Primer{
								Seq: "atgcatgctgac",
							},
							Primer{
								Seq: "gactgatcgatct",
							},
						},
					},
				},
				0.01,
				synthCostMock,
			},
			25 * 0.01, // 25 total primer bps
		},
		{
			"synth only",
			args{
				[]Fragment{
					Fragment{
						Type: Synthetic,
						Seq:  "agtgcatgcatgcatgctagctagctagctagctacg",
					},
				},
				0.01,
				synthCostMock,
			},
			37 * 0.1, // 37 total synthetic bps * 0.1
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotCost := solutionCost(tt.args.frags, tt.args.primerBP, tt.args.synthCostFunc); gotCost != tt.wantCost {
				t.Errorf("solutionCost() = %v, want %v", gotCost, tt.wantCost)
			}
		})
	}
}
