package defrag

import (
	"reflect"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

func Test_assembleFragments(t *testing.T) {
	c := config.New()
	c.Fragments.MinHomology = 8
	c.Fragments.MaxHomology = 20

	type args struct {
		inputFragments []Frag
		conf           *config.Config
	}
	tests := []struct {
		name             string
		args             args
		wantTargetVector Frag
		wantFragments    []Frag
	}{
		{
			"fragments with enough overlap",
			args{
				[]Frag{
					Frag{
						Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
					},
					Frag{
						Seq: "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
					},
					Frag{
						Seq: "TCGACTAGCTAGAACTGATCTAGACGTGCTAGCTACA",
					},
				},
				c,
			},
			Frag{
				Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCGACTGATCACTAGCATCGACTAGCTAGAACTGATCTAG",
			},
			[]Frag{
				Frag{
					Seq:   "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
					Type:  existing,
					start: -14,
					end:   22,
				},
				Frag{
					Seq:  "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
					Type: existing,
				},
				Frag{
					Seq:  "TCGACTAGCTAGAACTGATCTAGACGTGCTAGCTACA",
					Type: existing,
				},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotTargetVector, gotFragments := assembleFragments(tt.args.inputFragments, tt.args.conf)

			if !reflect.DeepEqual(gotTargetVector, tt.wantTargetVector) {
				t.Errorf("assembleFWD() gotTargetVector = %v, want %v", gotTargetVector, tt.wantTargetVector)
			}

			for i, wantF := range tt.wantFragments {
				if wantF.Seq != gotFragments[i].Seq {
					t.Errorf("assembleFWD() gotFragment.Seq = %v, want %v", gotFragments[i].Seq, wantF.Seq)
				}

				if wantF.Type != gotFragments[i].Type {
					t.Errorf("assembleFWD() gotFragment.Type = %v, want %v", gotFragments[i].Type, wantF.Type)
				}
			}
		})
	}
}
