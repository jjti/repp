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
		inputFragments []Fragment
		conf           *config.Config
	}
	tests := []struct {
		name             string
		args             args
		wantTargetVector Fragment
		wantFragments    []Fragment
	}{
		{
			"fragments with enough overlap",
			args{
				[]Fragment{
					Fragment{
						Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
					},
					Fragment{
						Seq: "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
					},
					Fragment{
						Seq: "TCGACTAGCTAGAACTGATCTAGACGTGCTAGCTACA",
					},
				},
				&c,
			},
			Fragment{
				Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCGACTGATCACTAGCATCGACTAGCTAGAACTGATCTAG",
			},
			[]Fragment{
				Fragment{
					Seq:  "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
					Type: existing,
				},
				Fragment{
					Seq:  "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
					Type: existing,
				},
				Fragment{
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
			if !reflect.DeepEqual(gotFragments, tt.wantFragments) {
				t.Errorf("assembleFWD() gotFragments = %v, want %v", gotFragments, tt.wantFragments)
			}
		})
	}
}
