package defrag

import (
	"testing"

	"github.com/jjtimmons/defrag/config"
)

func Test_hairpin(t *testing.T) {
	c := config.New()

	type args struct {
		seq  string
		conf *config.Config
	}
	tests := []struct {
		name     string
		args     args
		wantMelt int
	}{
		{
			"find hairpin of ~86 degrees",
			args{
				"TGTGCACTCATCATCATCATCGGGGGGGGGGGGTGAACACTATCCCCCCCCCCCCCCA",
				c,
			},
			86,
		},
		{
			"return 0 when no hairpin found",
			args{
				"TGTGcactcatcatcCCCA",
				c,
			},
			0,
		},
		{
			"return the right-most hairpin when >60bp",
			args{
				"TGTGcactcatcatcaacacaactacgtcgatcagctacgatcgatcgatgctgatcgatatttatatcgagctagctacggatcatcGGGGGGGGGGGGTGAACACTATCCCCCCCCCCCCCCA",
				c,
			},
			86,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotMelt := hairpin(tt.args.seq, tt.args.conf); gotMelt != tt.wantMelt {
				t.Errorf("hairpin() = %v, want %v", gotMelt, tt.wantMelt)
			}
		})
	}
}
