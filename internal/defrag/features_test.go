package defrag

import (
	"path/filepath"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

func TestNewFeatureDB(t *testing.T) {
	db := NewFeatureDB()

	if len(db.features) < 1 {
		t.Fail()
	}
}

func Test_features(t *testing.T) {
	test1, conf := NewFlags(
		"SV40 origin, p10 promoter, mEGFP, T7 terminator",
		filepath.Join("..", "..", "test", "output", "features.json"),
		"pSB1A3",
		"EcoRI",
		"",
		[]string{},
		true,
		true,
	)

	type args struct {
		flags *Flags
		conf  *config.Config
	}
	tests := []struct {
		name string
		args args
	}{
		{
			"test end to end features creation",
			args{
				flags: test1,
				conf:  conf,
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			features(tt.args.flags, tt.args.conf)

			t.Fail()
		})
	}
}
