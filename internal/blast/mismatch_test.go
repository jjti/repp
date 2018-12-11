package blast

import (
	"path"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/jjtimmons/defrag/config"
	"github.com/jjtimmons/defrag/internal/defrag"
)

func Test_isMismatch(t *testing.T) {
	type args struct {
		match defrag.Match
	}
	tests := []struct {
		name string
		args args
		want bool
	}{
		{
			"find a mismatching primer",
			args{
				match: defrag.Match{
					Seq:      "atgacgacgacgcggac",
					Mismatch: 0,
				},
			},
			true,
		},
		{
			"no false positive mistmatch",
			args{
				match: defrag.Match{
					Seq:      "atgacgacgacgac",
					Mismatch: 0,
				},
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := isMismatch(tt.args.match); got != tt.want {
				t.Errorf("isMismatch() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestMismatch(t *testing.T) {
	testDB, _ := filepath.Abs(path.Join("..", "..", "test", "blast", "db"))

	vendors := config.New().Vendors()

	type args struct {
		primer string
		parent string
	}
	tests := []struct {
		name         string
		args         args
		wantMismatch bool
		wantMatch    defrag.Match
		wantErr      bool
	}{
		{
			"avoids false positive",
			args{
				"GTTGGAGTCCACGTTCTTT",
				"gnl|addgene|113726(circular)",
			},
			false,
			defrag.Match{},
			false,
		},
		// I intentionally added another off-target seq to 107006, AGTATAGTAGGTAGTCATTCTT
		{
			"finds mismatch",
			args{
				"AGTATAGGATAGGTAGTCATTCTT",
				"gnl|addgene|107006(circular)",
			},
			true,
			defrag.Match{
				Entry:    "addgene:107006(circular)",
				Seq:      "AGTATAGTAGGTAGTCATTCTT",
				Start:    0,
				End:      23,
				Circular: true,
				Mismatch: 0,
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotMismatch, gotMatch, err := Mismatch(tt.args.primer, tt.args.parent, testDB, vendors)
			if (err != nil) != tt.wantErr {
				t.Errorf("Mismatch() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if gotMismatch != tt.wantMismatch {
				t.Errorf("Mismatch() gotMismatch = %v, want %v", gotMismatch, tt.wantMismatch)
			}
			if !reflect.DeepEqual(gotMatch, tt.wantMatch) {
				t.Errorf("Mismatch() gotMatch = %v, want %v", gotMatch, tt.wantMatch)
			}
		})
	}
}
