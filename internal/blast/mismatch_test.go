package blast

import (
	"path"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/jjtimmons/decvec/internal/dvec"
)

func Test_isMismatch(t *testing.T) {
	type args struct {
		match dvec.Match
	}
	tests := []struct {
		name string
		args args
		want bool
	}{
		{
			"catches a mismatching primer",
			args{
				match: dvec.Match{
					Seq:      "atgacgacgacgcggac",
					Mismatch: 0,
				},
			},
			true,
		},
		{
			"no false positive mismatch",
			args{
				match: dvec.Match{
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
	conf.DB = testDB

	type args struct {
		primer string
		parent string
	}
	tests := []struct {
		name         string
		args         args
		wantMismatch bool
		wantMatch    dvec.Match
		wantErr      bool
	}{
		{
			"finds mismatch",
			args{
				"AGTATAGGATAGGTAGTCATTCTT",
				"gnl|addgene|107006",
			},
			true,
			dvec.Match{
				Entry:    "addgene:107006",
				Seq:      "AGTATAGGATAGGTAGTCATTCTT",
				Start:    0,
				End:      23,
				Circular: false,
				Mismatch: 0,
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotMismatch, gotMatch, err := Mismatch(tt.args.primer, tt.args.parent)
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
