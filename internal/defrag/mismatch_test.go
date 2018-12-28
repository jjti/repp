package defrag

import (
	"path"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

func Test_isMismatch(t *testing.T) {
	c := config.New()

	type args struct {
		match match
	}
	tests := []struct {
		name string
		args args
		want bool
	}{
		{
			"find a mismatching primer",
			args{
				match: match{
					seq:         "atgacgacgacgcggac",
					mismatching: 0,
				},
			},
			true,
		},
		{
			"no false positive mistmatch",
			args{
				match: match{
					seq:         "atgacgacgacgac",
					mismatching: 0,
				},
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := isMismatch(tt.args.match, &c); got != tt.want {
				t.Errorf("isMismatch() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestQueryMismatch(t *testing.T) {
	testDB, _ := filepath.Abs(path.Join("..", "..", "test", "db", "db"))

	conf := config.New()

	type args struct {
		primer string
		parent string
	}
	tests := []struct {
		name         string
		args         args
		wantMismatch bool
		wantMatch    match
		wantErr      bool
	}{
		{
			"avoids false positive",
			args{
				"GTTGGAGTCCACGTTCTTT",
				"gnl|addgene|113726(circular)",
			},
			false,
			match{},
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
			match{
				entry:       "addgene:107006(circular)",
				seq:         "AGTATAGTAGGTAGTCATTCTT",
				start:       0,
				end:         23,
				circular:    true,
				mismatching: 0,
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotMismatch, gotMatch, err := parentMismatch([]Primer{
				Primer{
					Seq: tt.args.primer,
				},
			}, tt.args.parent, testDB, &conf)
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
