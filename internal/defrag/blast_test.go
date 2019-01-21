package defrag

import (
	"path"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

// test the ability to find test fragments in a mock database
// see test/blast/README.md for a description of where the subfragments
// in this test fragment's sequence came from (pieces from the 5 fragments)
// that make up the mock BLAST db
func Test_BLAST(t *testing.T) {
	// make path to test db
	testDB, _ := filepath.Abs(path.Join("..", "..", "test", "db", "db"))

	// create mock test fragment
	f := Frag{
		ID:  "test_target",
		Seq: "GGCCGCAATAAAATATCTTTATTTTCATTACATCTGTGTGTTGGTTTTTTGTGTGAATCGATAGTACTAACATGACCACCTTGATCTTCATGGTCTGGGTGCCCTCGTAGGGCTTGCCTTCGCCCTCGGATGTGCACTTGAAGTGGTGGTTGTTCACGGTGCCCTCCATGTACAGCTTCATGTGCATGTTCTCCTTGATCAGCTCGCTCATAGGTCCAGGGTTCTCCTCCACGTCTCCAGCCTGCTTCAGCAGGCTGAAGTTAGTAGCTCCGCTTCCGGATCCCCCGGGGAGCATGTCAAGGTCAAAATCGTCAAGAGCGTCAGCAGGCAGCATATCAAGGTCAAAGTCGTCAAGGGCATCGGCTGGGAgCATGTCTAAgTCAAAATCGTCAAGGGCGTCGGCCGGCCCGCCGCTTTcgcacGCCCTGGCAATCGAGATGCTGGACAGGCATCATACCCACTTCTGCCCCCTGGAAGGCGAGTCATGGCAAGACTTTCTGCGGAACAACGCCAAGTCATTCCGCTGTGCTCTCCTCTCACATCGCGACGGGGCTAAAGTGCATCTCGGCACCCGCCCAACAGAGAAACAGTACGAAACCCTGGAAAATCAGCTCGCGTTCCTGTGTCAGCAAGGCTTCTCCCTGGAGAACGCACTGTACGCTCTGTCCGCCGTGGGCCACTTTACACTGGGCTGCGTATTGGAGGATCAGGAGCATCAAGTAGCAAAAGAGGAAAGAGAGACACCTACCACCGATTCTATGCCTGACTGTGGCGGGTGAGCTTAGGGGGCCTCCGCTCCAGCTCGACACCGGGCAGCTGCTGAAGATCGCGAAGAGAGGGGGAGTAACAGCGGTAGAGGCAGTGCACGCCTGGCGCAATGCGCTCACCGGGGCCCCCTTGAACCTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAGCAATGGGGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTTTAACAAAATATTAACGCTTACAATTTAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAA",
	}

	// run blast
	matches, err := blast(&f, []string{testDB}, 10) // any match over 10 bp

	// check if it fails
	if err != nil {
		t.Errorf("failed to run BLAST: %v", err)
		return
	}

	// make sure matches are found
	if len(matches) < 1 {
		t.Error("failed to find any matches")
		return
	}

	matchesContain := func(targ match) {
		for _, m := range matches {
			if targ.entry == m.entry && targ.start == m.start && targ.end == m.end {
				return
			}
		}

		t.Errorf("failed to find match %v in fragment matches", targ)
	}

	matchesContain(match{
		entry: "gnl|addgene|107006(circular)",
		start: 0,
		end:   72,
	})
}

// test that we can filter out overlapping regions from blast results
// and those that are up against the edge of the fragment
func Test_filter(t *testing.T) {
	// test fragment with 3 matches that should be removed
	matches := []match{
		// shouldn't be removed
		match{
			entry: "m1",
			start: 15,
			end:   19,
		},
		// should be removed because it fits within m3
		match{
			entry: "m2",
			start: 29,
			end:   34,
		},
		// shouldn't be removed
		match{
			entry: "m3",
			start: 29,
			end:   35,
		},
		// shouldn't be removed
		match{
			entry: "m4",
			start: 31,
			end:   72,
		},
	}

	newMatches := filter(matches, 72, 3) // keep all fragments larger than 3bp (all of them)

	// make sure they're gone
	if len(newMatches) != 3 {
		t.Errorf("%d filtered matches found on test fragment, 3 expected: %v", len(newMatches), newMatches)
	}

	// make sure m2 has been removed
	for _, m := range newMatches {
		if m.entry == "m2" {
			t.Error("m2 found in resulting matches, should have been removed")
		}
	}
}

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
			if got := isMismatch(tt.args.match, c); got != tt.want {
				t.Errorf("isMismatch() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_parentMismatch(t *testing.T) {
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
			}, tt.args.parent, testDB, conf)
			if (err != nil) != tt.wantErr {
				t.Errorf("parentMismatch() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if gotMismatch != tt.wantMismatch {
				t.Errorf("parentMismatch() gotMismatch = %v, want %v", gotMismatch, tt.wantMismatch)
			}
			if !reflect.DeepEqual(gotMatch, tt.wantMatch) {
				t.Errorf("parentMismatch() gotMatch = %v, want %v", gotMatch, tt.wantMatch)
			}
		})
	}
}
