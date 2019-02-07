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
	id := "test_target"
	seq := "GGCCGCAATAAAATATCTTTATTTTCATTACATCTGTGTGTTGGTTTTTTGTGTGAATCGATAGTACTAACATGACCACCTTGATCTTCATGGTCTGGGTGCCCTCGTAGGGCTTGCCTTCGCCCTCGGATGTGCACTTGAAGTGGTGGTTGTTCACGGTGCCCTCCATGTACAGCTTCATGTGCATGTTCTCCTTGATCAGCTCGCTCATAGGTCCAGGGTTCTCCTCCACGTCTCCAGCCTGCTTCAGCAGGCTGAAGTTAGTAGCTCCGCTTCCGGATCCCCCGGGGAGCATGTCAAGGTCAAAATCGTCAAGAGCGTCAGCAGGCAGCATATCAAGGTCAAAGTCGTCAAGGGCATCGGCTGGGAgCATGTCTAAgTCAAAATCGTCAAGGGCGTCGGCCGGCCCGCCGCTTTcgcacGCCCTGGCAATCGAGATGCTGGACAGGCATCATACCCACTTCTGCCCCCTGGAAGGCGAGTCATGGCAAGACTTTCTGCGGAACAACGCCAAGTCATTCCGCTGTGCTCTCCTCTCACATCGCGACGGGGCTAAAGTGCATCTCGGCACCCGCCCAACAGAGAAACAGTACGAAACCCTGGAAAATCAGCTCGCGTTCCTGTGTCAGCAAGGCTTCTCCCTGGAGAACGCACTGTACGCTCTGTCCGCCGTGGGCCACTTTACACTGGGCTGCGTATTGGAGGATCAGGAGCATCAAGTAGCAAAAGAGGAAAGAGAGACACCTACCACCGATTCTATGCCTGACTGTGGCGGGTGAGCTTAGGGGGCCTCCGCTCCAGCTCGACACCGGGCAGCTGCTGAAGATCGCGAAGAGAGGGGGAGTAACAGCGGTAGAGGCAGTGCACGCCTGGCGCAATGCGCTCACCGGGGCCCCCTTGAACCTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAGCAATGGGGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTTTAACAAAATATTAACGCTTACAATTTAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAA"

	// run blast
	matches, err := blast(id, seq, true, []string{testDB}, []string{}, 10, blastWriter()) // any match over 10 bp

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
			if targ.entry == m.entry && targ.queryStart == m.queryStart && targ.queryEnd == m.queryEnd {
				return
			}
		}

		t.Errorf("failed to find match %v in fragment matches", targ)
	}

	matchesContain(match{
		entry:      "gnl|addgene|107006(circular)",
		queryStart: 0,
		queryEnd:   72,
	})
}

// test that we can filter out overlapping regions from blast results
// and those that are up against the edge of the fragment
func Test_filter(t *testing.T) {
	// test fragment with 3 matches that should be removed
	matches := []match{
		// shouldn't be removed
		match{
			entry:      "m1",
			queryStart: 15,
			queryEnd:   19,
		},
		// should be removed because it fits within m3
		match{
			entry:      "m2",
			queryStart: 29,
			queryEnd:   34,
		},
		// shouldn't be removed
		match{
			entry:      "m3",
			queryStart: 29,
			queryEnd:   35,
		},
		// shouldn't be removed
		match{
			entry:      "m4",
			queryStart: 31,
			queryEnd:   72,
		},
	}

	newMatches := filter(matches, 72, 3) // keep all fragments larger than 3bp (all of them)

	if len(newMatches) != 5 {
		t.Errorf("%d filtered matches found on test fragment, 5 expected: %v", len(newMatches), newMatches)
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
	c.PCRMaxOfftargetTm = 40.0

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
	conf.PCRMaxOfftargetTm = 40.0

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
				uniqueID:    "addgene:107006(circular)0",
				seq:         "AGTATAGTAGGTAGTCATTCTT",
				queryStart:  0,
				queryEnd:    23,
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
				t.Errorf("parentMismatch() error = %+v, wantErr %+v", err, tt.wantErr)
				return
			}
			if gotMismatch != tt.wantMismatch {
				t.Errorf("parentMismatch() gotMismatch = %+v, want %+v", gotMismatch, tt.wantMismatch)
			}
			if !reflect.DeepEqual(gotMatch, tt.wantMatch) {
				t.Errorf("parentMismatch() gotMatch = %+v, want %+v", gotMatch, tt.wantMatch)
			}
		})
	}
}
