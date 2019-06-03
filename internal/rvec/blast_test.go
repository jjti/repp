package rvec

import (
	"path"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/jjtimmons/rvec/config"
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
		entry:      "gnl|addgene|107006",
		queryStart: 0,
		queryEnd:   72,
	})
}

// test that we can filter out overlapping regions from blast results
// and those that are up against the edge of the fragment
func Test_cull(t *testing.T) {
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

	newMatches := cull(matches, 72, 3)

	// make sure m2 has been removed
	for _, m := range newMatches {
		if m.entry == "m2" {
			t.Error("m2 found in resulting matches, should have been removed")
		}
	}

	if len(newMatches) != 3 {
		t.Errorf("%d filtered matches found on test fragment, 3 expected: %v", len(newMatches), newMatches)
	}
}

func Test_isMismatch(t *testing.T) {
	c := config.New()
	c.PCRMaxOfftargetTm = 40.0

	type args struct {
		sequence string
		match    match
	}
	tests := []struct {
		name string
		args args
		want bool
	}{
		{
			"find a mismatching primer",
			args{
				sequence: "gtccgcgtcgtcgtcat",
				match: match{
					seq:     "atgacgacgacgcggac",
					forward: false,
				},
			},
			true,
		},
		{
			"no false positive mistmatch",
			args{
				sequence: "gtccgcgtcgtcgtcat",
				match: match{
					seq:     "atgacgacgacgac",
					forward: false,
				},
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := isMismatch(tt.args.sequence, tt.args.match, c); got != tt.want {
				t.Errorf("isMismatch() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_parentMismatch(t *testing.T) {
	testDB, _ := filepath.Abs(path.Join("..", "..", "test", "db", "db"))

	conf := config.New()
	conf.PCRMaxOfftargetTm = 35.0

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
				"gnl|addgene|107006",
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
				"gnl|addgene|107006",
			},
			true,
			match{
				entry:       "addgene:107006",
				uniqueID:    "addgene:1070060",
				seq:         "AGTATAGTAGGTAGTCATTCTT",
				querySeq:    "AGTATAGGATAGGTAGTCATTCTT",
				queryStart:  0,
				queryEnd:    23,
				mismatching: 2,
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			mismatchResult := parentMismatch([]Primer{Primer{Seq: tt.args.primer}}, tt.args.parent, testDB, conf)
			gotMismatch := mismatchResult.wasMismatch
			gotMatch := mismatchResult.m
			err := mismatchResult.err

			if (err != nil) != tt.wantErr {
				t.Errorf("parentMismatch() error = %+v, wantErr %+v", err, tt.wantErr)
				return
			}
			if gotMismatch != tt.wantMismatch {
				t.Errorf("parentMismatch() gotMismatch = %+v, want %+v", gotMismatch, tt.wantMismatch)
			}

			// have to mutate the fields not included in expected set
			gotMatch.circular = false
			gotMatch.title = ""
			gotMatch.subjectStart = 0
			gotMatch.subjectEnd = 0
			gotMatch.forward = false

			if !reflect.DeepEqual(gotMatch, tt.wantMatch) {
				t.Errorf("parentMismatch() gotMatch = %+v, want %+v", gotMatch, tt.wantMatch)
			}
		})
	}
}

func Test_queryDatabases(t *testing.T) {
	type args struct {
		entry string
		dbs   []string
	}
	tests := []struct {
		name    string
		args    args
		wantF   Frag
		wantErr bool
	}{
		{
			"query pSB1A3",
			args{
				entry: "pSB1A3",
				dbs:   []string{config.IGEMDB, config.AddgeneDB},
			},
			Frag{
				ID: "pSB1A3",
				db: config.IGEMDB,
			},
			false,
		},
		{
			"fail at nonsense frag",
			args{
				entry: "jahf9a8f9",
				dbs:   []string{config.IGEMDB, config.AddgeneDB},
			},
			Frag{},
			true,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotF, err := queryDatabases(tt.args.entry, tt.args.dbs)
			if (err != nil) != tt.wantErr {
				t.Errorf("queryDatabases() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if gotF.ID != tt.wantF.ID || gotF.db != tt.wantF.db {
				t.Errorf("queryDatabases() = %v, want %v", gotF, tt.wantF)
			}
		})
	}
}

func Test_blastdbcmd(t *testing.T) {
	testDB, _ := filepath.Abs(path.Join("..", "..", "test", "db", "db"))

	type args struct {
		entry string
		db    string
	}
	tests := []struct {
		name    string
		args    args
		wantErr bool
	}{
		{
			"find 107006",
			args{
				entry: "gnl|addgene|107006",
				db:    testDB,
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			_, _, err := blastdbcmd(tt.args.entry, tt.args.db)
			if (err != nil) != tt.wantErr {
				t.Errorf("blastdbcmd() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
		})
	}
}
