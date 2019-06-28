package repp

import (
	"path"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/jjtimmons/repp/config"
)

func Test_sequence_e2e(test *testing.T) {
	c := config.New()

	type testFlags struct {
		in       string
		out      string
		backbone string
		enzymes  []string
		filters  string
		dbs      []string
		addgene  bool
		igem     bool
	}

	tests := []testFlags{
		testFlags{
			path.Join("..", "..", "test", "input", "backbone.fa"),
			path.Join("..", "..", "test", "output", "backbone.json"),
			"pSB1A3",
			[]string{"PstI"},
			"2018,2019",
			[]string{},
			false,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2224001.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2224001.json"),
			"pSB1A3",
			[]string{"PstI"},
			"2017,2018,2019",
			[]string{},
			false,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "110056.fa"),
			path.Join("..", "..", "test", "output", "110056.json"),
			"",
			[]string{},
			"2019,2018",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2602025.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2602025.json"),
			"pSB1A3",
			[]string{"PstI"},
			"",
			[]string{},
			false,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2779020.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2779020.json"),
			"pSB1A3",
			[]string{"PstI"},
			"BBa_K277",
			[]string{},
			false,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_E0610.fa"),
			path.Join("..", "..", "test", "output", "BBa_E0610.json"),
			"pSB1C3",
			[]string{"EcoRI"},
			"2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,BBa_E061",
			[]string{},
			false,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_I5310.fa"),
			path.Join("..", "..", "test", "output", "BBa_I5310.json"),
			"pSB1C3",
			[]string{"EcoRI"},
			"",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2651001.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2651001.json"),
			"pSB1C3",
			[]string{"EcoRI"},
			"BBa_K265",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2779020.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2779020.json"),
			"pSB1A3",
			[]string{"PstI"},
			"BBa_K277", // no year filters needed
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K1085023.fa"),
			path.Join("..", "..", "test", "output", "BBa_K1085023.json"),
			"pSB1C3",
			[]string{"EcoRI"},
			"2009,2010,2011,2012,BBa_K108",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "113490.fa"),
			path.Join("..", "..", "test", "output", "113490.json"),
			"",
			[]string{},
			"2018,2019",
			[]string{},
			true,
			false,
		},
	}

	for _, t := range tests {
		sols := Sequence(NewFlags(t.in, t.out, t.backbone, t.filters, t.enzymes, t.dbs, t.addgene, t.igem, false))

		if len(sols) < 1 {
			test.Errorf("no solutions for %s", t.in)
		}

		for _, s := range sols {
			e := validateJunctions(s, c)
			if e != nil {
				test.Logf("failed making %s\n", t.in)
				test.Error(e)
			}
		}
	}
}

func Test_features(t *testing.T) {
	test1, conf := NewFlags(
		"p10 promoter, mEGFP, T7 terminator",
		filepath.Join("..", "..", "test", "output", "features.json"),
		"pSB1A3",
		"EcoRI",
		[]string{},
		[]string{},
		true,
		true,
		false,
	)

	test2, _ := NewFlags(
		"BBa_R0062,BBa_B0034,BBa_C0040,BBa_B0010,BBa_B0012",
		filepath.Join("..", "..", "test", "output", "igem.features.json"),
		"pSB1C3",
		"PstI",
		[]string{},
		[]string{},
		false,
		true,
		false,
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
		{
			"test end to end features creation using iGEM parts",
			args{
				flags: test2,
				conf:  conf,
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			sols := Features(tt.args.flags, tt.args.conf)

			if len(sols) < 1 {
				t.Failed()
			}

			for _, s := range sols {
				e := validateJunctions(s, conf)
				if e != nil {
					t.Error(e)
				}
			}
		})
	}
}

func Test_fragments(t *testing.T) {
	c := config.New()
	c.PCRMinLength = 10
	c.FragmentsMinHomology = 8
	c.FragmentsMaxHomology = 20

	type args struct {
		inputFragments []*Frag
		conf           *config.Config
	}
	tests := []struct {
		name             string
		args             args
		wantTargetVector *Frag
		wantFragments    []*Frag
	}{
		{
			"fragments with linear overlap",
			args{
				[]*Frag{
					&Frag{
						Seq:  "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
						conf: c,
					},
					&Frag{
						Seq:  "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
						conf: c,
					},
					&Frag{
						Seq:  "TCGACTAGCTAGAACTGATCTAGACGTGCTAGCTACA",
						conf: c,
					},
				},
				c,
			},
			&Frag{
				Seq: "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCGACTGATCACTAGCATCGACTAGCTAGAACTGATCTAG",
			},
			[]*Frag{
				&Frag{
					Seq:      "ACGTGCTAGCTACATCGATCGTAGCTAGCTAGCATCG",
					fragType: linear,
				},
				&Frag{
					Seq:      "AGCTAGCATCGACTGATCACTAGCATCGACTAGCTAG",
					fragType: linear,
				},
				&Frag{
					Seq:      "TCGACTAGCTAGAACTGATCTAGACGTGCTAGCTACA",
					fragType: linear,
				},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotTargetVector, gotFragments := fragments(tt.args.inputFragments, tt.args.conf)

			if !reflect.DeepEqual(gotTargetVector.Seq, tt.wantTargetVector.Seq) {
				t.Errorf("fragments() gotTargetVector = %v, want %v", gotTargetVector, tt.wantTargetVector)
			}

			if len(gotFragments) != len(tt.wantFragments) {
				t.Errorf("fragments() got %d fragments, expected %d", len(gotFragments), len(tt.wantFragments))
				return
			}

			for i, wantF := range tt.wantFragments {
				if wantF.Seq != gotFragments[i].Seq {
					t.Errorf("fragments() gotFragment.Seq = %v, want %v", gotFragments[i].Seq, wantF.Seq)
				}

				if wantF.fragType != gotFragments[i].fragType {
					t.Errorf("fragments() gotFragment.Type = %v, want %v", gotFragments[i].fragType, wantF.fragType)
				}
			}
		})
	}
}

// if an input fragment being built is exactly the same as one in a DB, it should be used
// as is and without PCR or any preparation
func Test_vector_single_vector(t *testing.T) {
	fs, c := NewFlags(
		path.Join("..", "..", "test", "input", "109049.addgene.fa"),
		path.Join("..", "..", "test", "output", "109049.output.json"),
		"",
		"",
		[]string{},
		[]string{},
		true,
		false,
		false,
	)

	assemblies := Sequence(fs, c) // use addgene database

	if !strings.Contains(assemblies[0][0].URL, "109049") {
		t.Fatal("failed to use 109049 to build the vector")
	}
}
