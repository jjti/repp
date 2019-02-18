package defrag

import (
	"path"
	"strings"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

// if an input fragment being built is exactly the same as one in a DB, it should be used
// as is and without PCR or any preparation
func Test_vector_single_vector(t *testing.T) {
	fs, c := NewFlags(
		path.Join("..", "..", "test", "input", "109049.addgene.fa"),
		path.Join("..", "..", "test", "output", "109049.output.json"),
		"",
		"",
		"",
		[]string{},
		true,
		false,
		false,
	)

	_, _, assemblies, err := sequence(fs, c) // use addgene database
	if err != nil {
		t.Error(err)
	}

	if !strings.Contains(assemblies[0][0].ID, "109049") && !strings.Contains(assemblies[1][0].ID, "109049") {
		t.Fatal("failed to use 109049 to build the vector")
	}
}

func Test_sequence(test *testing.T) {
	c := config.New()

	type testFlags struct {
		in       string
		out      string
		backbone string
		enzyme   string
		filters  string
		dbs      []string
		addgene  bool
		igem     bool
	}

	tests := []testFlags{
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2602025.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2602025.json"),
			"pSB1A3",
			"PstI",
			"",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2779020.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2779020.json"),
			"pSB1A3",
			"PstI",
			"BBa_K277",
			[]string{},
			false,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_E0610.fa"),
			path.Join("..", "..", "test", "output", "BBa_E0610.json"),
			"pSB1C3",
			"EcoRI",
			"2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,BBa_E061",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_I5310.fa"),
			path.Join("..", "..", "test", "output", "BBa_I5310.json"),
			"pSB1C3",
			"EcoRI",
			"",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K077557.fa"),
			path.Join("..", "..", "test", "output", "BBa_K077557.json"),
			"pSB1C3",
			"EcoRI",
			"2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,BBa_K0077",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2651001.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2651001.json"),
			"pSB1C3",
			"EcoRI",
			"BBa_K265",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2779020.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2779020.json"),
			"pSB1A3",
			"PstI",
			"BBa_K277", // no year filters needed
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K1085023.fa"),
			path.Join("..", "..", "test", "output", "BBa_K1085023.json"),
			"pSB1C3",
			"EcoRI",
			"2009,2010,2011,2012,BBa_K108",
			[]string{},
			true,
			true,
		},
	}

	for _, t := range tests {
		sols := Sequence(NewFlags(t.in, t.out, t.backbone, t.enzyme, t.filters, t.dbs, t.addgene, t.igem, false))

		if len(sols) < 1 {
			test.Fail()
		}

		for _, s := range sols {
			ValidateJunctions(t.in, s, c, test)
		}
	}
}
