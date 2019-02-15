package test

import (
	"path"
	"testing"

	"github.com/jjtimmons/defrag/config"
	"github.com/jjtimmons/defrag/internal/defrag"
)

func Test_Sequence(t *testing.T) {
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
		// testFlags{
		// 	path.Join("input", "BBa_E0610.fa"),
		// 	path.Join("output", "BBa_E0610.json"),
		// 	"pSB1C3",
		// 	"EcoRI",
		// 	"2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,BBa_E061",
		// 	[]string{},
		// 	true,
		// 	true,
		// },
		testFlags{
			path.Join("input", "BBa_K2602025.fa"),
			path.Join("output", "BBa_K2602025.json"),
			"pSB1A3",
			"PstI",
			"",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("input", "BBa_I5310.fa"),
			path.Join("output", "BBa_I5310.json"),
			"pSB1C3",
			"EcoRI",
			"",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("input", "BBa_K077557.fa"),
			path.Join("output", "BBa_K077557.json"),
			"pSB1C3",
			"EcoRI",
			"2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,BBa_K0077",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("input", "BBa_K2651001.fa"),
			path.Join("output", "BBa_K2651001.json"),
			"pSB1C3",
			"EcoRI",
			"BBa_K265",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("input", "BBa_K2779020.fa"),
			path.Join("output", "BBa_K2779020.json"),
			"pSB1A3",
			"PstI",
			"BBa_K277", // no year filters needed
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("input", "BBa_K1085023.fa"),
			path.Join("output", "BBa_K1085023.json"),
			"pSB1C3",
			"EcoRI",
			"2009,2010,2011,2012,BBa_K108",
			[]string{},
			true,
			true,
		},
	}

	c := config.New()

	for _, t := range tests {
		sols := defrag.Sequence(defrag.NewFlags(t.in, t.out, t.backbone, t.enzyme, t.filters, t.dbs, t.addgene, t.igem))

		for _, s := range sols {
			defrag.ValidateJunctions(s, c)
		}
	}

	t.Fatal("fail (dev)")
}

func Test_Features(t *testing.T) {
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
			"lacI promoter, EGFP, rrnBT1 terminator:rev",
			path.Join("output", "lacI-EGFP-rrnBT1.features.json"),
			"pSB1A3",
			"PstI",
			"",
			[]string{},
			true,
			true,
		},
	}

	c := config.New()

	for _, t := range tests {
		sols := defrag.Features(defrag.NewFlags(t.in, t.out, t.backbone, t.enzyme, t.filters, t.dbs, t.addgene, t.igem))

		for _, s := range sols {
			defrag.ValidateJunctions(s, c)
		}
	}

	t.Fatal("fail (dev)")
}
