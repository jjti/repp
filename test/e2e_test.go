package test

import (
	"path"
	"testing"

	"github.com/jjtimmons/defrag/config"
	"github.com/jjtimmons/defrag/internal/defrag"
)

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
