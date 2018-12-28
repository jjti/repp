package defrag

import (
	"path"
	"testing"
)

// Test reading of a FASTA file
func Test_read(t *testing.T) {
	type fileRead struct {
		name      string
		file      string
		fragCount int
	}

	files := []fileRead{
		{
			"113726(circular)",
			path.Join("..", "..", "test", "113726(circular).parent"),
			1,
		},
		{
			"multi.fasta",
			path.Join("..", "..", "test", "multi.fasta"),
			5,
		},
	}

	for _, f := range files {
		fragments, err := read(f.file)

		if err != nil {
			t.Error(err)
		}

		if len(fragments) != f.fragCount {
			t.Errorf("failed to load fragments, len=%d, expected=%d", len(fragments), f.fragCount)
		}

		for _, f := range fragments {
			// ensure we got an ID
			if len(f.ID) < 1 {
				t.Error("failed to load an ID for a Fragment from FASTA")
			}

			// ensure we got a Seq
			if len(f.Seq) < 1 {
				t.Errorf("failed to parse a sequence for Fragment %s", f.ID)
			}
		}
	}
}
