package defrag

import (
	"path"
	"testing"
)

// Test reading of a FASTA file
func TestRead(t *testing.T) {
	file := path.Join("..", "..", "test", "multi.fasta")
	fragments, err := read(file)

	if err != nil {
		t.Error(err)
		return
	}

	if len(fragments) != 5 {
		t.Errorf("failed to load fragments, len=%d, slice=%v", len(fragments), fragments)
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
