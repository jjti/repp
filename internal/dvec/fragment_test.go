package dvec

import "testing"

// Test reading of a FASTA file
func TestRead(t *testing.T) {
	fragments, err := Read("../../test/frag/fasta.fa")

	if err != nil {
		t.Errorf("failed in Read: %s", err.Error())
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
