package blast

import (
	"testing"

	"github.com/jjtimmons/decvec/internal/frag"
)

// test that we can filter out overlapping regions from blast results
// and those that are up against the enge of the fragment
func TestFilter(t *testing.T) {
	// test fragment with 3 matches that should be removed
	matches := []frag.Match{
		// should be removed for ending before 0-index
		frag.Match{
			ID:    "m1",
			Start: 15,
			End:   9,
		},
		// should be removed because it fits within m3
		frag.Match{
			ID:    "m2",
			Start: 29,
			End:   34,
		},
		// shouldn't be removed
		frag.Match{
			ID:    "m3",
			Start: 29,
			End:   35,
		},
		// shouldn't be removed
		frag.Match{
			ID:    "m4",
			Start: 31,
			End:   72,
		},
		// should be removed because it starts past 2x the
		// target's sequence length
		frag.Match{
			ID:    "m5",
			Start: 40,
			End:   79,
		},
	}

	newMatches := filter(matches, 24, 48)

	// make sure they're gone
	if len(newMatches) != 3 {
		t.Errorf("%d matches found on test fragment, 3 expected: %v", newMatches, matches)
	}

	// make sure m2 has been removed
	for _, m := range newMatches {
		if m.ID == "m2" {
			t.Error("m2 found in resulting matches, should have been removed")
		}
	}
}
