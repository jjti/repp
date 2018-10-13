package blast

import (
	"testing"

	"github.com/jjtimmons/decvec/internal/dvec"
)

// test that we can filter out overlapping regions from blast results
// and those that are up against the edge of the fragment
func TestFilter(t *testing.T) {
	// test fragment with 3 matches that should be removed
	matches := []dvec.Match{
		// shouldn't be removed
		dvec.Match{
			Entry: "m1",
			Start: 15,
			End:   9,
		},
		// should be removed because it fits within m3
		dvec.Match{
			Entry: "m2",
			Start: 29,
			End:   34,
		},
		// shouldn't be removed
		dvec.Match{
			Entry: "m3",
			Start: 29,
			End:   35,
		},
		// shouldn't be removed
		dvec.Match{
			Entry: "m4",
			Start: 31,
			End:   72,
		},
	}

	newMatches := filter(matches)

	// make sure they're gone
	if len(newMatches) != 3 {
		t.Errorf("%d matches found on test fragment, 3 expected: %v", len(newMatches), newMatches)
	}

	// make sure m2 has been removed
	for _, m := range newMatches {
		if m.Entry == "m2" {
			t.Error("m2 found in resulting matches, should have been removed")
		}
	}
}
