package blast

import (
	"testing"

	"github.com/jjtimmons/decvec/internal/frag"
)

// test that we can filter out overlapping regions from blast results
// and those that are up against the enge of the fragment
func TestFilter(t *testing.T) {
	// test fragment with 3 matches that should be removed
	f := frag.Fragment{
		ID:  "testFragment",
		Seq: "atgctagctagctagctagctagc",
		Matches: []frag.Match{
			// should be removed for starting at 0-index
			frag.Match{
				ID:    "m1",
				Start: 0,
				End:   5,
			},
			// should be removed because it fits within m3
			frag.Match{
				ID:    "m2",
				Start: 5,
				End:   10,
			},
			frag.Match{
				ID:    "m3",
				Start: 5,
				End:   11,
			},
			// should be removed because it ends near end
			frag.Match{
				ID:    "m4",
				Start: 7,
				End:   48,
			},
			// should be removed because it starts past 2x the
			// target's sequence length
			frag.Match{
				ID:    "m5",
				Start: 50,
				End:   55,
			},
		},
	}

	filter(&f)

	// make sure they're gone
	if len(f.Matches) != 3 {
		t.Errorf("%d matches found on test fragment, 3 expected", len(f.Matches))
	}

	// make sure m2 has been removed
	for _, m := range f.Matches {
		if m.ID == "m2" {
			t.Error("m2 found in resulting matches, should have been removed")
		}
	}
}
