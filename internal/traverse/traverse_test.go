package traverse

import (
	"fmt"
	"testing"

	"github.com/jjtimmons/decvec/internal/frag"
)

// test the ability to calculate distances from each fragment to the end
// of a mock fragment sequence range
//
// we want to know how far each fragment is the from the "end of the vector",
// so, looking at each fragment from the start, we can find which will be
// included in assemblies with less than a certain number of fragments
func Test_CalcFragDistance(t *testing.T) {
	f1 := frag.Match{
		ID:    "1",
		Start: 8,
		End:   13,
	}
	f2 := frag.Match{
		ID:    "2",
		Start: 12,
		End:   15,
	}
	f3 := frag.Match{
		ID:    "3",
		Start: 13,
		End:   18,
	}
	f4 := frag.Match{
		ID:    "4",
		Start: 15,
		End:   20,
	}
	f5 := frag.Match{
		ID:    "5",
		Start: 17,
		End:   21,
	}
	f := frag.Fragment{
		Seq:     "atgatctagc", // 10 bp, search range will be seq[10:21]
		Matches: []frag.Match{f1, f2, f3, f4, f5},
	}

	dists := calcFragDistance(&f, 3)

	// should only be 5 keys
	if len(dists) != 5 {
		t.Errorf("failed, distance map has %d keys, should have 5", len(dists))
	}

	fmt.Printf("%v", dists)

	// make sure the distance in the distance map is only as long as
	// we expect it to be
	checkDist := func(key frag.Match, dist int) {
		if d, _ := dists[key]; d != dist {
			t.Errorf("failed, distance map's value for %s is %d, not %d", key.ID, d, dist)
		}
	}
	checkDist(f1, 3)
	checkDist(f5, 1)
}
