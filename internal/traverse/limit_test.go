package traverse

import (
	"testing"
)

// test the ability to calculate distances from each node to the end
// of a mock node sequence range
//
// we want to know how far each node is the from the "end of the vector",
// so, looking at each node from the start, we can find which will be
// included in assemblies with less than a certain number of fragments
func Test_distance(t *testing.T) {
	conf.Synthesis.MaxLength = 2

	f1 := node{
		start: 8,
		end:   13,
		entry: true,
	}
	f2 := node{
		start: 12,
		end:   17,
	}
	f3 := node{
		start: 13,
		end:   18,
	}
	f4 := node{
		start:    15,
		end:      20,
		terminal: true,
	}
	f5 := node{
		start:    17,
		end:      21,
		terminal: true,
	}
	nodes := []node{f1, f2, f3, f4, f5}

	// using a target vector sequence length of 10
	dists := distance(nodes, 10)

	// should only be 5 keys (one for each node)
	if len(dists) != 5 {
		t.Errorf("failed, distance map has %d keys, should have 5", len(dists))
	}

	// make sure the distance in the distance map is only as long as
	// we expect it to be
	testIndex := 0
	checkDist := func(key node, dist int) {
		testIndex++
		if d, _ := dists[key]; d != dist {
			t.Errorf("failed, distance for test %d is %d, not %d", testIndex, d, dist)
		}
	}
	checkDist(f1, 3)
	checkDist(f2, 2)
	checkDist(f3, 2)
	checkDist(f4, 1)
	checkDist(f5, 1)
	return
}
