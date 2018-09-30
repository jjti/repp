package traverse

import "testing"

// test ability to calculate synthetic fragment count distance
// between two fragmennts
func Test_synthDist(t *testing.T) {
	conf.Synthesis.MaxLength = 2
	n := node{
		start: 2,
		end:   4,
	}

	testIndex := 0
	assert := func(dist1, dist2 int) {
		testIndex++
		if dist1 != dist2 {
			t.Errorf("got %d, expected %d in test %d", dist1, dist2, testIndex)
		}
	}

	assert(n.synthDist(node{start: 3, end: 5}), 0)
	assert(n.synthDist(node{start: 5, end: 10}), 1)
	assert(n.synthDist(node{start: 1, end: 3}), 0)
	assert(n.synthDist(node{start: 0, end: 1}), 1)
	return
}
