package assemble

import (
	"testing"
)

// test the ability to find the minimum cost node
func Test_minCost(t *testing.T) {
	// create the test cost map
	n1 := node{
		id:    "n1",
		entry: false,
	}
	n2 := node{
		id:    "n2",
		entry: true,
	}
	n3 := node{
		id:    "n3",
		entry: true,
	}
	n4 := node{
		id:    "n4",
		entry: false,
	}

	costsWithEntries := map[node]float32{
		n1: 0.9,
		n2: 1.5,
		n3: 2.0,
		n4: 0.2,
	}
	n, err := minCost(costsWithEntries)
	if err != nil {
		t.Error(err)
	}

	if n != n2 {
		t.Errorf("expected starting node to be %v, got %v", n2, n)
	}

	// in a map of costs where there's no compatible entry nodes,
	// we should have a synthetic vector returned instead that will
	// bridge the region from the start of the vector to the first
	// actual building fragment
	costsWithoutEntries := map[node]float32{
		n1: 0.9,
		n4: 0.2,
	}
	n, err = minCost(costsWithoutEntries)
	// should be a null/undeifned node object
	if err == nil {
		t.Errorf("expected minCost to throw an error without a valid starting node, none thrown")
	}
}
