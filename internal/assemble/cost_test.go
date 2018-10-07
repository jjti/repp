package assemble

import "testing"

// test the ability to find the distance from each node to the end of the vector
func Test_cost(t *testing.T) {
	// start = 10, end = 20
	n1 := node{
		start:    10,
		end:      15,
		terminus: false,
	}
	n2 := node{
		start:    14,
		end:      20,
		terminus: true,
	}
	nodes := []node{n1, n2}

	costs := cost(nodes)

	if len(costs) < len(nodes) {
		t.Errorf("finding costs for %d nodes, expected %d", len(costs), len(nodes))
	}

	c1, _ := costs[n1]
	c2, _ := costs[n2]

	if c2 > c1 {
		t.Errorf("n2 costs %f which is more than n1, %f", c2, c1)
	}
}
