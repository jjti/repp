package assemble

import (
	"fmt"
	"math"
)

// traverse is for finding the minimum cost estimated assembly and traversing
// the assembly while "filling in the nodes." This means creating primers on
// PCR fragments and creating synthetic fragments
//
// if the built-up assembly is too large, or if any of the nodes are deemed
// invalid, the engine should start over again with a new starting node
func traverse(nodes []node) error {
	// get costs of the nodes
	costs := cost(nodes)

	// get minimum cost starting point
	start, err := minCost(costs)
	if err != nil {
		return err
	}

	// begin to traverse the potential assembly, creating
	// synthetic fragments and primers to build it
	return nil
}

// minCost returns the cheapest available entry node
// returns an error if there's no valid entry node
func minCost(costs map[node]float32) (node, error) {
	// find the assembly with the minimum estimated cost
	// first set the starting node as one that begins on the first entry bp
	// and is synthesized to the next fragment after it
	minCost := float32(math.MaxFloat32)
	var minCostNode node
	emptyNode := node{}

	for n, cost := range costs {
		if n.entry && cost < minCost {
			minCost = cost
			minCostNode = n
		}
	}

	if minCostNode == emptyNode {
		// we were unable to find a valid entry node, error out
		return minCostNode, fmt.Errorf("no valid entry node found in cost map with %d nodes", len(costs))
	}
	return minCostNode, nil
}
