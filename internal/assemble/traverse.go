package assemble

// traverse is for finding the minimum cost estimated assembly and traversing
// the assembly while "filling in the nodes." This means creating primers on
// PCR fragments and creating synthetic fragments
//
// if the built-up assembly is too large, or if any of the nodes are deemed
// invalid, the engine should start over again with a new starting node
func traverse(nodes []node) {
	// get costs of the nodes
	costs := cost(nodes)

	// find the assembly with the minimum estimated cost
	// first set the starting node as one that begins on the first entry bp
	// and is synthesized to the next fragment after it
	var minCost node
	for n, cost := range costs {

	}
}
