package assemble

// assembly is a slice of nodes ordered by the nodes
// distance from the end of the target vector
type assembly struct {
	// estimated cost of making this assembly
	cost float32

	// nodes, ordered by distance from the "end" of the vector
	nodes []node
}

// add a node to the start of this assembly.
// Update the cost of the assembly to include the link between the new first node and the one after it.
// store the node's id in the list of node ids.
// Complete if a node has matched up onto itself across the zero-index.
// TODO: incorporate cost estimate of the last node in an assembly
func (a *assembly) add(n node) (newAssembly assembly, created, complete bool) {
	// do not create a new assembly if the total number of nodes
	// would be less than the upper limit
	if a.len()+1 > conf.Fragments.MaxCount {
		return newAssembly, false, false
	}

	// we have not exceeded the fragment count limit, will be making new assembly
	created = true

	// add to list of nodes, update cost, and return
	if a.len() > 0 {
		// we've complete an assembly if the id of a node is the same
		complete = n.uniqueID == a.nodes[a.len()-1].uniqueID

		return assembly{
			nodes: append([]node{n}, a.nodes...),
			cost:  n.costTo(a.nodes[0]),
		}, created, complete
	}

	// create the start of this assembly, no other nodes
	return assembly{
		nodes: []node{n},
		cost:  0,
	}, created, false
}

// contains returns if the id of the node has already been seen in this assembly
func (a *assembly) contains(n node) (isContained bool) {
	for _, otherN := range a.nodes {
		// they're the same if they have the same id and start index
		// id isn't enough by itself because there may be multiple with the same
		// entry id in the BLAST db
		if otherN.id == n.id && otherN.start == n.start {
			return true
		}
	}
	return false
}

// len returns len(nodes)
func (a *assembly) len() int {
	return len(a.nodes)
}
