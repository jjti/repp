package assemble

// assembly is a slice of nodes ordered by the nodes
// distance from the end of the target vector
type assembly struct {
	// estimated cost of making this assembly
	cost float32

	// nodes, ordered by distance from the "end" of the vector
	nodes []node

	// ids of nodes in this assembly
	// used for determining whether this assembly holds a given
	// node later (for avoiding BLACKLISTED_NODES)
	ids map[string]bool
}

// add a node to the start of this assembly.
// update the cost of the assembly to include the link between
// the new first fragment and the one after it.
// store the node's id in the list of node ids.
// TODO: incorporate cost estimate of the last node in an assembly
func (a *assembly) add(n node) *assembly {
	if len(a.nodes) > 0 {
		a.cost += n.costTo(a.nodes[0])
	}
	a.ids[n.id] = true
	a.nodes = append([]node{n}, a.nodes...)
	return a
}

// has returns if the id of the node has already been seen in this assembly.
func (a *assembly) has(n node) bool {
	_, hasNode := a.ids[n.id]
	return hasNode
}
