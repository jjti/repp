package assemble

// build is for building up an ASSEMBLIES_LIST
//
// ASSEMBLIES_LIST is a map from the number of fragments in an assembly to
// a list of assemblies that have that number of fragments within them
//
// It is created by traversing a DAG in reverse order:
//
// number of additional nodes to consider for each node
// (ie additional number of fragments to try synthing to)
// synth_count = math.max(10, 0.05 * len(nodes))
//
// // traverse the nodes
// foreach this.node (sorted in reverse order):
// 	foreach that.node that this node overlaps with + synth_count:
//	  foreach assembly on that.node:
//      add this.node to the assembly to create a new assembly, store on this.node
// after the last node is traversed, add an empty node and synth to synth_count other nodes
//
// // store in a map
// create a map from the number of fragments in each assembly to a list with assemblies
// 		containing that many assemblies
func build(nodes []node) (assemblies map[int]assembly) {

	return assemblies
}
