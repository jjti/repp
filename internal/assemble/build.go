package assemble

// build is for building up an ASSEMBLIES_LIST
//
// ASSEMBLIES_LIST is a map from the number of fragments in an assembly to
// a list of assemblies that have that number of fragments within them
//
// It is created by traversing a DAG in reverse, looking for nodes that overlap
// with the current node. An "edge" (limited meaning here tbh) is either PCR
// or synthesis. Ie the nodes in the DAG are joined by either PCR and
// homology between the nodes through PCR or synthesis (new fragments are
// created that create the homology between the nodes)
func build() {

}
