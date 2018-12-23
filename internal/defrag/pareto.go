package defrag

import (
	"sort"
)

// pareto returns a map from the number of fragments in a build
// to a slice of builds with that number of fragments, sorted by their cost
//
// It also removes builds that cost more and have more fragments than
// cheaper assemblies with fewer fragments (hence name).
// So if it's $50 to make a vector with 2 fragments, this will prune out
// the $75 vector with 3 fragments. It's not worth making an assembly with 3
// fragments at $75, because it would be cheaper (and better) to do it with 2
func pareto(assemblies []assembly) map[int][]assembly {
	paretos := make(map[int][]assembly)
	for _, a := range assemblies {
		// The "-1"s here are because the assemblies are circular and
		// their last node is the same as the first. The total number
		// of nodes/fragments in the assembly is actually one less than
		// the assembly's length
		if as, ok := paretos[a.len()-1]; ok {
			paretos[a.len()-1] = append(as, a)
		} else {
			paretos[a.len()-1] = []assembly{a}
		}
	}

	for count := range paretos {
		sort.Slice(paretos[count], func(i, j int) bool {
			return paretos[count][i].cost < paretos[count][j].cost
		})
	}

	return paretos
}
