package assemble

import (
	"sort"
)

// pareto returns a map from the number of fragments in a build
// to a slice of builds with that number of fragments, sorted by their cost
func pareto(assemblies []assembly) map[int][]assembly {
	paretos := make(map[int][]assembly)
	for _, a := range assemblies {
		if as, ok := paretos[a.len()]; ok {
			paretos[a.len()] = append(as, a)
		} else {
			paretos[a.len()] = []assembly{a}
		}
	}

	for count := range paretos {
		sort.Slice(paretos[count], func(i, j int) bool {
			return paretos[count][i].cost < paretos[count][j].cost
		})
	}

	return paretos
}
