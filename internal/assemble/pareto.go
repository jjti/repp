package assemble

import (
	"math"
	"sort"
)

// pareto returns a list of assemblies, each of which is a member of
// the pareto set for assemblies when minimizing for cost and
// fragment assembly count
//
// for example: return one assembly with 2 fragments if it's cheaper
// than all other 2 fragment assemblies. And, if it's cheaper than all
// other 3 fragment assemblies, there's no 3-fragment pareto optimal
// solution, the 2-fragment assembly sufficies
func pareto(assemblies []assembly) (paretoSet []assembly) {
	// build up the fragment count to assemblies list map
	countMap := make(map[int][]assembly)
	for _, a := range assemblies {
		if as, ok := countMap[a.len()]; ok {
			countMap[a.len()] = append(as, a)
		} else {
			countMap[a.len()] = []assembly{a}
		}
	}

	// minimum cost seen so far
	minCostOverall := float32(math.MaxFloat32)

	// frag counts sorted in increasing order
	var countKeys []int
	for count := range countMap {
		countKeys = append(countKeys, count)
	}
	sort.Ints(countKeys)

	// now check each fragment count to see if it's cheapest overall
	// and find the cheapest assembly within it
	for _, count := range countKeys {
		minAss := countMap[count][0]

		// check each assembly to see if it's cheaper that local min
		for _, assembly := range countMap[count] {
			if assembly.cost < minAss.cost {
				minAss = assembly
			}
		}

		// check if the cheapest assembly here is the cheapest overall.
		// keep if it is, ignore if it's not (there's another assembly
		// with fewer overall fragments and a lower overall cost).
		if minAss.cost < minCostOverall {
			paretoSet = append(paretoSet, minAss)
			minCostOverall = minAss.cost
		}
	}

	return
}
