package defrag

import (
	"log"
	"math"
	"sort"

	"github.com/jjtimmons/defrag/config"
)

// assemble the BLAST matches into assemblies that span the target sequence
//
// First build up assemblies, creating all possible assemblies that are
// beneath the upper-bound limit on the number of fragments
//
// Then find the pareto optimal solutions that minimize either the cost
// of the assembly or the number of fragments (relative to the other
// assembly plans)
//
// Then, for each pareto optimal solution, traverse the assembly and
// "fill-in" the nodes. Create primers on the node if it's a PCR Fragment
// or create a sequence to be synthesized if it's a synthetic fragment.
// Error out and repeat the build stage if a node fails to be filled
func assemble(matches []Match, seq string, conf *config.Config) [][]Fragment {
	// map fragment Matches to nodes
	var nodes []node
	for _, m := range matches {
		nodes = append(nodes, new(m, len(seq), conf))
	}

	// build up slice of assemblies that could, within the upper-limit on
	// fragment count, be assembled to make the target vector
	assemblies := build(nodes, conf.Fragments.MaxCount)

	// build up a map from fragment count to a sorted list of assemblies with that number
	sortedAssemblies := countMap(assemblies)

	// sort the keys of the assemblies (their fragment count) so we start with the smaller
	// assemblies first and work out way up
	assemblyCounts := []int{}
	for count := range sortedAssemblies {
		assemblyCounts = append(assemblyCounts, count)
	}
	sort.Ints(assemblyCounts)

	// "fill-in" the fragments (create synthetic fragments and primers)
	// skip assemblies that wouldn't be less than the current cheapest assembly
	minCostAssembly := math.MaxFloat64
	filled := make(map[int][]Fragment)
	for _, count := range assemblyCounts {
		// get the first assembly that fills properly (cheapest workable solution)
		for _, singleAssembly := range sortedAssemblies[count] {
			if singleAssembly.cost > minCostAssembly {
				// skip this and the rest with this count, there's another
				// cheaper option with fewer fragments
				break
			}

			filledFragments, err := singleAssembly.fill(seq, conf)
			if err != nil {
				log.Fatal(err)
			}

			// if a node in the assembly fails to be prepared,
			// remove all assemblies with the node and try again
			if filledFragments != nil {
				assemblyCost := 0.0
				for _, f := range filledFragments {
					assemblyCost += f.Cost
				}

				if assemblyCost > minCostAssembly {
					continue // wasn't actually cheaper, keep trying
				}

				// add this list of fragments to the list of such
				filled[len(filledFragments)] = filledFragments

				// store this as the new cheapest assembly
				minCostAssembly = assemblyCost

				// break, don't look at other possible assemblies because this one worked
				break
			}
		}
	}

	var found [][]Fragment
	for _, frags := range filled {
		found = append(found, frags)
	}
	return found
}

// countMap returns a map from the number of fragments in a build
// to a slice of builds with that number of fragments, sorted by their cost
//
// It also removes builds that cost more and have more fragments than
// cheaper assemblies with fewer fragments (hence name).
// So if it's $50 to make a vector with 2 fragments, this will prune out
// the $75 vector with 3 fragments. It's not worth making an assembly with 3
// fragments at $75, because it would be cheaper (and better) to do it with 2
func countMap(assemblies []assembly) map[int][]assembly {
	countToAssemblies := make(map[int][]assembly)
	for _, a := range assemblies {
		// The "-1"s here are because the assemblies are circular and
		// their last node is the same as the first. The total number
		// of nodes/fragments in the assembly is actually one less than
		// the assembly's length
		if as, ok := countToAssemblies[a.len()-1]; ok {
			countToAssemblies[a.len()-1] = append(as, a)
		} else {
			countToAssemblies[a.len()-1] = []assembly{a}
		}
	}

	for count := range countToAssemblies {
		sort.Slice(countToAssemblies[count], func(i, j int) bool {
			return countToAssemblies[count][i].cost < countToAssemblies[count][j].cost
		})
	}

	return countToAssemblies
}
