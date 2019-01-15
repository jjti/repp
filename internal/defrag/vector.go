package defrag

import (
	"fmt"
	"log"
	"math"
	"os"
	"sort"
	"strings"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
)

// VectorCmd is the root of the `defrag vector` functionality
//
// the goal is to find an "optimal" assembly vector with:
// 	1. the fewest fragments
// 	2. the lowest overall assembly cost ($)
// and, secondarily:
//	3. no duplicate end regions between Gibson fragments
// 	4. no inverted repeats in the junctions
// 	5. no off-target binding sites in the parent vectors
//	6. low primer3 penalty scores
func VectorCmd(cmd *cobra.Command, args []string) {
	conf := config.New()

	input, err := parseCmdFlags(cmd, conf)
	if err != nil {
		log.Fatalln(err)
	}

	vector(input, conf)

	os.Exit(0)
}

// vector builds a vector using reverse engineering
//
// TODO: redesign this for a linear input stretch
func vector(input *flags, conf *config.Config) [][]Frag {
	// read the target sequence (the first in the slice is used)
	fragments, err := read(input.in)
	if err != nil {
		log.Fatalf("failed to read in fasta files at %s: %v", input.in, err)
	}

	// set target fragment
	if len(fragments) > 1 {
		log.Printf(
			"warning: %d fragments were in %s. Only targeting the sequence of the first: %s",
			len(fragments),
			input.in,
			fragments[0].ID,
		)
	}
	targetFrag := fragments[0]

	// if a backbone was specified, add it to the sequence of the target frag
	if input.backbone.ID != "" {
		targetFrag.Seq += input.backbone.Seq
	}

	// get all the matches against the fragment
	matches, err := blast(&targetFrag, input.dbs, conf.Fragments.MinHomology)
	if err != nil {
		fmt.Printf("%+v", input.dbs)
		dbMessage := strings.Join(input.dbs, ", ")
		log.Fatalf("failed to blast %s against the dbs %s: %v", targetFrag.ID, dbMessage, err)
	}

	// build up the assemblies
	builds := buildVector(matches, targetFrag.Seq, conf)

	// try to write the JSON to the output file path
	if err := write(input.out, targetFrag, builds); err != nil {
		log.Fatal(err)
	}

	// return the builds (for e2e testing)
	return builds
}

// buildVector converts BLAST matches into assemblies spanning the target sequence
//
// First build up assemblies, creating all possible assemblies that are
// beneath the upper-bound limit on the number of fragments fully covering
// the target sequence
//
// Then find the pareto optimal solutions that minimize either the cost
// of the assembly or the number of fragments (relative to the other
// assembly plans)
//
// Then, for each pareto optimal solution, traverse the assembly and
// "fill-in" the nodes. Create primers on the Frag if it's a PCR Frag
// or create a sequence to be synthesized if it's a synthetic fragment.
// Error out and repeat the build stage if a Frag fails to be filled
func buildVector(matches []match, seq string, conf *config.Config) [][]Frag {
	// map fragment Matches to nodes
	var nodes []*Frag
	for _, m := range matches {
		nodes = append(nodes, newFrag(m, len(seq), conf))
	}

	// build up a slice of assemblies that could, within the upper-limit on
	// fragment count, be assembled to make the target vector
	assemblies := build(nodes, conf.Fragments.MaxCount, seq)

	// build up a map from fragment count to a sorted list of assemblies with that number
	sortedAssemblies := countMap(assemblies)

	// sort the keys of the assemblies (their fragment count). Start with the smaller
	// assemblies first and work out way up
	assemblyCounts := []int{}
	for count := range sortedAssemblies {
		assemblyCounts = append(assemblyCounts, count)
	}
	sort.Ints(assemblyCounts)

	// "fill-in" the fragments (create synthetic fragments and primers)
	// skip assemblies that wouldn't be less than the current cheapest assembly
	minCostAssembly := math.MaxFloat64
	filled := make(map[int][]Frag)
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

			// if a Frag in the assembly fails to be prepared,
			// remove all assemblies with the Frag and try again
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

	var found [][]Frag
	for _, frags := range filled {
		found = append(found, frags)
	}
	return found
}

// build builds up circular assemblies with less fragments than the build limit
//
// maxNodes is the maximum number of frags in a single assembly
// seq is the target sequence we're trying to build up
//
// It is created by traversing a DAG in forward order:
// 	foreach this.Frag (sorted in forward order):
// 	  foreach that.Frag that this Frag overlaps with + synthCount:
//	 	foreach assembly on that.Frag:
//    	    add this.Frag to the assembly to create a new assembly, store on this.Frag
func build(frags []*Frag, maxNodes int, seq string) (assemblies []assembly) {
	// number of additional frags try synthesizing to, in addition to those that
	// already have enough homology for overlap without any modifications for each Frag
	synthCount := int(math.Max(5, 0.05*float64(len(frags)))) // 5 of 5%, whichever is greater

	// sort with increasing start index
	sort.Slice(frags, func(i, j int) bool {
		return frags[i].start < frags[j].start
	})

	// create a starting assembly on each Frag including just itself
	for i, f := range frags {
		// edge case where the Frag spans the entire target vector... 100% match
		// it is the target vector. just return that as the assembly
		if len(f.Seq) >= len(seq) {
			assemblies = append(assemblies, assembly{
				frags:  []*Frag{f}, // just self
				cost:   f.Cost,     // just cost of procurement
				synths: 0,
			})
			break // can't top that
		}

		frags[i].assemblies = []assembly{
			assembly{
				frags:  []*Frag{f},  // just self
				cost:   f.costTo(f), // just PCR,
				synths: 0,           // no synthetic frags at start
			},
		}
	}

	// for every Frag in the list of increasing start index frags
	for i, f := range frags {
		// for every overlapping fragment + synthCount more
		for _, j := range f.reach(frags, i, synthCount) {
			// for every assembly on the reaching fragment
			for _, a := range frags[i].assemblies {
				// see if we can create a new assembly with this Frag included
				if newAssembly, created, complete := a.add(frags[j], maxNodes); created {
					if complete {
						// we've completed a circlular plasmid assembly
						// it has wrapped back onto itself
						assemblies = append(assemblies, newAssembly)
						// TODO: check if we can break here
					} else {
						// add to this Frag's list of assemblies
						frags[j].assemblies = append(frags[j].assemblies, newAssembly)
					}
				}
			}
		}
	}

	return
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
		// their last Frag is the same as the first. The total number
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
