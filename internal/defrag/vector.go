package defrag

import (
	"fmt"
	"log"
	"math"
	"os"
	"sort"
	"strings"
	"time"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
)

// VectorCmd is the interface for building a vector, from a sequence, from
// the command line
func VectorCmd(cmd *cobra.Command, args []string) {
	conf := config.New()
	start := time.Now()

	input, err := parseCmdFlags(cmd, conf)
	if err != nil {
		log.Fatalln(err)
	}

	target, builds, err := vector(input, conf)
	if err != nil {
		log.Fatalln(err)
	}

	// write the results to the filesystem at the out location
	fragments, err := read(input.in)
	if err != nil {
		log.Fatalln(err)
	}
	insertLength := len(fragments[0].Seq)

	if _, err := write(input.out, target, builds, insertLength, conf); err != nil {
		log.Fatalln(err)
	}

	elapsed := time.Since(start)
	fmt.Printf("%s\n\n", elapsed)

	os.Exit(0)
}

// VectorJSON is an interface for designing a vector from a JSON
// representation of the flags to pass to `vector`
func VectorJSON(json []byte) (output []byte, err error) {
	conf := config.New()

	input, err := parseJSONFlags(json, conf)
	defer os.Remove(input.in) // they're temporary files
	defer os.Remove(input.out)
	if err != nil {
		return
	}

	target, builds, err := vector(input, conf)

	// write the results to a JSON formatted []byte slice, return
	return write(input.out, target, builds, 0, conf)
}

// VectorFlags is for running an end to end vector design using pre-built flags
func VectorFlags(flags *Flags, conf *config.Config) {
	start := time.Now()

	target, builds, err := vector(flags, conf)
	if err != nil {
		log.Fatalln(err)
	}

	// write the results to the filesystem at the out location
	fragments, err := read(flags.in)
	if err != nil {
		log.Fatalln(err)
	}
	insertLength := len(fragments[0].Seq)

	if _, err := write(flags.out, target, builds, insertLength, conf); err != nil {
		log.Fatalln(err)
	}

	elapsed := time.Since(start)

	fmt.Printf("%s\n\n", elapsed)
}

// vector builds a vector using reverse engineering
//
// the goal is to find an "optimal" assembly vector with:
// 	1. the fewest fragments
// 	2. the lowest overall assembly cost ($)
// and, secondarily:
//	3. no duplicate end regions between Gibson fragments
// 	4. no inverted repeats in the junctions
// 	5. no off-target binding sites in the parent vectors
//	6. low primer3 penalty scores
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
func vector(input *Flags, conf *config.Config) (Frag, [][]*Frag, error) {
	// read the target sequence (the first in the slice is used)
	fragments, err := read(input.in)
	if err != nil {
		return Frag{}, nil, fmt.Errorf("failed to read in fasta files at %s: %v", input.in, err)
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
	fmt.Printf("Building %s\n", targetFrag.ID)

	// if a backbone was specified, add it to the sequence of the target frag
	if input.backbone.ID != "" {
		targetFrag.Seq += input.backbone.Seq
	}

	// get all the matches against the fragment
	matches, err := blast(&targetFrag, input.dbs, input.filters, conf.PCRMinLength, 100, targetFrag.Seq)
	if err != nil {
		dbMessage := strings.Join(input.dbs, ", ")
		return Frag{}, nil, fmt.Errorf("failed to blast %s against the dbs %s: %v", targetFrag.ID, dbMessage, err)
	}

	// for _, m := range matches {
	// 	fmt.Println(m.entry, m.start, m.end)
	// }

	// map fragment Matches to nodes
	var nodes []*Frag
	for _, m := range matches {
		nodes = append(nodes, newFrag(m, len(targetFrag.Seq), conf))
	}

	// build up a slice of assemblies that could, within the upper-limit on
	// fragment count, be assembled to make the target vector
	assemblies := createAssemblies(nodes, targetFrag.Seq, conf)

	fmt.Printf("%d assemblies created\n", len(assemblies))

	// build up a map from fragment count to a sorted list of assemblies with that number
	groupedAssemblies := groupAssemblies(assemblies)

	// sort the keys of the assemblies (their fragment count). Start with the smaller
	// assemblies first and work our way up
	assemblyCounts := []int{}
	for count := range groupedAssemblies {
		assemblyCounts = append(assemblyCounts, count)
	}
	sort.Ints(assemblyCounts)

	// "fill-in" the fragments (create synthetic fragments and primers)
	// skip assemblies that wouldn't be less than the current cheapest assembly
	filled := make(map[int][]*Frag)

	// append a fully synthetic solution at first, nothing added should cost more than this
	minCostAssembly := addFullySyntheticVector(filled, targetFrag.Seq, conf)

	for _, count := range assemblyCounts {
		// get the first assembly that fills properly (cheapest workable solution)
		for _, testAssembly := range groupedAssemblies[count] {
			if testAssembly.cost > minCostAssembly {
				// skip this and the rest with this count, there's another
				// cheaper option with fewer fragments
				break
			}

			filledFragments, err := testAssembly.fill(targetFrag.Seq, conf)
			if err != nil || filledFragments == nil {
				// write the console for debugging, continue looking
				// fmt.Println(testAssembly.log(), "error", err.Error())
				continue
			}

			// fmt.Println(testAssembly.log(), fragsCost(filledFragments))

			// if a Frag in the assembly fails to be prepared,
			// remove all assemblies with the Frag and try again
			newAssemblyCost := fragsCost(filledFragments)

			if newAssemblyCost >= minCostAssembly {
				continue // wasn't actually cheaper, keep trying
			}

			// store this as the new cheapest assembly
			minCostAssembly = newAssemblyCost

			// check whether there's already a cheaper solution
			if existingCost := fragsCost(filled[len(filledFragments)]); existingCost > 0 && existingCost < newAssemblyCost {
				continue
			}

			// set this is as the new cheapest of this length
			filled[len(filledFragments)] = filledFragments

			// delete all assemblies with more fragments that cost more
			for filledCount, existingFilledFragments := range filled {
				if filledCount <= len(filledFragments) {
					continue
				}

				existingCost := fragsCost(existingFilledFragments)

				if existingCost >= newAssemblyCost {
					delete(filled, filledCount)
				}
			}

			// break, don't look at other possible assemblies because this one worked
			break
		}
	}

	var found [][]*Frag
	for _, frags := range filled {
		found = append(found, frags)
	}

	return targetFrag, found, nil
}

// createAssemblies builds up circular assemblies (unfilled lists of fragments that should be combinable)
//
// maxNodes is the maximum number of fragments in a single assembly
// target is the target sequence we're trying to createAssemblies up
//
// It is created by traversing a DAG in forward order:
// 	foreach thisFragment (sorted in increasing start index order):
// 	  foreach otherFragment that thisFragment overlaps with + synthCount more:
//	 	foreach assembly on thisFragment:
//    	    add otherFragment to the assembly to create a new assembly, store on otherFragment
func createAssemblies(frags []*Frag, target string, conf *config.Config) (assemblies []assembly) {
	// number of additional frags try synthesizing to, in addition to those that
	// already have enough homology for overlap without any modifications for each Frag
	maxNodes := conf.FragmentsMaxCount
	synthCount := int(math.Max(5, 0.05*float64(len(frags)))) // 5 of 5%, whichever is greater

	// create a starting assembly on each Frag including just itself
	for i, f := range frags {
		// edge case where the Frag spans the entire target vector... 100% match
		// it is the target vector. just return that as the assembly
		if len(f.Seq) >= len(target) {
			return []assembly{
				assembly{
					frags:  []*Frag{f.copy()},
					synths: 0,
				},
			}
		}

		// create a starting assembly for each fragment just containing it
		frags[i].assemblies = []assembly{
			assembly{
				frags:  []*Frag{f.copy()}, // just self
				cost:   f.costTo(f),       // just PCR,
				synths: 0,                 // no synthetic frags at start
			},
		}
	}

	// for every Frag in the list of increasing start index frags
	for i, f := range frags {
		// for every overlapping fragment + synthCount more
		for _, j := range f.reach(frags, i, synthCount) {
			// for every assembly on the reaching fragment
			for _, a := range f.assemblies {
				newAssembly, created, circularized := a.add(frags[j], maxNodes)

				// if a new assembly wasn't created, move on
				if !created {
					continue
				}

				if circularized {
					// we've circularized a plasmid assembly, it's ready for filling
					assemblies = append(assemblies, newAssembly)
				} else {
					// add to the other fragment's list of assemblies
					frags[j].assemblies = append(frags[j].assemblies, newAssembly)
				}
			}
		}
	}

	return assemblies
}

// groupAssemblies returns a map from the number of fragments in a build
// to a slice of builds with that number of fragments, sorted by their cost
func groupAssemblies(assemblies []assembly) map[int][]assembly {
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

// addFullySyntheticVector adds a new fully synthetic vector to the built map if it's cheaper
// that any other solution of that many fragments.
// Returns the cost of synthesis for the fragments
func addFullySyntheticVector(built map[int][]*Frag, seq string, conf *config.Config) float64 {
	start := &Frag{start: 0, end: 0, conf: conf}
	end := &Frag{start: len(seq), end: len(seq), conf: conf}

	syntheticFrags := start.synthTo(end, seq)
	fCount := len(syntheticFrags)

	if _, filled := built[fCount]; filled {
		syntheticCost := fragsCost(syntheticFrags)
		existingCost := fragsCost(built[fCount])

		if syntheticCost < existingCost {
			built[fCount] = syntheticFrags
		}

		return syntheticCost // return the total cost of synthesis
	}
	built[fCount] = syntheticFrags

	return fragsCost(syntheticFrags)
}
