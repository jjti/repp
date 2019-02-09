package defrag

import (
	"fmt"
	"os"
	"strings"
	"time"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
)

// SequenceCmd takes a cobra command (with its flags) and runs Vector
func SequenceCmd(cmd *cobra.Command, args []string) {
	Sequence(parseCmdFlags(cmd, args))
}

// Sequence is for running an end to end vector design using a target sequence
func Sequence(flags *Flags, conf *config.Config) {
	handleErr := func(err error) {
		if err != nil {
			stderr.Fatalf("failed to assemble the vector sequence in %s: %v", flags.in, err)
		}
	}
	start := time.Now()

	target, builds, err := sequence(flags, conf) // build up the assemblies that make the sequence
	handleErr(err)

	// write the results to a file
	elapsed := time.Since(start)
	_, err = write(flags.out, target.ID, target.Seq, builds, len(target.Seq), conf, elapsed.Seconds())
	handleErr(err)

	fmt.Printf("%s\n\n", elapsed)
	os.Exit(0)
}

// sequence builds a vector using a simple cost optimization scheme
//
// the goal is to find an "optimal" assembly sequence with:
// 	1. the fewest fragments
// 	2. the lowest overall assembly cost ($)
// and, secondarily:
//	3. no duplicate end regions between Gibson fragments
// 	4. no hairpins in the junctions
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
func sequence(input *Flags, conf *config.Config) (Frag, [][]*Frag, error) {
	// read the target sequence (the first in the slice is used)
	fragments, err := read(input.in, false)
	if err != nil {
		return Frag{}, nil, fmt.Errorf("failed to read fragments from %s: %v", input.in, err)
	}

	if len(fragments) > 1 {
		stderr.Printf(
			"warning: %d fragments were in %s. Only targeting the sequence of the first: %s\n",
			len(fragments),
			input.in,
			fragments[0].ID,
		)
	}

	target := fragments[0]
	fmt.Printf("Building %s\n", target.ID)

	// if a backbone was specified, add it to the sequence of the target frag
	if input.backbone.ID != "" {
		target.Seq += input.backbone.Seq
	}

	// get all the matches against the target vector
	tw := blastWriter()
	matches, err := blast(target.ID, target.Seq, true, input.dbs, input.filters, 100, tw)
	tw.Flush()

	if err != nil {
		dbMessage := strings.Join(input.dbs, ", ")
		return Frag{}, nil, fmt.Errorf("failed to blast %s against the dbs %s: %v", target.ID, dbMessage, err)
	}
	if len(matches) < 1 {
		return Frag{}, nil, fmt.Errorf("did not find any matches for %s", target.ID)
	}

	// keep only "proper" arcs (non-self-contained)
	matches = filter(matches, len(target.Seq), conf.PCRMinLength)
	fmt.Printf("%d matches after filtering\n", len(matches))

	// map fragment Matches to nodes
	frags := newFrags(matches, conf)

	// build up a slice of assemblies that could, within the upper-limit on
	// fragment count, be assembled to make the target vector
	assemblies := createAssemblies(frags, len(target.Seq), conf, cmdSequence)

	// build up a map from fragment count to a sorted list of assemblies with that number
	assemblyCounts, countToAssemblies := groupAssembliesByCount(assemblies)

	// append a fully synthetic solution at first, nothing added should cost more than this (single vector)
	filled := make(map[int][]*Frag)
	minCostAssembly := addSyntheticVector(filled, target.Seq, conf)

	for _, count := range assemblyCounts {
		for _, assemblyToFill := range countToAssemblies[count] {
			if assemblyToFill.cost > minCostAssembly {
				// skip this and the rest with this count, there's another
				// cheaper option with fewer fragments
				break
			}

			filledFragments, err := assemblyToFill.fill(target.Seq, conf)
			if err != nil || filledFragments == nil {
				// write the console for debugging, continue looking
				// logger.Println(assemblyToFill.log(), "error", err.Error())
				continue
			}
			// fmt.Println(assemblyToFill.log(), fragsCost(filledFragments))

			// if a Frag in the assembly fails to be prepared,
			// remove all assemblies with the Frag and try again
			newAssemblyCost := fragsCost(filledFragments)

			if newAssemblyCost >= minCostAssembly {
				continue // wasn't actually cheaper, keep trying
			}

			// store this as the new cheapest assembly
			minCostAssembly = newAssemblyCost

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

			// don't look at other possible assemblies, assume this will be the cheapest of this length
			break
		}
	}

	var solutions [][]*Frag
	for _, frags := range filled {
		solutions = append(solutions, frags) // flatten
	}

	return target, solutions, nil
}

// addSyntheticVector adds a new fully synthetic vector to the built map if it's cheaper
// that any other solution of that many fragments.
// Returns the cost of synthesis for the fragments
func addSyntheticVector(built map[int][]*Frag, seq string, conf *config.Config) float64 {
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
