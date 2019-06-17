package plade

import (
	"fmt"
	"math"
	"sort"
	"strings"

	"github.com/jjtimmons/plade/config"
)

// assembly is a slice of nodes ordered by the nodes
// distance from the end of the target vector.
type assembly struct {
	// frags, ordered by distance from the "end" of the vector
	frags []*Frag

	// estimated cost of making this assembly
	cost float64

	// total number of synthetic nodes that will be needed to make this
	synths int
}

// add Frag to the end of an assembly. Return a new assembly and whether it circularized
func (a *assembly) add(f *Frag, maxCount, targetLength int, features bool) (newAssembly assembly, created, circularized bool) {
	firstStart := a.frags[0].start
	start := f.start
	end := f.end
	lastEnd := a.frags[len(a.frags)-1].end

	if features {
		firstStart = a.frags[0].featureStart
		start = f.featureStart
		end = f.featureEnd
		lastEnd = a.frags[len(a.frags)-1].featureEnd
	}

	// check if we could complete an assembly with this new Frag
	circularized = end >= firstStart+targetLength-1

	// check if this is the first fragment annealing to itself
	selfAnnealing := f.uniqueID == a.frags[0].uniqueID

	// the last fragment in the assembly
	last := a.frags[len(a.frags)-1]

	// calc the number of synthesis fragments needed to get to this next Frag
	synths := last.synthDist(f)
	if features && start > lastEnd {
		synths = start - lastEnd - 1
	}

	newCount := a.len() + synths
	if !selfAnnealing {
		newCount++
	}

	assemblyEnd := lastEnd
	if newCount > maxCount || (end-assemblyEnd < f.conf.PCRMinLength && !features) {
		return assembly{}, false, false
	}

	created = true

	// calc the estimated dollar cost of getting to the next Frag
	annealCost := last.costTo(f)
	if features {
		annealCost = 10.0
	}
	if selfAnnealing && synths == 0 {
		annealCost = 0 // does not cost extra to anneal to the first fragment
	}

	// check whether the Frag is already contained in the assembly
	// if so, the cost of procurement is not incurred twice
	nodeContained := false
	for _, included := range a.frags {
		if included.ID == f.ID && included.fragType == f.fragType {
			nodeContained = true
			break
		}
	}

	if nodeContained {
		// don't double count the cost of procuring this Frag to the total assembly cost
		annealCost += f.cost(false)
	} else {
		annealCost += f.cost(true)
	}

	// copy over all the fragments, need to avoid referencing same frags
	newFrags := []*Frag{}
	for _, frag := range a.frags {
		newFrags = append(newFrags, frag.copy())
	}
	if !selfAnnealing {
		newFrags = append(newFrags, f.copy())
	}

	return assembly{
		frags:  newFrags,
		cost:   a.cost + annealCost,
		synths: a.synths + synths,
	}, created, circularized
}

// len returns len(assembly.nodes) + the synthesis fragment count.
func (a *assembly) len() int {
	return len(a.frags) + a.synths
}

// log logs a description of the assembly (the entires in it and its cost).
func (a *assembly) log() {
	logString := ""
	if len(a.frags) >= 1 {
		for _, f := range a.frags {
			logString += f.ID + ":" + f.fragType.String() + " "
		}
	}

	fmt.Println(fmt.Sprintf("%s- $%.2f", logString, a.cost))
}

// fill traverses frags in an assembly and adds primers or makes syntheic fragments where necessary.
// It can fail. For example, a PCR Frag may have off-targets in the parent vector.
func (a *assembly) fill(target string, conf *config.Config) (frags []*Frag, err error) {
	min := conf.FragmentsMinHomology
	max := conf.FragmentsMaxHomology

	// check for and error out if there are duplicate ends between fragments,
	// ie unintended junctions between fragments that shouldn't be annealing
	if hasDuplicate, left, right, dupSeq := a.duplicates(a.frags, min, max); hasDuplicate {
		return nil, fmt.Errorf("duplicate junction between %s and %s: %s", left, right, dupSeq)
	}

	// edge case where a single Frag fills the whole target vector. Return just a single
	// "fragment" (of circular type... it is misnomer) that matches the target sequence 100%
	if a.len() == 1 && len(a.frags[0].Seq) >= len(target) {
		f := a.frags[0]

		return []*Frag{
			&Frag{
				ID:       f.ID,
				Seq:      strings.ToUpper(f.Seq)[0:len(target)], // it may be longer
				fragType: circular,
				URL:      f.URL,
				conf:     conf,
			},
		}, nil
	}

	// copy all the fragments. needed because ranges are mutated in assembly.fill,
	// so distance to neightbor estimates become invalid after a neighbor is mutated
	var origFrags []*Frag
	for _, f := range a.frags {
		origFrags = append(origFrags, f.copy())
	}

	// fill in primers. let each Frag create primers for itself that
	// will span it to the last and next fragments (if reachable)
	for i, f := range a.frags {
		// try and make primers for the fragment (need last and next nodes)
		var last *Frag
		if i == 0 {
			// mock up a last fragment that's to the left of this starting Frag
			last = &Frag{
				start: origFrags[len(origFrags)-1].start - len(target),
				end:   origFrags[len(origFrags)-1].end - len(target),
				conf:  conf,
			}
		} else {
			last = origFrags[i-1]
		}

		next := a.mockNext(origFrags, i, target, conf)

		// create primers for the Frag and add them to the Frag if it needs them
		// to anneal to the adjacent fragments
		lastPCR := !last.overlapsViaHomology(f) && last.overlapsViaPCR(f)
		nextPCR := !f.overlapsViaHomology(next) && f.overlapsViaPCR(next)
		needsPCR := f.fragType == circular || f.fragType == pcr || lastPCR || nextPCR

		// if the Frag has a full target from upload or
		if needsPCR {
			if err := f.setPrimers(last, next, target, conf); err != nil || len(f.Primers) < 2 {
				return nil, fmt.Errorf("failed to pcr %s: %v", f.ID, err)
			}
			f.fragType = pcr // is now a pcr type
		}

		// accumulate the prepared fragment
		frags = append(frags, f)
	}

	// second loop to fill in gaps between fragments that need to be filled via synthesis
	fragsWithSynth := []*Frag{}
	for i, f := range frags {
		if f.Seq != "" {
			fragsWithSynth = append(fragsWithSynth, f)
		}

		// add synthesized fragments between this Frag and the next (if necessary)
		next := a.mockNext(frags, i, target, conf)
		if synthedFrags := f.synthTo(next, target); synthedFrags != nil {
			fragsWithSynth = append(fragsWithSynth, synthedFrags...)
		}
	}
	frags = fragsWithSynth

	// validate that fragments will anneal to one another
	if err := validateJunctions(frags, conf); err != nil {
		return nil, err
	}

	return frags, nil
}

// mockNext returns the fragment that's one beyond the one passed.
// If there is none, it mocks one using the first fragment and changing
// its start and end index.
func (a *assembly) mockNext(frags []*Frag, i int, target string, conf *config.Config) *Frag {
	if i < len(frags)-1 {
		return frags[i+1]
	}

	// mock up a next fragment that's to the right of this terminal Frag
	return &Frag{
		ID:    frags[0].ID,
		start: frags[0].start + len(target),
		end:   frags[0].end + len(target),
		conf:  conf,
	}
}

// duplicates runs through all the nodes in an assembly and checks whether any of
// them have unintended homology, or "duplicate homology".
func (a *assembly) duplicates(frags []*Frag, min, max int) (isDup bool, first, second, dup string) {
	c := len(frags) // Frag count
	for i, f := range frags {
		// check to make sure the fragment doesn't anneal to itself
		if c > 1 {
			if selfJ := f.junction(f, min, max); selfJ != "" && len(selfJ) < len(f.Seq) {
				return true, f.ID, f.ID, selfJ
			}
		}

		for j := 2; j < c; j++ { // skip next Frag, i+1 is supposed to anneal to i
			junc := f.junction(frags[(j+i)%c], min, max)
			if junc != "" {
				return true, f.ID, frags[(j+i)%c].ID, junc
			}
		}
	}

	return false, "", "", ""
}

// createAssemblies builds up circular assemblies (unfilled lists of fragments that should be combinable)
//
// It is created by traversing a DAG in forward order:
//
// foreach fragment (sorted in increasing start index order):
//   foreach otherFragment that fragment overlaps with + reachSynthCount more:
//	   foreach assembly on fragment:
//       add otherFragment to the assembly to create a new assembly, store on otherFragment
func createAssemblies(frags []*Frag, vectorLength, targetLength int, features bool, conf *config.Config) (assemblies []assembly) {
	// number of additional frags try synthesizing to, in addition to those that
	// already have enough homology for overlap without any modifications for each Frag
	maxNodes := conf.FragmentsMaxCount

	// sort by start index again
	sort.Slice(frags, func(i, j int) bool {
		return frags[i].start < frags[j].start
	})

	// create a starting assembly on each Frag including just itself
	for i, f := range frags {
		// edge case where the Frag spans the entire target vector... 100% match
		// it is the target vector. just return that as the assembly
		if len(f.Seq) >= vectorLength {
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

	for i, f := range frags { // for every Frag in the list of increasing start index frags
		for _, j := range f.reach(frags, i) { // for every overlapping fragment + reach more
			for _, a := range f.assemblies { // for every assembly on the reaching fragment
				newAssembly, created, circularized := a.add(frags[j], maxNodes, targetLength, features)

				if !created { // if a new assembly wasn't created, move on
					continue
				}

				if circularized { // we've circularized a vector, it's ready for filling
					assemblies = append(assemblies, newAssembly)
				} else {
					// add to the other fragment's list of assemblies
					frags[j].assemblies = append(frags[j].assemblies, newAssembly)
				}
			}
		}
	}

	if conf.Verbose {
		fmt.Printf("%d assemblies made\n", len(assemblies))
	}

	return assemblies
}

// groupAssembliesByCount returns a map from the number of fragments in a build
// to a slice of builds with that number of fragments, sorted by their cost.
func groupAssembliesByCount(assemblies []assembly) ([]int, map[int][]assembly) {
	countToAssemblies := make(map[int][]assembly)
	for _, a := range assemblies {
		if as, ok := countToAssemblies[a.len()]; ok {
			countToAssemblies[a.len()] = append(as, a)
		} else {
			countToAssemblies[a.len()] = []assembly{a}
		}
	}

	// sort the fragment counts of assemblies and the assemblies within each
	// assembly count, so we're trying the shortest assemblies first, and the cheapest
	// assembly within each fragment count before the others
	var counts []int
	for count := range countToAssemblies {
		counts = append(counts, count)
		sort.Slice(countToAssemblies[count], func(i, j int) bool {
			return countToAssemblies[count][i].cost < countToAssemblies[count][j].cost
		})
	}
	sort.Ints(counts)

	return counts, countToAssemblies
}

// fillAssemblies fills in assemblies and returns the pareto optimal solutions.
func fillAssemblies(target string, counts []int, countToAssemblies map[int][]assembly, conf *config.Config) (solutions [][]*Frag) {
	// append a fully synthetic solution at first, nothing added should cost more than this (single vector)
	filled := make(map[int][]*Frag)
	minCostAssembly := math.MaxFloat64

	for _, count := range counts {
		for _, assemblyToFill := range countToAssemblies[count] {
			if assemblyToFill.cost > minCostAssembly {
				// skip this and the rest with this count, there's another
				// cheaper option with the same number or fewer fragments (estimated)
				// fmt.Println("too expensive break")
				break
			}

			filledFragments, err := assemblyToFill.fill(target, conf)
			if err != nil || filledFragments == nil {
				// assemblyToFill.log()
				// fmt.Println("error", err.Error())
				continue
			}

			newAssemblyCost := fragsCost(filledFragments)

			if newAssemblyCost >= minCostAssembly || len(filledFragments) > conf.FragmentsMaxCount {
				// fmt.Println("too expensive continue")
				continue // wasn't actually cheaper, keep trying
			}
			minCostAssembly = newAssemblyCost // store this as the new cheapest assembly

			// delete all assemblies with more fragments that cost more
			for filledCount, existingFilledFragments := range filled {
				if filledCount < len(filledFragments) {
					continue
				}

				existingCost := fragsCost(existingFilledFragments)
				if existingCost >= newAssemblyCost {
					delete(filled, filledCount)
				}
			}

			// set this is as the new cheapest of this length
			filled[len(filledFragments)] = filledFragments
		}
	}

	for _, frags := range filled {
		solutions = append(solutions, frags) // flatten
	}

	return solutions
}
