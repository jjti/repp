package defrag

import (
	"fmt"
	"math"
	"sort"
	"strings"

	"github.com/jjtimmons/defrag/config"
)

// assembly is a slice of nodes ordered by the nodes
// distance from the end of the target vector
type assembly struct {
	// frags, ordered by distance from the "end" of the vector
	frags []*Frag

	// estimated cost of making this assembly
	cost float64

	// total number of synthetic nodes that will be needed to make this
	synths int
}

// add adds Frag to the end of a sequence assembly. Used when building up a target sequence.
func (a *assembly) add(f *Frag, maxCount, targetLength int) (newAssembly assembly, created, circularized bool) {
	// check if we could complete an assembly with this new Frag
	circularized = f.end > a.frags[0].start+targetLength

	// check if this is the first fragment annealing to itself
	selfAnnealing := f.uniqueID == a.frags[0].uniqueID

	// calc the number of synthesis fragments needed to get to this next Frag
	synths := a.frags[len(a.frags)-1].synthDist(f)
	newCount := a.len() + synths
	if !selfAnnealing {
		newCount++
	}

	if newCount > maxCount {
		return assembly{}, false, false
	}

	created = true

	// calc the estimated dollar cost of getting to the next Frag
	annealCost := a.frags[len(a.frags)-1].costTo(f)
	if selfAnnealing && synths == 0 {
		annealCost = 0 // does not cost extra to anneal to the first fragment
	}

	// check whether the Frag is already contained in the assembly
	// if so, the cost of procurement is not incurred twice
	nodeContained := false
	for _, included := range a.frags {
		if included.ID == f.ID {
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

// len returns len(assembly.nodes) + the synthesis fragment count
func (a *assembly) len() int {
	return len(a.frags) + a.synths
}

// log returns a description of the assembly (the entires in it and its cost)
func (a *assembly) log() string {
	logString := ""
	if len(a.frags) >= 1 {
		for _, f := range a.frags {
			logString += f.ID + " "
		}
	}

	return fmt.Sprintf("%s- $%.2f", logString, a.cost)
}

// fill traverses frags in an assembly and adds primers or makes syntheic fragments where necessary
//
// - Vector fragments if a Frag spans the whole target range (entirely contains it)
//
// - PCR fragments if the Frag subselects a region (BLAST match)
//
// - Synthetic fragments if there wasn't a match for a region
//
// It can fail. For example, a PCR Frag may have off-targets in the parent vector
func (a *assembly) fill(target string, conf *config.Config) (frags []*Frag, err error) {
	minHomology := conf.FragmentsMinHomology
	maxHomology := conf.FragmentsMaxHomology

	// check for and error out if there are duplicate ends between fragments,
	// ie unintended junctions between fragments that shouldn't be annealing
	if hasDuplicate, left, right, dupSeq := a.duplicates(a.frags, minHomology, maxHomology); hasDuplicate {
		return nil, fmt.Errorf("failed to fill: duplicate junction sequence in %s and %s: %s", left, right, dupSeq)
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

	// do two loops. the first is to fill in primers. let each Frag create primers for
	// itself that will span it to the last and next fragments (if reachable). this has to be done
	// in two loops because synthesis depends on nodes' ranges, and nodes' ranges
	// may change during setPrimers (we let primer3_core pick from a range of primer options)
	for i, f := range a.frags {
		// try and make primers for the fragment (need last and next nodes)
		var last *Frag
		if i == 0 {
			// mock up a last fragment that's to the left of this starting Frag
			last = &Frag{
				start: a.frags[len(a.frags)-1].start - len(target),
				end:   a.frags[len(a.frags)-1].end - len(target),
				conf:  conf,
			}
		} else {
			last = a.frags[i-1]
		}

		next := a.mockNext(i, target, conf)

		// create primers for the Frag and add them to the Frag if it needs them
		// to anneal to the adjacent fragments
		needsPCR := f.fragType == circular ||
			f.fullSeq == "" ||
			!last.overlapsViaHomology(f) && last.overlapsViaPCR(f) ||
			!f.overlapsViaHomology(next) && f.overlapsViaPCR(next)

		// if the Frag has a full target from upload or
		if needsPCR {
			if err := f.setPrimers(last, next, target, conf); err != nil || len(f.Primers) < 2 {
				return nil, fmt.Errorf("failed to pcr %s: %v", f.ID, err)
			}
			f.fragType = pcr // is now a pcr type
		} else {
			f.fragType = existing // no change needed
		}

		// accumulate the prepared fragment
		frags = append(frags, f)

		// add synthesized fragments between this Frag and the next (if necessary)
		if synthedFrags := f.synthTo(a.mockNext(i, target, conf), target); synthedFrags != nil {
			frags = append(frags, synthedFrags...)
		}
	}

	return
}

// mockNext returns the fragment that's one beyond the one passed
// if there is none, it mocks one using the first fragment and changing
// its start and end index
func (a *assembly) mockNext(i int, target string, conf *config.Config) *Frag {
	if i < len(a.frags)-1 {
		return a.frags[i+1]
	}

	// mock up a next fragment that's to the right of this terminal Frag
	return &Frag{
		start: a.frags[0].start + len(target),
		end:   a.frags[0].end + len(target),
		conf:  conf,
	}
}

// duplicates runs through all the nodes in an assembly and checks whether any of
// them have unintended homology, or "duplicate homology"
func (a *assembly) duplicates(nodes []*Frag, minHomology, maxHomology int) (isDup bool, first, second, dup string) {
	c := len(nodes) // Frag count
	for i, f := range nodes {
		for j := 2; j <= c; j++ { // skip next Frag, i+1 is supposed to anneal to i
			junc := f.junction(nodes[(j+i)%c], minHomology, maxHomology)
			if junc != "" && f.uniqueID != nodes[(j+i)%c].uniqueID {
				return true, f.ID, nodes[(j+i)%c].ID, junc
			}
		}
	}

	return false, "", "", ""
}

// createAssemblies builds up circular assemblies (unfilled lists of fragments that should be combinable)
//
// maxNodes is the maximum number of fragments in a single assembly
// target is the target sequence we're trying to createAssemblies up
//
// It is created by traversing a DAG in forward order:
// 	foreach thisFragment (sorted in increasing start index order):
// 	  foreach otherFragment that thisFragment overlaps with + reachSynthCount more:
//	 	foreach assembly on thisFragment:
//    	    add otherFragment to the assembly to create a new assembly, store on otherFragment
func createAssemblies(frags []*Frag, targetLength int, conf *config.Config) (assemblies []assembly) {
	// number of additional frags try synthesizing to, in addition to those that
	// already have enough homology for overlap without any modifications for each Frag
	maxNodes := conf.FragmentsMaxCount
	reachSynthCount := int(math.Max(5, 0.05*float64(len(frags)))) // 5 of 5%, whichever is greater

	// create a starting assembly on each Frag including just itself
	for i, f := range frags {
		// edge case where the Frag spans the entire target vector... 100% match
		// it is the target vector. just return that as the assembly
		if len(f.Seq) >= targetLength {
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
		// for every overlapping fragment + reachSynthCount more
		for _, j := range f.reach(frags, i, reachSynthCount) {
			// for every assembly on the reaching fragment
			for _, a := range f.assemblies {
				newAssembly, created, circularized := a.add(frags[j], maxNodes, targetLength)

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

	fmt.Printf("%d assemblies made\n", len(assemblies))

	return assemblies
}

// groupAssembliesByCount returns a map from the number of fragments in a build
// to a slice of builds with that number of fragments, sorted by their cost
func groupAssembliesByCount(assemblies []assembly) ([]int, map[int][]assembly) {
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

	// sort the fragment counts of assemblies and the assemblies within each
	// assembly count, so we're trying the shortest assemblies first, and the cheapest
	// assembly within each fragment count before the others
	assemblyCounts := []int{}
	for count := range countToAssemblies {
		assemblyCounts = append(assemblyCounts, count)
		sort.Slice(countToAssemblies[count], func(i, j int) bool {
			return countToAssemblies[count][i].cost < countToAssemblies[count][j].cost
		})
	}
	sort.Ints(assemblyCounts)

	return assemblyCounts, countToAssemblies
}
