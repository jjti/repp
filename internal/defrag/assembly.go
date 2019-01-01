package defrag

import (
	"fmt"
	"strings"

	"github.com/jjtimmons/defrag/config"
)

// assembly is a slice of nodes ordered by the nodes
// distance from the end of the target vector
type assembly struct {
	// nodes, ordered by distance from the "end" of the vector
	nodes []*Frag

	// estimated cost of making this assembly
	cost float64

	// total number of synthetic nodes that will be needed to make this
	synths int
}

// add a Frag to the start of this assembly.
//
// Update the cost of the assembly to include the link between the new first Frag and the one after it.
// Store the Frag's ID in the list of Frag ids.
// Complete  an assembly if a Frag has matched up onto itself across the zero-index.
func (a *assembly) add(f *Frag, maxCount int) (newAssembly assembly, created, complete bool) {
	// check if we could complete an assembly with this new Frag
	complete = f.uniqueID == a.nodes[0].uniqueID
	// calc the number of synthesis fragments needed to get to this next Frag
	synths := a.nodes[len(a.nodes)-1].synthDist(f)

	newCount := a.len() + synths
	if !complete {
		// only adding a new Frag if not annealing to starting Frag
		newCount++
	}

	// stay beneath upper limit
	if newCount > maxCount {
		return assembly{}, false, false
	}

	// we will create a new assembly
	created = true

	// calc the estimated dollar cost of getting to the next Frag
	annealCost := a.nodes[len(a.nodes)-1].costTo(f)

	// check whether the Frag is already contained in the assembly
	nodeContained := false
	for _, included := range a.nodes {
		if included.ID == f.ID {
			nodeContained = true
			break
		}
	}

	if !nodeContained {
		// add the cost of procuring this Frag to the total assembly cost
		annealCost += f.Cost
	}

	if complete {
		if synths < 1 {
			// costs nothing to anneal Frag to self, already been PCR'ed
			annealCost = 0
		}

		return assembly{
			nodes:  append(a.nodes, f.copy()),
			cost:   a.cost + annealCost,
			synths: a.synths + synths,
		}, created, complete
	}

	return assembly{
		nodes:  append(a.nodes, f.copy()),
		cost:   a.cost + annealCost,
		synths: a.synths + synths,
	}, created, false
}

// contains returns if the ID of the Frag has already been seen in this assembly
func (a *assembly) contains(f Frag) (hasNode bool) {
	for _, otherN := range a.nodes {
		// they're the same if they have the same ID and start index
		// ID isn't enough by itself because there may be multiple with the same
		// entry ID in the BLAST db
		if otherN.uniqueID == f.uniqueID {
			return true
		}
	}
	return false
}

// len returns len(assembly.nodes) + the synthesis fragment count
func (a *assembly) len() int {
	return len(a.nodes) + a.synths
}

// fill traverses the nodes in an assembly and converts them into fragments
//
// - Vector fragments if a Frag spans the whole target range (entirely contains it)
//
// - PCR fragments if the Frag subselects a region (BLAST match)
//
// - Synthetic fragments if there wasn't a match for a region
//
// It can fail. For example, a PCR Frag may have off-targets in the parent vector
func (a *assembly) fill(seq string, conf *config.Config) (frags []Frag, err error) {
	minHomology := conf.Fragments.MinHomology
	maxHomology := conf.Fragments.MaxHomology

	// check for and error out if there are duplicate ends between fragments
	// ie unintended junctions between fragments that shouldn't be annealing
	if hasDuplicate, duplicateSeq := a.duplicates(a.nodes, minHomology, maxHomology); hasDuplicate {
		return nil, fmt.Errorf("cannot fill, has duplicate junction sequence %s", duplicateSeq)
	}

	// edge case where a single Frag fills the whole target vector. Return just a single
	// "fragment" (of Vector type... misnomer) that matches the target sequence 100%
	if a.len() == 1 && len(a.nodes[0].Seq) >= len(seq) {
		f := a.nodes[0]
		return []Frag{
			Frag{
				ID:    f.ID,
				Seq:   strings.ToUpper(f.Seq)[0:len(seq)], // it may be longer
				Entry: f.ID,
				Type:  circular,
				URL:   f.URL,
				Cost:  f.Cost, // only the ordering cost, no PCR/Synth etc
			},
		}, nil
	}

	// do two loops. the first is to fill in primers. let each Frag create primers for
	// itself that will span it to the last and next fragments (if reachable). this has to be done
	// in two loops because synthesis depends on nodes' ranges, and nodes' ranges
	// may change during setPrimers (we let primer3_core pick from a range of primer options)
	for i, f := range a.nodes {
		// last Frag, do nothing
		// here only to allow for vector "circularization" if we need to synthesize
		// from a.nodes[len(a.nodes)-2] to a.nodes[len(a.nodes)-1]
		if i > 0 && f.uniqueID == a.nodes[0].uniqueID {
			break
		}

		// try and make primers for the fragment (need last and next nodes)
		var last *Frag
		if i == 0 {
			// mock up a last fragment that's to the left of this starting Frag
			final := a.nodes[len(a.nodes)-1]
			if f.uniqueID != "" && f.uniqueID == final.uniqueID {
				// -2 is if the first and last are the same and is just there for circularization
				final = a.nodes[len(a.nodes)-2]
			}

			last = &Frag{
				start: final.start - len(seq),
				end:   final.end - len(seq),
				conf:  conf,
			}
		} else {
			last = a.nodes[i-1]
		}

		var next *Frag
		if i < len(a.nodes)-1 {
			next = a.nodes[i+1]
		} else {
			// mock up a next fragment that's to the right of this terminal Frag
			first := a.nodes[0]
			next = &Frag{
				start: first.start + len(seq),
				end:   first.end + len(seq),
				conf:  conf,
			}
		}

		// create primers for the Frag and add them to the Frag if it needs them
		// to anneal to the adjacent fragments
		distLeft := last.distTo(f)
		distRight := f.distTo(next)

		// if the Frag has a full seq from upload and already has enough overlap with
		// the last and next nodes, we don't have to add anything to it (PCR it)
		if f.fullSeq == "" || (distLeft < -minHomology && distRight < -minHomology) {
			if err := f.setPrimers(last, next, seq, conf); err != nil || len(f.Primers) < 2 {
				return nil, fmt.Errorf("failed to fill %s: %v", f.ID, err)
			}
		}
	}

	// do another loop to covert the nodes with primers to fragments and
	// synthesize the sequence between nodes
	for i, f := range a.nodes {
		// last Frag, do nothing
		// here only to allow for vector "circularization" if we need to synthesize
		// from a.nodes[len(a.nodes)-2] to a.nodes[len(a.nodes)-1]
		if i > 0 && f.uniqueID != "" && f.uniqueID == a.nodes[0].uniqueID {
			break
		}

		// convert to a fragment from a Frag, store this to the list of building fragments
		// cost is calculated here as the summed cost of both primers (based on length)
		frags = append(frags, *f)

		// add synthesized fragments between this Frag and the next (if necessary)
		if synthedFrags := f.synthTo(a.nodes[(i+1)%len(a.nodes)], seq); synthedFrags != nil {
			frags = append(frags, synthedFrags...)
		}
	}
	return
}

// duplicates runs through all the nodes in an assembly and checks whether any of
// them have unintended homology, or "duplicate homology"
// Need to cheack each Frag's junction() with nodes other than the one after it
// (which it's supposed to anneal to)
func (a *assembly) duplicates(nodes []*Frag, minHomology, maxHomology int) (bool, string) {
	c := len(nodes) // Frag count
	for i, f := range nodes {
		for j := 2; j <= c; j++ { // skip next Frag, this is supposed to anneal to that
			if junc := f.junction(nodes[(j+i)%c], minHomology, maxHomology); junc != "" {
				return true, junc
			}
		}
	}
	return false, ""
}
