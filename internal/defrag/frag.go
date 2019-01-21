package defrag

import (
	"fmt"
	"math"
	"regexp"
	"strconv"
	"strings"

	"github.com/jjtimmons/defrag/config"
)

// Type is the Frag building type to be used in the assembly
type Type int

const (
	// circular is a circular sequence of DNA, e.g.: many of Addgene's plasmids
	circular Type = 0

	// pcr fragments are those prepared by pcr, often a subselection of their parent vector
	pcr Type = 1

	// synthetic fragments are those that will be fully synthesized (ex: gBlocks)
	synthetic Type = 2

	// linear fragment, ie the type of a fragment as it was uploaded submitted and without PCR/synthesis
	existing Type = 3
)

// Frag is a single building block stretch of DNA for assembly
type Frag struct {
	// ID is a unique identifier for this fragment
	ID string `json:"-"`

	// uniqueID of a match, the start index % seq-length + ID
	// used identified that catches nodes that cross the zero-index
	uniqueID string

	// URL, eg link to a vector's addgene page
	URL string `json:"url,omitempty"`

	// Cost to make the fragment
	Cost float64 `json:"costDollars"`

	// fragment's sequence (linear)
	Seq string `json:"seq,omitempty"`

	// fullSeq is the entire seq of the Frag/fragment as it was read in (for forward engineering)
	fullSeq string

	// db that the frag came from
	db string

	// circular if it was submitted/created as a circular fragment (marked as such in the DB)
	circular bool

	// start of this Frag on the target vector (which has been 3x'ed for BLAST)
	start int

	// end of this Frag on the target vector
	end int

	// assemblies that span from this Frag to the end of the vector
	assemblies []assembly

	// primers necessary to create this (if pcr fragment)
	Primers []Primer `json:"primers,omitempty"`

	// Entry of this fragment In the DB that it came from
	// Used to look for off-targets
	Entry string `json:"-"`

	// Type of this fragment
	Type Type `json:"-"`

	// Assembly configuration
	conf *config.Config
}

// newFrag creates a Frag from a match
func newFrag(m match, seqL int, conf *config.Config) *Frag {
	cost := 0.0
	url := ""
	if strings.Contains(m.entry, "gnl|addgene") {
		// create a link to the source Addgene page
		re := regexp.MustCompile("^.*addgene\\|(\\d*)")
		match := re.FindStringSubmatch(m.entry)

		cost = conf.AddGeneVectorCost // cost of addgene procurement
		url = fmt.Sprintf("https://www.addgene.org/%s/", match[1])
	}
	if strings.Contains(m.entry, "gnl|igem") {
		// create a source to the source iGEM page
		re := regexp.MustCompile("^.*igem\\|(\\w*)")
		match := re.FindStringSubmatch(m.entry)

		cost = conf.IGEMPartCost // cost of igem part procurement
		url = fmt.Sprintf("http://parts.igem.org/Part:%s", match[1])
	}

	return &Frag{
		ID:       m.entry,
		uniqueID: strconv.Itoa(m.start%seqL) + m.entry,
		Seq:      strings.ToUpper(m.seq),
		circular: m.circular,
		start:    m.start,
		end:      m.end,
		db:       m.db,
		Cost:     cost,
		URL:      url,
		conf:     conf,
	}
}

// copy returns a deep dopy of a Frag. used because nodes are mutated
// during assembly filling, and we don't want primers being shared between
// nodes in different assemblies
func (f *Frag) copy() *Frag {
	return &Frag{
		ID:       f.ID,
		Seq:      f.Seq,
		uniqueID: f.uniqueID,
		start:    f.start,
		end:      f.end,
		db:       f.db,
		conf:     f.conf,
		Cost:     f.Cost,
		URL:      f.URL,
		Type:     f.Type,
	}
}

// distTo returns the distance between the start of this Frag and the end of the other.
// assumes that this Frag starts before the other
// will return a negative number if this Frag overlaps with the other and positive otherwise
func (f *Frag) distTo(other *Frag) (bpDist int) {
	return other.start - f.end
}

// synthDist returns the number of synthesized fragments that would need to be created
// between one Frag and another if the two were to be joined, with no existing
// fragments/nodes in-between, in an assembly
func (f *Frag) synthDist(other *Frag) (synthCount int) {
	dist := f.distTo(other)

	if dist <= f.conf.PCR.MaxEmbedLength {
		// if the dist is <MaxEmbedLength, we can PCR our way there
		// and add the mutated bp between the nodes with PCR
		return 0
	}

	// if it's > 5, we have to synthesize to get there (arbitatry, TODO: move to settings)
	// can't be negative before the ceil
	floatDist := math.Max(1.0, float64(dist))

	// split up the distance between them by the max synthesized fragment size
	return int(math.Ceil(floatDist / float64(f.conf.Synthesis.MaxLength)))
}

// costTo calculates the $ amount needed to get from this fragment
// to the other Frag passed, either by PCR or synthesis
//
// PCR is, of course, preferred if we can add homology between this and
// the other fragment with PCR homology arms alone
//
// Otherwise we find the total synthesis distance between this and
// the other fragment and divide that by the cost per bp of synthesized DNA
//
// This does not add in the cost of procurement, which is added to the assembly cost
// in assembly.add()
func (f *Frag) costTo(other *Frag) (cost float64) {
	dist := f.distTo(other)

	if dist <= f.conf.PCR.MaxEmbedLength {
		if dist < -(f.conf.Fragments.MinHomology) {
			// there's already enough overlap between this Frag and the one being tested
			// estimating two primers, 20bp each
			return 46 * f.conf.PCR.BPCost
		}

		// we have to create some additional primer sequence to reach the next fragment
		// guessing 46bp (2x primer3 target primer length) plus half MinHomology on each primer
		return float64(46+f.conf.Fragments.MinHomology) * f.conf.PCR.BPCost
	}

	// we need to create a new synthetic fragment to get from this fragment to the next
	// to account for both the bps between them as well as the additional bps we need to add
	// for homology between the two
	fragLength := f.conf.Fragments.MinHomology + dist
	return f.conf.SynthCost(fragLength)
}

// reach returns a slice of Frag indexes that overlap with, or are the first synth_count nodes
// away from this one within a slice of ordered nodes
func (f *Frag) reach(nodes []*Frag, i, synthCount int) (reachable []int) {
	reachable = []int{}

	// accumulate the nodes that overlap with this one
	for true {
		i++

		// we've run out of nodes
		if i >= len(nodes) {
			return reachable
		}

		// these nodes overlap by enough for assembly without PCR
		if f.distTo(nodes[i]) <= -(f.conf.Fragments.MinHomology) {
			reachable = append(reachable, i)
		} else if synthCount > 0 {
			// there's not enough existing overlap, but we can synthesize to it
			synthCount--
			reachable = append(reachable, i)
		} else {
			break
		}
	}
	return
}

// junction checks for and returns any 100% identical homology between the end of this
// Frag and the start of the other. returns an empty string if there's no junction between them
func (f *Frag) junction(other *Frag, minHomology, maxHomology int) (junction string) {
	thisSeq := strings.ToUpper(f.Seq)
	otherSeq := strings.ToUpper(other.Seq)

	//      v-maxHomology from end    v-minHomology from end
	// ------------------------------------
	//                    -----------------------------
	start := len(thisSeq) - maxHomology
	end := len(thisSeq) - minHomology

	if start < 0 {
		start = 0
	}

	// for every possible start index
	for i := start; i <= end; i++ {
		// traverse from that index to the end of the seq
		j := 0
		for k := i; k < len(thisSeq); k++ {
			if thisSeq[k] != otherSeq[j] {
				break
			} else {
				j++
			}

			// we made it to the end of the sequence, there's a junction
			if k == len(thisSeq)-1 {
				return thisSeq[i:]
			}
		}
	}
	return
}

// synthTo returns synthetic fragments to get this Frag to the next.
// It creates a slice of building fragments that have homology against
// one another and are within the upper and lower synthesis bounds.
// seq is the vector's sequence. We need it to build up the target
// vector's sequence
func (f *Frag) synthTo(next *Frag, seq string) (synthedFrags []Frag) {
	// check whether we need to make synthetic fragments to get
	// to the next fragment in the assembly
	fragC := f.synthDist(next)
	if fragC == 0 {
		return nil
	}

	// make the slice
	synthedFrags = []Frag{}

	// length of each synthesized fragment
	fragL := f.distTo(next) / fragC
	if f.conf.Synthesis.MinLength > fragL {
		// need to synthesize at least Synthesis.MinLength bps
		fragL = f.conf.Synthesis.MinLength
	}

	// account for homology on either end of each synthetic fragment
	fragL += f.conf.Fragments.MinHomology * 2
	seq += seq // double to account for sequence across the zero-index

	// slide along the range of sequence to create synthetic fragments
	// and create one at each point, each w/ MinHomology for the fragment
	// before and after it
	for fragIndex := 0; fragIndex < int(fragC); fragIndex++ {
		start := f.end - f.conf.Fragments.MinHomology // start w/ homology
		start += fragIndex * fragL                    // slide along the range to cover
		end := start + fragL + f.conf.Fragments.MinHomology

		sFrag := Frag{
			ID:   fmt.Sprintf("%s-synthetic-%d", f.ID, fragIndex+1),
			Seq:  strings.ToUpper(seq[start:end]),
			Cost: f.conf.SynthCost(len(f.Seq)) + f.Cost,
			Type: synthetic,
		}
		synthedFrags = append(synthedFrags, sFrag)
	}
	return
}
