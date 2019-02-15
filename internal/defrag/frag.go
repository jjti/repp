package defrag

import (
	"fmt"
	"math"
	"os"
	"regexp"
	"strings"

	"github.com/jinzhu/copier"
	"github.com/jjtimmons/defrag/config"
)

// fragType is the Frag building type to be used in the assembly
type fragType int

const (
	// circular is a circular sequence of DNA, e.g.: many of Addgene's plasmids
	circular fragType = iota

	// pcr fragments are those prepared by pcr, often a subselection of their parent vector
	pcr

	// synthetic fragments are those that will be fully synthesized (eg: gBlocks)
	synthetic

	// linear fragment, ie the type of a fragment as it was uploaded submitted and without PCR/synthesis
	existing
)

// Frag is a single building block stretch of DNA for assembly
type Frag struct {
	// ID is a unique identifier for this fragment
	ID string `json:"-"`

	// type of the fragment in string representation for export
	Type string `json:"type"`

	// Cost to make the fragment
	Cost float64 `json:"cost"`

	// URL, eg link to a vector's addgene page
	URL string `json:"url,omitempty"`

	// fragment/vector's sequence
	Seq string `json:"seq,omitempty"`

	// sequence of a pcr fragment after PCR's addition of bp
	PCRSeq string `json:"pcrSeq,omitempty"`

	// primers necessary to create this (if pcr fragment)
	Primers []Primer `json:"primers,omitempty"`

	// fragType of this fragment. circular | pcr | synthetic | existing
	fragType fragType

	// uniqueID of a match, ID + the start index % seq-length
	// used identified that catches nodes that cross the zero-index
	uniqueID string

	// fullSeq is the entire seq of the Frag/fragment as it was read in (for forward engineering)
	fullSeq string

	// db that the frag came from
	db string

	// start of this Frag on the target vector
	start int

	// end of this Frag on the target vector
	end int

	// assemblies that span from this Frag to the end of the vector
	assemblies []assembly

	// build configuration
	conf *config.Config
}

// ranged: a stretch that a primer spans relative to the target sequence
type ranged struct {
	// start of the primer's range
	start int

	// end of the primer's range
	end int
}

// Primer is a single Primer used to create a PCR fragment
type Primer struct {
	// Seq of the primer (In 5' to 3' direction)
	Seq string `json:"seq"`

	// Strand of the primer; true if template, false if complement
	Strand bool `json:"strand"`

	// Penalty score
	Penalty float64 `json:"penalty"`

	// PairPenalty score from primer3
	PairPenalty float64 `json:"pairPenalty"`

	// Tm of the primer
	Tm float64 `json:"tm"`

	// GC % max
	GC float64 `json:"gc"`

	// Range that the primer spans on the
	Range ranged `json:"-"`
}

// newFrag creates a Frag from a match
func newFrag(m match, conf *config.Config) *Frag {
	fType := existing
	if m.circular {
		fType = circular
	}

	return &Frag{
		ID:       m.entry,
		uniqueID: m.uniqueID,
		Seq:      strings.ToUpper(m.seq),
		start:    m.queryStart,
		end:      m.queryEnd,
		db:       m.db,
		URL:      parseURL(m.entry),
		conf:     conf,
		fragType: fType,
	}
}

// parseURL turns a fragment identifier into a URL to its repository
func parseURL(id string) string {
	if strings.Contains(id, "addgene") {
		// create a link to the source Addgene page
		re := regexp.MustCompile("^.*addgene\\|(\\d*)")
		match := re.FindStringSubmatch(id)
		if len(match) > 0 {
			return fmt.Sprintf("https://www.addgene.org/%s/", match[1])
		}
	}

	if strings.Contains(id, "igem") {
		// create a source to the source iGEM page
		re := regexp.MustCompile("^.*igem\\|(\\w*)")
		match := re.FindStringSubmatch(id)
		if len(match) > 0 {
			return fmt.Sprintf("http://parts.igem.org/Part:%s", match[1])
		}
	}

	return ""
}

// newFlags is the plural of newFlag
func newFrags(matches []match, conf *config.Config) []*Frag {
	var frags []*Frag
	for _, m := range matches {
		frags = append(frags, newFrag(m, conf))
	}
	return frags
}

// copy returns a deep dopy of a Frag. used because nodes are mutated
// during assembly filling, and we don't want primers being shared between
// nodes in different assemblies
func (f *Frag) copy() (newFrag *Frag) {
	newFrag = &Frag{}
	copier.Copy(newFrag, f)

	return
}

// cost returns the estimated cost of a fragment. Combination of source and preparation
func (f *Frag) cost(procure bool) (c float64) {
	if procure {
		if strings.Contains(f.URL, "addgene") {
			c += f.conf.AddGeneVectorCost
		} else if strings.Contains(f.URL, "igem") {
			c += f.conf.IGEMPartCost
		}
	}

	if f.fragType == pcr && f.Primers != nil {
		c += float64(len(f.Primers[0].Seq)+len(f.Primers[1].Seq)) * f.conf.PCRBPCost
	} else if f.fragType == synthetic {
		c += f.conf.SynthFragmentCost(len(f.Seq))
	}

	return
}

// distTo returns the distance between the start of this Frag and the end of the other.
// assumes that this Frag starts before the other
// will return a negative number if this Frag overlaps with the other and positive otherwise
func (f *Frag) distTo(other *Frag) (bpDist int) {
	return other.start - f.end
}

// overlapsViaPCR returns whether this Frag could overlap the other Frag through homology
// created via PCR
func (f *Frag) overlapsViaPCR(other *Frag) bool {
	return f.distTo(other) <= f.conf.PCRMaxEmbedLength
}

// overlapsViaHomology returns whether this Frag already has sufficient overlap with the
// other Frag without any preparation like PCR
func (f *Frag) overlapsViaHomology(other *Frag) bool {
	return f.distTo(other) < -(f.conf.FragmentsMinHomology)
}

// synthDist returns the number of synthesized fragments that would need to be created
// between one Frag and another if the two were to be joined, with no existing
// fragments/nodes in-between, in an assembly
func (f *Frag) synthDist(other *Frag) (synthCount int) {
	dist := f.distTo(other)

	if f.overlapsViaPCR(other) {
		// if the dist is <MaxEmbedLength, we can PCR our way there
		// and add the mutated bp between the nodes with PCR
		return 0
	}

	floatDist := math.Max(1.0, float64(dist))

	// split up the distance between them by the max synthesized fragment size
	return int(math.Ceil(floatDist / float64(f.conf.SynthesisMaxLength)))
}

// costTo estimates the $ amount needed to get from this fragment
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

	if f.overlapsViaPCR(other) {
		if f.overlapsViaHomology(other) {
			// there's already enough overlap between this Frag and the one being tested
			// estimating two primers, primer length assumed to be 25bp
			return 50 * f.conf.PCRBPCost
		}

		// we have to create some additional primer sequence to reach the next fragment
		// estimating here that we'll add half of minHomology to both sides
		return float64(50+f.conf.FragmentsMinHomology) * f.conf.PCRBPCost
	}

	// we need to create a new synthetic fragment to get from this fragment to the next
	// to account for both the bps between them as well as the additional bps we need to add
	// for homology between the two
	fragLength := f.conf.FragmentsMinHomology + dist
	return f.conf.SynthFragmentCost(fragLength)
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

		// these nodes overlap by enough (via PCR or existing homology)
		if f.overlapsViaPCR(nodes[i]) {
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
	s1 := f.Seq
	if f.PCRSeq != "" {
		s1 = f.PCRSeq
	}

	s2 := other.Seq
	if other.PCRSeq != "" {
		s2 = other.PCRSeq
	}

	thisSeq := strings.ToUpper(s1)
	otherSeq := strings.ToUpper(s2)

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
			if j == len(otherSeq) || thisSeq[k] != otherSeq[j] {
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
// target is the vector's full sequence. We need it to build up the target
// vector's sequence
func (f *Frag) synthTo(next *Frag, target string) (synths []*Frag) {
	minHomology := f.conf.FragmentsMinHomology

	// check whether we need to make synthetic fragments to get
	// to the next fragment in the assembly
	fragCount := f.synthDist(next)
	if fragCount == 0 {
		return nil
	}

	// length of each synthesized fragment
	fragLength := f.distTo(next) / fragCount
	// account for homology on either end of each synthetic fragment
	fragLength += minHomology * 2
	if f.conf.SynthesisMinLength > fragLength {
		// need to synthesize at least Synthesis.MinLength bps
		fragLength = f.conf.SynthesisMinLength
	}

	// add to self to account for sequence across the zero-index (when sequence subselecting)
	targetLength := len(target)
	target = strings.ToUpper(target + target + target + target)

	// slide along the range of sequence to create synthetic fragments
	// and create one at each point, each w/ MinHomology for the fragment
	// before and after it
	synths = []*Frag{}
	start := f.end - minHomology // start w/ homology, move left
	for len(synths) < int(fragCount) {
		end := start + fragLength + 1
		seq := target[start+targetLength : end+targetLength]
		for hairpin(seq[len(seq)-minHomology:], f.conf) > f.conf.FragmentsMaxHairpinMelt {
			end += minHomology / 2
			seq = target[start+targetLength : end+targetLength]
		}

		synths = append(synths, &Frag{
			ID:       fmt.Sprintf("%s-synthetic-%d", f.ID, len(synths)+1),
			Seq:      seq,
			fragType: synthetic,
			conf:     f.conf,
		})

		start = end - minHomology
	}

	return
}

// setPrimers creates primers against a Frag and returns an error if:
//	1. the primers have an unacceptably high primer3 penalty score
//	2. the primers have off-targets in their source vector/fragment
func (f *Frag) setPrimers(last, next *Frag, seq string, conf *config.Config) (err error) {
	psExec := newPrimer3(last, f, next, seq, conf)

	// make input file and write to the fs
	// find how many bp of additional sequence need to be added
	// to the left and right primers (too large for primer3_core)
	addLeft, addRight, err := psExec.input(
		conf.FragmentsMinHomology,
		conf.FragmentsMaxHomology,
		conf.PCRMaxEmbedLength,
		conf.PCRMinLength,
	)
	if err != nil {
		return
	}

	if err = psExec.run(); err != nil {
		return
	}

	if err = psExec.parse(seq); err != nil {
		return
	}

	// update Frag's range, and add additional bp to the left and right primer if it wasn't included in the primer3 output
	mutatePrimers(f, seq, addLeft, addRight)

	// make sure the fragment's length is still long enough for PCR
	if f.end-f.start < conf.PCRMinLength {
		return fmt.Errorf(
			"failed to execute primer3: %s is %dbp, needs to be > %dbp",
			f.ID,
			f.end-f.start,
			conf.PCRMinLength,
		)
	}

	// 1. check for whether the primers have too have a pair penalty score
	if f.Primers[0].PairPenalty > conf.PCRP3MaxPenalty {
		errMessage := fmt.Sprintf(
			"primers have pair primer3 penalty score of %f, should be less than %f:\f%+v\f%+v",
			f.Primers[0].PairPenalty,
			conf.PCRP3MaxPenalty,
			f.Primers[0],
			f.Primers[1],
		)
		f.Primers = nil
		return fmt.Errorf(errMessage)
	}

	// 2. check for whether either of the primers have an off-target/mismatch
	var mismatchExists bool
	var mm match

	if f.fullSeq != "" {
		// we have the full sequence (it was included in the forward design)
		mismatchExists, mm, err = seqMismatch(f.Primers, f.ID, f.fullSeq, conf)
	} else {
		// otherwise, query the fragment from the DB (try to find it) and then check for mismatches
		mismatchExists, mm, err = parentMismatch(f.Primers, f.ID, f.db, conf)
	}

	if err != nil {
		f.Primers = nil
		return err
	}
	if mismatchExists {
		err = fmt.Errorf(
			"found a mismatching sequence %s for primers: %s, %s",
			mm.seq,
			f.Primers[0].Seq,
			f.Primers[1].Seq,
		)
		f.Primers = nil
		return err
	}

	os.Remove(psExec.in.Name()) // delete the temporary input and output files
	os.Remove(psExec.out.Name())

	return
}

// mutatePrimers adds additional bp to the sides of a Frag
// if there was additional homology bearing sequence that we were unable
// to add through primer3 alone
//
// it also updates the range, start + end, of the Frag to match that of the primers
//
// returning Frag for testing
func mutatePrimers(f *Frag, seq string, addLeft, addRight int) (mutated *Frag) {
	sl := len(seq)
	template := strings.ToUpper(seq + seq + seq)

	// update fragment sequence
	f.Seq = template[f.Primers[0].Range.start+len(seq) : f.Primers[1].Range.end+len(seq)+1]

	// change the Frag's start and end index to match those of the start and end index
	// of the primers, since the range may have shifted to get better primers
	f.start = f.Primers[0].Range.start
	f.end = f.Primers[1].Range.end

	// add bp to the left/FWD primer to match the fragment to the left
	if addLeft > 0 {
		oldStart := f.Primers[0].Range.start + sl
		f.Primers[0].Seq = template[oldStart-addLeft:oldStart] + f.Primers[0].Seq
		f.Primers[0].Range.start -= addLeft
	}

	// add bp to the right/REV primer to match the fragment to the right
	if addRight > 0 {
		oldEnd := f.Primers[1].Range.end + sl
		f.Primers[1].Seq = revComp(template[oldEnd+1:oldEnd+addRight+1]) + f.Primers[1].Seq
		f.Primers[1].Range.end += addRight
	}

	// update fragment sequence
	f.PCRSeq = template[f.Primers[0].Range.start+sl : f.Primers[1].Range.end+sl+1]

	return f
}

// ToString returns a string representation of a fragment's type
func (t fragType) String() string {
	return []string{"vector", "pcr", "synthetic", "existing"}[t]
}

// fragsCost returns the total cost of a slice of frags. Just the summation of their costs
func fragsCost(frags []*Frag) (cost float64) {
	for _, f := range frags {
		cost += f.cost(true)
	}

	return
}
