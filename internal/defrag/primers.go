package defrag

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"os/exec"
	"path"
	"strconv"
	"strings"

	"github.com/jjtimmons/defrag/config"
)

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

// p3Exec is a utility struct for executing primer3 to create primers for a part
type p3Exec struct {
	// node that we're trying to create primers for
	n *node

	// the node before this one
	last *node

	// the node after this one
	next *node

	// the target sequence
	seq string

	// input file
	in string

	// output file
	out string

	// path to primer3 executable
	p3Path string

	// path to primer3 config folder (with trailing separator)
	p3Conf string

	// path to the primer3 io output
	p3Dir string
}

// setPrimers creates primers against a node and returns an error if:
//	1. the primers have an unacceptably high primer3 penalty score
//	2. the primers have off-targets in their source vector/fragment
func (n *node) setPrimers(last, next *node, seq string, conf *config.Config) (err error) {
	exec := newP3Exec(last, n, next, seq, conf)

	// make input file and write to the fs
	// find how many bp of additional sequence need to be added
	// to the left and right primers (too large for primer3_core)
	minHomology := conf.Fragments.MinHomology
	maxHomology := conf.Fragments.MaxHomology
	maxEmbedLength := conf.PCR.MaxEmbedLength
	addLeft, addRight, err := exec.input(minHomology, maxHomology, maxEmbedLength)
	if err != nil {
		return
	}

	if err = exec.run(); err != nil {
		return
	}

	if err = exec.parse(seq); err != nil {
		return
	}

	// update node's range, and add additional bp to the left and right primer if it wasn't included in the primer3 output
	mutateNodePrimers(n, seq, addLeft, addRight)

	// 1. check for whether the primers have too have a pair penalty score
	if n.primers[0].PairPenalty > conf.PCR.P3MaxPenalty {
		n.primers = nil
		return fmt.Errorf(
			"Primers have pair primer3 penalty score of %f, should be less than %f:\n%+v\n%+v",
			n.primers[0].PairPenalty,
			conf.PCR.P3MaxPenalty,
			n.primers[0],
			n.primers[1],
		)
	}

	// 2. check for whether either of the primers have an off-target/mismatch
	var mismatchExists bool
	var mm match

	if n.fullSeq != "" {
		mismatchExists, mm, err = seqMismatch(n.primers, n.id, n.fullSeq, conf)
	} else {
		// otherwise, query the fragment from the DB (try to find it) and then check for mismatches
		mismatchExists, mm, err = parentMismatch(n.primers, n.id, n.db, conf)
	}

	if err != nil {
		n.primers = nil
		return err
	}
	if mismatchExists {
		n.primers = nil
		return fmt.Errorf(
			"Found a mismatching sequence, %s, against the primer",
			mm.seq,
		)
	}
	return
}

// newP3Exec creates a p3Exec from a fragment
func newP3Exec(last, this, next *node, seq string, conf *config.Config) p3Exec {
	vendorConf := conf.Vendors()

	return p3Exec{
		n:      this,
		last:   last,
		next:   next,
		seq:    strings.ToUpper(seq),
		in:     path.Join(vendorConf.Primer3dir, this.id+".in"),
		out:    path.Join(vendorConf.Primer3dir, this.id+".out"),
		p3Path: vendorConf.Primer3core,
		p3Conf: vendorConf.Primer3config,
		p3Dir:  vendorConf.Primer3dir,
	}
}

// input makes a primer3 input settings file and writes it to the filesystem
//
// the primers on this node should account for creating homology
// against the last node and the next node if there isn't enough
// existing homology to begin with (the two nodes should share ~50/50)
//
// returning settings for unit testing only
func (p *p3Exec) input(minHomology, maxHomology, maxEmbedLength int) (bpAddLeft, bpAddRight int, err error) {
	// adjust the node's start and end index in the event that there's too much homology
	// with the neighboring fragment
	p.shrink(p.last, p.n, p.next, maxHomology) // could skip passing as a param, but this is a bit easier to test imo

	// calc the bps to add on the left and right side of this node
	addLeft := p.bpToAdd(p.last, p.n, minHomology)
	addRight := p.bpToAdd(p.n, p.next, minHomology)
	growPrimers := addLeft
	if growPrimers < addRight {
		growPrimers = addRight
	}

	start := p.n.start
	length := p.n.end - start

	// determine whether we need to add additional bp to the primers. if there's too much to
	// add, or if we're adding largely different amounts to the FWD and REV primer, we set
	// the number of bp to add as bpAddLeft and bpAddRight and do so after creating primers
	// that anneal to the template
	primerDiff := addLeft - addRight
	if primerDiff > 5 || primerDiff < 5 {
		// if one sides has a lot to add but the other doesn't, don't increase
		// the primer generation size in primer3. will instead concat the sequence on later
		// because we do not want to throw off the annealing temp for the primers
		bpAddLeft = addLeft
		bpAddRight = addRight
		growPrimers = 0
	} else if growPrimers > 36-26 {
		// we can't exceed 36 bp here (primer3 upper-limit), just create primers for the portion that
		// anneals to the template and add the other portion/seqs on later (in mutateNodePrimers)
		bpAddLeft = addLeft
		bpAddRight = addRight
		growPrimers = 0
	} else {
		start -= addLeft
		length = p.n.end - start
		length += addRight
	}

	// sizes to make the primers and target size (min, opt, and max)
	primerMin := 18 + growPrimers // defaults to 18
	primerOpt := 20 + growPrimers
	primerMax := 26 + growPrimers // defaults to 23

	// check whether we have wiggle room on the left or right hand sides to move the
	// primers inward (let primer3 pick better primers)
	//
	// also adjust start and length in case there's TOO large an overhang and we need
	// to trim it in one direction or the other
	leftBuffer := p.buffer(p.last.distTo(p.n), minHomology, maxEmbedLength)
	rightBuffer := p.buffer(p.n.distTo(p.next), minHomology, maxEmbedLength)

	// create the settings map from all instructions
	file := p.settings(
		p.seq,
		p.p3Conf,
		start,
		length,
		primerMin,
		primerOpt,
		primerMax,
		leftBuffer,
		rightBuffer,
	)

	if err := ioutil.WriteFile(p.in, file, 0666); err != nil {
		return 0, 0, fmt.Errorf("failed to write primer3 input file %v: ", err)
	}
	return
}

// shrink adjusts the start and end of a node in the scenario where
// it excessively overlaps a neighboring fragment. For example, if there's
// 700bp of overlap, this will trim it back so we just PCR a subselection of
// the node and keep the overlap beneath the upper limit
// TODO: consider giving minLength for PCR fragments precedent, this may subvert that rule
func (p *p3Exec) shrink(last, n, next *node, maxHomology int) *node {
	var shiftInLeft int
	var shiftInRight int

	if distLeft := last.distTo(n); distLeft < -maxHomology {
		// there's too much homology on the left side, we should move the node's start inward
		shiftInLeft = (-distLeft) - maxHomology
	}

	if distRight := n.distTo(next); distRight < -maxHomology {
		// there's too much homology on the right side, we should move the node's end inward
		shiftInRight = (-distRight) - maxHomology
	}

	// update the seq to slice for the new modified range
	n.start += shiftInLeft
	n.end -= shiftInRight

	if n.seq != "" {
		n.seq = n.seq[shiftInLeft : len(n.seq)-shiftInRight-shiftInLeft]
	}

	return n
}

// bpToAdd calculates the number of bp to add the end of a left node to create a junction
// with the rightmost node
func (p *p3Exec) bpToAdd(left, right *node, minHomology int) (bpToAdd int) {
	if synthDist := left.synthDist(right); synthDist == 0 {
		// we're not going to synth our way here, check that there's already enough homology
		if bpDist := left.distTo(right); bpDist > -minHomology {
			// this node will add half the homology to the last fragment
			// ex: 5 bp distance leads to 2.5bp + ~10bp additonal
			// ex: -10bp distance leads to ~0 bp additional:
			// 		other node is responsible for all of it
			bpToAdd = bpDist + (minHomology / 2)
		}
	}
	return
}

// buffer takes the dist from a one fragment to another and
// returns the length of the "buffer" in which the primers can be optimized (let primer3 pick)
//
// dist is positive if there's a gap between the start/end of a fragment and the start/end of
// the other and negative if they overlap
func (p *p3Exec) buffer(dist, minHomology, maxEmbedLength int) (buffer int) {
	if dist > maxEmbedLength {
		// we'll synthesize because the gap is so large, add 100bp of buffer
		return 100 // TODO: move "100" to settings
	}
	if dist < -minHomology {
		// there's enough additonal overlap that we can move this FWD primer inwards
		// but only enough to ensure that there's still minHomology bp overlap
		// and only enough so we leave the neighbor space for primer optimization too
		return (-dist - minHomology) / 2
	}
	return 0
}

// settingsMap returns a new settings map for the primer3 config files
// can either use pick_cloning_primers mode, if the start and end primers' locations
// are fixed, or pick_primer_list mode if we're letting the primers shift and allowing
// primer3 to pick the best ones. One side may be free to move and the other not
func (p *p3Exec) settings(
	seq, p3conf string,
	start, length, primerMin, primerOpt, primerMax int,
	leftBuffer, rightBuffer int) (file []byte) {
	// fmt.Println(start, length, leftBuffer, rightBuffer)

	// see primer3 manual or /vendor/primer3-2.4.0/settings_files/p3_th_settings.txt
	settings := map[string]string{
		"PRIMER_THERMODYNAMIC_PARAMETERS_PATH": p3conf,
		"PRIMER_NUM_RETURN":                    "1",
		"PRIMER_PICK_ANYWAY":                   "1",
		"SEQUENCE_TEMPLATE":                    seq + seq + seq,         // triple sequence
		"PRIMER_MIN_SIZE":                      strconv.Itoa(primerMin), // default 18
		"PRIMER_OPT_SIZE":                      strconv.Itoa(primerOpt),
		"PRIMER_MAX_SIZE":                      strconv.Itoa(primerMax),
		"PRIMER_EXPLAIN_FLAG":                  "1",
	}

	// if there is room to optimize, we let primer3 pick the best primers available in the space
	// http://primer3.sourceforge.net/primer3_manual.htm#SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
	if leftBuffer > 0 || rightBuffer > 0 {
		// V-ls  v-le               v-rs   v-re
		// ---------------------------------
		leftStart := start + len(seq)
		leftEnd := start + len(seq) + primerMax + leftBuffer
		rightStart := start + len(seq) + length - rightBuffer - primerMax
		excludeLength := rightStart - leftEnd

		settings["PRIMER_TASK"] = "generic"
		settings["PRIMER_PICK_LEFT_PRIMER"] = "1"
		settings["PRIMER_PICK_INTERNAL_OLIGO"] = "0"
		settings["PRIMER_PICK_RIGHT_PRIMER"] = "1"
		settings["PRIMER_MIN_TM"] = "54.0" // defaults to 57.0
		settings["PRIMER_MAX_TM"] = "67.0" // defaults to 63.0

		// ugly undoing of the above in case only one side has buffer
		if leftBuffer == 0 {
			settings["SEQUENCE_FORCE_LEFT_START"] = strconv.Itoa(leftStart)
		} else if rightBuffer == 0 {
			settings["SEQUENCE_FORCE_RIGHT_START"] = strconv.Itoa(rightStart)
		}

		settings["PRIMER_PRODUCT_SIZE_RANGE"] = fmt.Sprintf("%d-%d", excludeLength, length)
		settings["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = fmt.Sprintf("%d,%d,%d,%d", leftStart, leftBuffer+primerMax, rightStart, rightBuffer+primerMax)
	} else {
		// otherwise force the start and end of the PCR range
		settings["PRIMER_TASK"] = "pick_cloning_primers"
		settings["SEQUENCE_INCLUDED_REGION"] = fmt.Sprintf("%d,%d", start+len(seq), length)
	}

	var fileBuffer bytes.Buffer
	for key, val := range settings {
		fmt.Fprintf(&fileBuffer, "%s=%s\n", key, val)
	}
	fileBuffer.WriteString("=") // required at file's end

	return fileBuffer.Bytes()
}

// run the primer3 executable against the input file
func (p *p3Exec) run() (err error) {
	p3Cmd := exec.Command(
		p.p3Path,
		p.in,
		"-output", p.out,
		"-strict_tags",
	)

	// execute primer3 and wait on it to finish
	if output, err := p3Cmd.CombinedOutput(); err != nil {
		return fmt.Errorf("failed to execute primer3: %s: %v", string(output), err)
	}
	return
}

// parse the output into primers
//
// input is the target sequence we're building for. We need it to modulo
// the primer ranges
func (p *p3Exec) parse(input string) (err error) {
	file, err := ioutil.ReadFile(p.out)
	if err != nil {
		return
	}
	fileContents := string(file)

	// read in results into map, they're all 1:1
	results := make(map[string]string)
	for _, line := range strings.Split(fileContents, "\n") {
		keyVal := strings.Split(line, "=")
		if len(keyVal) > 1 {
			results[strings.TrimSpace(keyVal[0])] = strings.TrimSpace(keyVal[1])
		}
	}

	if p3Warnings := results["PRIMER_WARNING"]; p3Warnings != "" {
		return fmt.Errorf("Primer3 generated warnings: %s", p3Warnings)
	}

	if p3Error := results["PRIMER_ERROR"]; p3Error != "" {
		return fmt.Errorf("failed to execute primer3: %s", p3Error)
	}

	if numPairs := results["PRIMER_PAIR_NUM_RETURNED"]; numPairs == "0" {
		return fmt.Errorf("failed to create primers")
	}

	// read in a single primer from the output string file
	// side is either "LEFT" or "RIGHT"
	parsePrimer := func(side string) Primer {
		seq := results[fmt.Sprintf("PRIMER_%s_0_SEQUENCE", side)]
		tm := results[fmt.Sprintf("PRIMER_%s_0_TM", side)]
		gc := results[fmt.Sprintf("PRIMER_%s_0_GC_PERCENT", side)]
		penalty := results[fmt.Sprintf("PRIMER_%s_0_PENALTY", side)]
		pairPenalty := results["PRIMER_PAIR_0_PENALTY"]

		tmfloat, _ := strconv.ParseFloat(tm, 64)
		gcfloat, _ := strconv.ParseFloat(gc, 64)
		penaltyfloat, _ := strconv.ParseFloat(penalty, 64)
		pairfloat, _ := strconv.ParseFloat(pairPenalty, 64)

		primerRange := results[fmt.Sprintf("PRIMER_%s_0", side)]
		primerStart, _ := strconv.Atoi(strings.Split(primerRange, ",")[0])
		primerStart = primerStart - len(input)
		primerEnd := primerStart + len(seq)
		if side == "RIGHT" {
			primerStart -= len(seq)
			primerEnd = primerStart + len(seq)
		}

		return Primer{
			Seq:         seq,
			Strand:      side == "LEFT",
			Tm:          tmfloat,
			GC:          gcfloat,
			Penalty:     penaltyfloat,
			PairPenalty: pairfloat,
			Range: ranged{
				start: primerStart,
				end:   primerEnd,
			},
		}
	}

	p.n.primers = []Primer{
		parsePrimer("LEFT"),
		parsePrimer("RIGHT"),
	}
	return
}

// mutateNodePrimers adds additional bp to the sides of a node
// if there was additional homology bearing sequence that we were unable
// to add through primer3 alone
//
// it also updates the range, start + end, of the node to match that of the primers
//
// returning node for testing
func mutateNodePrimers(n *node, seq string, addLeft, addRight int) (mutated *node) {
	template := seq + seq

	// add bp to the left/FWD primer to match the fragment to the left
	if addLeft > 0 {
		oldStart := n.primers[0].Range.start
		n.primers[0].Seq = template[oldStart-addLeft:oldStart] + n.primers[0].Seq
		n.primers[0].Range.start -= addLeft
	}

	// add bp to the right/REV primer to match the fragment to the right
	if addRight > 0 {
		oldEnd := n.primers[1].Range.end
		n.primers[1].Seq = revComp(template[oldEnd+1:oldEnd+addRight+1]) + n.primers[1].Seq
		n.primers[1].Range.end += addRight
	}

	// change the node's start and end index to match those of the start and end index
	// of the primers, since the range may have shifted to get better primers
	n.start = n.primers[0].Range.start
	n.end = n.primers[1].Range.end

	// update the node's seq to reflect that change
	n.seq = template[n.start : n.end+1]

	return n
}

// revComp returns the reverse complement of a template sequence
func revComp(seq string) string {
	seq = strings.ToUpper(seq)

	revCompMap := map[rune]byte{
		'A': 'T',
		'T': 'A',
		'G': 'C',
		'C': 'G',
	}

	var revCompBuffer bytes.Buffer
	for _, c := range seq {
		revCompBuffer.WriteByte(revCompMap[c])
	}

	revCompBytes := revCompBuffer.Bytes()
	for i := 0; i < len(revCompBytes)/2; i++ {
		j := len(revCompBytes) - i - 1
		revCompBytes[i], revCompBytes[j] = revCompBytes[j], revCompBytes[i]
	}
	return string(revCompBytes)
}
