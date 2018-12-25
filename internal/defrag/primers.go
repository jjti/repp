package defrag

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"math"
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

// setPrimers creates primers against a node and return an error if
//	1. the primers have an unacceptably high primer3 penalty score
//	2. the primers have off-targets in their parent source
func (n *node) setPrimers(last, next *node, seq string, conf *config.Config) (err error) {
	exec := newP3Exec(last, n, next, seq, conf)

	// make input file, figure out how to create primers that share homology
	// with neighboring nodes
	addLeft, addRight, err := exec.input(conf.Fragments.MinHomology, conf.PCR.MaxEmbedLength)
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
	for _, primer := range n.primers {
		// the node's id is the same as the entry ID in the database
		mismatchExists, mm, err := mismatch(primer.Seq, n.id, n.db, conf)

		if err != nil {
			n.primers = nil
			return err
		}

		if mismatchExists {
			n.primers = nil
			return fmt.Errorf(
				"Found a mismatching sequence, %s, against the primer %s",
				mm.seq,
				primer.Seq,
			)
		}
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
func (p *p3Exec) input(minHomology, maxEmbedLength int) (bpAddLeft, bpAddRight int, err error) {
	// calc the bps to add on the left and right side of this node
	addLeft := bpToShare(p.last, p.n, minHomology)
	addRight := bpToShare(p.n, p.next, minHomology)
	growPrimers := addLeft
	if growPrimers < addRight {
		growPrimers = addRight
	}

	start := p.n.start
	length := p.n.end - start

	// if one sides has a lot to add but the other doesn't, don't increase
	// the primer generation size in primer3. will instead concat the sequence on later
	// because we do not want to throw off the annealing temp for the primers
	if math.Abs(float64(addLeft)-float64(addRight)) > 6.0 {
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
	distFromLast := p.last.distTo(*p.n)
	distFromNext := p.n.distTo(*p.next)

	leftBuffer := 0
	if distFromLast > maxEmbedLength {
		// we'll synthesize on left, add 100bp of buffer
		// TODO: move to settings
		leftBuffer = 100
	} else if distFromLast < -minHomology {
		// there's enough additonal overlap that we can move this FWD primer inwards
		// but only enough to ensure that there's still minHomology bp overlap
		// and only enough so we leave the neighbor space for primer optimization too
		leftBuffer = (-distFromLast - minHomology) / 2
	}

	rightBuffer := 0
	if distFromNext > maxEmbedLength {
		// will we have to synthesize on the right anyway, add 100bp of buffer
		rightBuffer = 100
	} else if distFromNext < -minHomology {
		// we can shift the REV primer inwards by enough to ensure there's still
		// minHomology bp overlap
		// and leave enough space for the neighboring fragment to optimize too
		rightBuffer = (-distFromNext - minHomology) / 2
	}

	// create the settings map from all instructions
	file := settingsFile(
		p.seq,
		p.p3Conf,
		start,
		length,
		primerMin,
		primerOpt,
		primerMax,
		leftBuffer,
		rightBuffer)

	if err := ioutil.WriteFile(p.in, file, 0666); err != nil {
		return 0, 0, fmt.Errorf("failed to write primer3 input file %v: ", err)
	}
	return
}

// bpToShare calculates the number of bp to add the end of a node to creat a junction between them
// returns the number of bp that needs to be added to each node
func bpToShare(left, right *node, minHomology int) (bpToAdd int) {
	if synthDist := left.synthDist(*right); synthDist == 0 {
		// we're not going to synth our way here, check that there's already enough homology
		if bpDist := left.distTo(*right); bpDist > -minHomology {
			// this node will add half the homology to the last fragment
			// ex: 5 bp distance leads to 2.5bp + ~10bp additonal
			// ex: -10bp distance leads to ~0 bp additional:
			// 		other node is responsible for all of it
			bpToAdd = bpDist + (minHomology / 2)
		}
	}
	return
}

// settingsMap returns a new settings map for the primer3 config files
// can either use pick_cloning_primers mode, if the start and end primers' locations
// are fixed, or pick_primer_list mode if we're letting the primers shift and allowing
// primer3 to pick the best ones. One side may be free to move and the other not
func settingsFile(
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
func (p *p3Exec) run() error {
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
	return nil
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

	return nil
}

// mutateNodePrimers adds additional bp to the sides of a node
// if there was additional homology bearing sequence that we were unable
// to add through primer3 alone
//
// it also updates the range, start + end, of the node to match that of the primers
//
// returning node for testing
func mutateNodePrimers(n *node, seq string, addLeft, addRight int) (mutated *node) {
	// fmt.Println(n.id, addLeft, addRight, n.primers[0].Range.start, n.primers[1].Range.end)

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
