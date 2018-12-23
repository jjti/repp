package defrag

import (
	"fmt"
	"io/ioutil"
	"os/exec"
	"path"
	"strconv"
	"strings"

	"github.com/jjtimmons/defrag/config"
)

// p3Exec is a utility struct for executing primer3 to create primers for a part
type p3Exec struct {
	// node that we're trying to create primers for
	n *node

	// the node before this one
	last node

	// the node after this one
	next node

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

// newP3Exec creates a p3Exec from a fragment
func newP3Exec(last, this, next node, seq string, conf *config.Config) p3Exec {
	vendorConf := conf.Vendors()

	return p3Exec{
		n:      &this,
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

// primers creates primers against a node and return an error if
//	1. the primers have an unacceptably high primer3 penalty score
//	2. the primers have off-targets in their parent source
func primers(last, this, next node, seq string, conf *config.Config) (primers []Primer, err error) {
	exec := newP3Exec(last, this, next, seq, conf)
	vendorConfig := conf.Vendors()

	// make input file, figure out how to create primers that share homology
	// with neighboring nodes
	if err = exec.input(conf.Fragments.MinHomology); err != nil {
		return
	}

	if err = exec.run(); err != nil {
		return
	}

	if primers, err = exec.parse(); err != nil {
		return
	}

	// 1. check for whether the primers have too have a pair penalty score
	if primers[0].PairPenalty > conf.PCR.P3MaxPenalty {
		return nil, fmt.Errorf(
			"Primers have pair primer3 penalty score of %f, should be less than %f:\n%+v\n%+v",
			primers[0].PairPenalty,
			conf.PCR.P3MaxPenalty,
			primers[0],
			primers[1],
		)
	}

	// 2. check for whether either of the primers have an off-target/mismatch
	for _, primer := range primers {
		// the node's id is the same as the entry ID in the database
		mismatchExists, mismatch, err := Mismatch(primer.Seq, this.id, conf.DBs, vendorConfig)

		if err != nil {
			return nil, err
		}

		if mismatchExists {
			return nil, fmt.Errorf(
				"Found a mismatching sequence, %s, against the primer %s",
				mismatch.Seq,
				primer.Seq,
			)
		}
	}

	return
}

// input makes the primer3 input settings file and writes it to the filesystem
//
// the primers on this node should account for creating homology
// against the last node and the next node if there isn't enough
// existing homology to begin with (the two nodes should share ~50/50)
func (p *p3Exec) input(minHomology int) error {
	// calc the # of bp this node shares with another
	bpToShare := func(left, right node) (bpToAdd int) {
		// calc the # of bp the left node is responsible with the right one
		bpToAdd = 0
		if synthDist := left.synthDist(right); synthDist == 0 {
			// we're not going to synth our way here, check that there's already enough homology
			if bpDist := left.distTo(right); bpDist > -(minHomology) {
				// this node will add half the homology to the last fragment
				// ex: 5 bp distance leads to 2.5bp + ~10bp additonal
				// ex: -10bp distance leads to ~0 bp additional:
				//	other node is responsible for all of it
				bpToAdd = bpDist + (minHomology / 2)
			}
		}
		return
	}

	// calc the bps to add on the left and right side of this node
	addLeft := bpToShare(p.last, *p.n)
	addRight := bpToShare(*p.n, p.next)
	maxAdded := addLeft
	if addRight > maxAdded {
		maxAdded = addRight
	}

	// the node's range plus the additional bp added because of adding homology
	start := p.n.start - addLeft
	length := p.n.end - start + 1
	length += addRight

	// sizes to make the primers and target size (min, opt, and max)
	primerMin := 18 // defaults to 18
	primerOpt := 20
	primerMax := 25 // defaults to 23
	if maxAdded > 0 {
		maxAdded += 2

		// we can't exceed 36 here
		if maxAdded > 36-primerMax {
			maxAdded = 36 - primerMax
		}

		primerMin += maxAdded
		primerOpt += maxAdded
		primerMax += maxAdded
	}

	// see primer3 manual or /vendor/primer3-2.4.0/settings_files/p3_th_settings.txt
	settings := map[string]string{
		"PRIMER_THERMODYNAMIC_PARAMETERS_PATH": p.p3Conf,
		"PRIMER_NUM_RETURN":                    "1",
		"PRIMER_TASK":                          "pick_cloning_primers",
		"PRIMER_PICK_ANYWAY":                   "1",
		"SEQUENCE_TEMPLATE":                    p.seq + p.seq + p.seq, // triple sequence
		"SEQUENCE_INCLUDED_REGION":             fmt.Sprintf("%d,%d", start+len(p.seq), length),
		"PRIMER_MIN_SIZE":                      strconv.Itoa(primerMin), // default 18
		"PRIMER_OPT_SIZE":                      strconv.Itoa(primerOpt),
		"PRIMER_MAX_SIZE":                      strconv.Itoa(primerMax),
		"PRIMER_EXPLAIN_FLAG":                  "1",
	}

	// check whether we have wiggle room on the left or right hand sides to move the
	// primers inward (let primer3 pick better primers)
	lastDist := p.last.distTo(*p.n)
	if 0 > lastDist {
		lastDist = 0
	}

	nextDist := p.n.distTo(p.next)
	if 0 > nextDist {
		nextDist = 0
	}

	// change 5 to the number of bp we're willing to "insert" via PCR
	// checking here for whether we're going to synthetisize the neighboring fragments
	//
	// if we are, there's more room to pick better primers at the ends of the fragment
	// (or the best combination of primers within the available space)
	//
	// http://primer3.sourceforge.net/primer3_manual.htm#SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
	if len(p.seq) > 300 && (lastDist > 5 || nextDist > 5) {
		buffer := 150

		leftStart := start + len(p.seq)
		leftEnd := start + len(p.seq) + buffer
		rightStart := start + len(p.seq) + length - buffer
		excludeLength := rightStart - leftEnd

		if excludeLength > 0 {
			delete(settings, "PRIMER_PICK_ANYWAY")
			settings["PRIMER_TASK"] = "pick_primer_list"
			settings["PRIMER_PICK_LEFT_PRIMER"] = "1"
			settings["PRIMER_PICK_INTERNAL_OLIGO"] = "0"
			settings["PRIMER_PICK_RIGHT_PRIMER"] = "1"

			// ugly undoing of the above in case only one side has buffer
			if lastDist <= 5 {
				settings["SEQUENCE_FORCE_LEFT_START"] = strconv.Itoa(start + len(p.seq))
			} else if nextDist <= 5 {
				settings["SEQUENCE_FORCE_RIGHT_START"] = strconv.Itoa(start + len(p.seq) + length)
			}

			settings["PRIMER_PRODUCT_SIZE_RANGE"] = fmt.Sprintf("%d-%d", excludeLength, length)
			settings["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = fmt.Sprintf("%d,%d,%d,%d", leftStart, buffer, rightStart, buffer)
		}
	}

	var fileContents string
	for key, val := range settings {
		fileContents += fmt.Sprintf("%s=%s\n", key, val)
	}
	fileContents += "=" // required at file's end

	if err := ioutil.WriteFile(p.in, []byte(fileContents), 0666); err != nil {
		return fmt.Errorf("Failed to create primer3 input file %v: ", err)
	}
	return nil
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
		return fmt.Errorf("Failed to execute primer3: %s: %v", string(output), err)
	}
	return nil
}

// parse the output into primers
func (p *p3Exec) parse() (primers []Primer, err error) {
	fmt.Println(p.out)
	file, err := ioutil.ReadFile(p.out)
	if err != nil {
		return nil, err
	}
	fileS := string(file)

	// read in results into map, they're all 1:1
	results := make(map[string]string)
	for _, line := range strings.Split(fileS, "\n") {
		keyVal := strings.Split(line, "=")
		if len(keyVal) > 1 {
			results[strings.TrimSpace(keyVal[0])] = strings.TrimSpace(keyVal[1])
		}
	}

	if p3Warnings := results["PRIMER_WARNING"]; p3Warnings != "" {
		return nil, fmt.Errorf("Primer3 generated warnings: %s", p3Warnings)
	}

	if p3Error := results["PRIMER_ERROR"]; p3Error != "" {
		return nil, fmt.Errorf("Failed to execute primer3: %s", p3Error)
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

		return Primer{
			Seq:         seq,
			Strand:      side == "LEFT",
			Tm:          tmfloat,
			GC:          gcfloat,
			Penalty:     penaltyfloat,
			PairPenalty: pairfloat,
		}
	}

	return []Primer{
		parsePrimer("LEFT"),
		parsePrimer("RIGHT"),
	}, nil
}
