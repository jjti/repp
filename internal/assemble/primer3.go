package assemble

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/jjtimmons/decvec/config"
	"github.com/jjtimmons/decvec/internal/blast"
	"github.com/jjtimmons/decvec/internal/dvec"
)

// path to the primer3 executable and config folder
var (
	p3Path string
	p3Conf string
	p3Dir  string
)

// setPrimers creates primers on a PCR fragment and returns an error if
//	1. the primers have an unacceptably high primer3 penalty score
//	2. there are off-targets in the primers
func setPrimers(p *dvec.PCR) error {
	maxPairP := config.NewConfig().PCR.P3MaxPenalty

	handleP3 := func(err error) {
		// we should fail totally on any primer3 errors -- shouldn't happen
		if err != nil {
			log.Panic(err)
		}
	}

	exec := p3exec{
		f:   p,
		in:  path.Join(p3Dir, p.ID+".in"),
		out: path.Join(p3Dir, p.ID+".out"),
	}

	// make input file
	err := exec.input()
	handleP3(err)

	// execute
	err = exec.run()
	handleP3(err)

	// parse the results into primers for storing on the fragment
	primers, err := exec.parse()
	handleP3(err)
	p.Primers = primers

	// 1. check for whether the primers have too have a pair penalty score
	if p.Primers[0].PairPenalty > maxPairP {
		return fmt.Errorf(
			"primers have pair primer3 penalty score of %f, should be less than %f:\n%+v\n%+v",
			p.Primers[0].PairPenalty,
			maxPairP,
			p.Primers[0],
			p.Primers[1],
		)
	}

	// 2. check for whether either of the primers have an off-target/mismatch
	for _, primer := range p.Primers {
		mismatchExists, mismatch, err := blast.Mismatch(primer.Seq, p.Entry)
		handleP3(err) // shouldn't be erroring here either
		if mismatchExists {
			return fmt.Errorf(
				"found a mismatching sequence, %s, against the primer %s",
				mismatch.Seq,
				primer.Seq,
			)
		}
	}

	return nil
}

// p3Exec is a utility struct for executing primer3 to create primers for a part
type p3exec struct {
	// fragment that we're trying to create primers for
	f *dvec.PCR

	// input file
	in string

	// output file
	out string
}

// input makes the primer3 input settings file
func (p *p3exec) input() error {
	// create primer3 settings
	settings := map[string]string{
		"PRIMER_THERMODYNAMIC_PARAMETERS_PATH": p3Conf,
		"PRIMER_NUM_RETURN":                    "1",
		"PRIMER_TASK":                          "pick_cloning_primers",
		"PRIMER_PICK_ANYWAY":                   "1",
		"SEQUENCE_TEMPLATE":                    p.f.Seq,
		"SEQUENCE_INCLUDED_REGION":             fmt.Sprintf("0,%d", len(p.f.Seq)),
		"PRIMER_PRODUCT_SIZE_RANGE":            fmt.Sprintf("%d-%d", len(p.f.Seq), len(p.f.Seq)+1),
	}

	var fileContents string
	for key, val := range settings {
		fileContents += fmt.Sprintf("%s=%s\n", key, val)
	}
	fileContents += "=" // required

	// write to the fs
	inputFile, err := os.Create(p.in)
	if err != nil {
		return fmt.Errorf("failed to create primer3 input file %v: ", err)
	}
	defer inputFile.Close()
	_, err = inputFile.WriteString(fileContents)
	if err != nil {
		return fmt.Errorf("failed to create primer3 input file %v: ", err)
	}
	return nil
}

// run the primer3 executable on the input file
func (p *p3exec) run() error {
	p3Cmd := exec.Command(
		p3Path,
		p.in,
		"-output", p.out,
		"-strict_tags",
	)

	var stderr bytes.Buffer
	p3Cmd.Stderr = &stderr

	// execute primer3 and wait on it to finish
	err := p3Cmd.Run()
	if err != nil {
		fmt.Fprintln(os.Stderr, err, stderr.String())
		return err
	}
	return nil
}

// parse the output into primers for the part
func (p *p3exec) parse() ([]dvec.Primer, error) {
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

	// read in a single primer from the output string file
	// side is either "LEFT" or "RIGHT"
	parsePrimer := func(side string) dvec.Primer {
		seq := results[fmt.Sprintf("PRIMER_%s_0_SEQUENCE", side)]
		tm := results[fmt.Sprintf("PRIMER_%s_0_TM", side)]
		gc := results[fmt.Sprintf("PRIMER_%s_0_GC_PERCENT", side)]
		penalty := results[fmt.Sprintf("PRIMER_%s_0_PENALTY", side)]
		pairPenalty := results["PRIMER_PAIR_0_PENALTY"]

		tmfloat, _ := strconv.ParseFloat(tm, 32)
		gcfloat, _ := strconv.ParseFloat(gc, 32)
		penaltyfloat, _ := strconv.ParseFloat(penalty, 32)
		pairfloat, _ := strconv.ParseFloat(pairPenalty, 32)

		return dvec.Primer{
			Seq:         seq,
			Strand:      side == "LEFT",
			Tm:          float32(tmfloat),
			GC:          float32(gcfloat),
			Penalty:     float32(penaltyfloat),
			PairPenalty: float32(pairfloat),
		}
	}

	return []dvec.Primer{
		parsePrimer("LEFT"),
		parsePrimer("RIGHT"),
	}, nil
}

// create the primer3 path, error Out if we can't find the executable or the config folder
func init() {
	// make sure the primer3 binary and settings folder are defined
	p3Path = filepath.Join("..", "..", "vendor", "primer3-2.4.0", "src", "primer3_core")
	// TODO: fix this forward slash at the end using an OS-specific path separator
	p3Conf = filepath.Join("..", "..", "vendor", "primer3-2.4.0", "src", "primer3_config") + "/"

	_, err := os.Stat(p3Path)
	if err != nil {
		log.Fatalf("failed to locate primer3 executable: %v", err)
	}

	_, err = os.Stat(p3Conf)
	if err != nil {
		log.Fatalf("failed to locate primer3 config folder: %v", err)
	}

	// make a folder for primer3 io
	p3Dir = filepath.Join("..", "..", "bin", "primer3")
	err = os.MkdirAll(p3Dir, os.ModePerm)
	if err != nil {
		log.Fatalf("failed to create a primer3 outut dir: %v", err)
	}
}
