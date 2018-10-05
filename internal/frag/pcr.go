package frag

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
)

// path to the primer3 executable and config folder
var (
	p3Path string
	p3Conf string
	p3Dir  string
)

// primer is a single primer used to ampligy a parent fragment
type primer struct {
	// seq of the primer (In 5' to 3' direction)
	seq string

	// strand of the primer; true if template, false if complement
	strand bool

	// Start of the fragment (0-index)
	Start int

	// End of the fragment (0-indexed)
	End int

	// Penalty score
	Penalty float32

	// PairPenalty score from primer3
	PairPenalty float32

	// Tm of the primer
	Tm float32

	// GC % max
	GC float32
}

// PCR is a PCR building fragment. It includes primers to create
// the fragment (using its template sequence/vector for addgene as a source)
// and the sequence of the resulting fragment
type PCR struct {
	// ID of this fragment
	ID string

	// Seq of this fragment
	Seq string

	// Entry of this fragment In the DB that it came from
	// Used to look for off-targets
	Entry string

	// Start of the fragment (0-indexed)
	Start int

	// End of the fragment (0-indexed)
	End int

	// Primers necessary to create this PCR Fragment from the template sequence
	Primers []primer
}

// p3Exec is a utility struct for executing primer3 to create primers for a part
type p3exec struct {
	// fragment that we're trying to create primers for
	Frag *PCR

	// input file
	In string

	// output file
	Out string
}

// input is for making the primer3 input settings file
func (p *p3exec) input() error {
	// create primer3 settings
	settings := map[string]string{
		"PRIMER_THERMODYNAMIC_PARAMETERS_PATH": p3Conf,
		"PRIMER_NUM_RETURN":                    "1",
		"PRIMER_TASK":                          "pick_cloning_primers",
		"PRIMER_PICK_ANYWAY":                   "1",
		"SEQUENCE_TEMPLATE":                    p.Frag.Seq,
		"SEQUENCE_INCLUDED_REGION":             fmt.Sprintf("0,%d", len(p.Frag.Seq)),
		"PRIMER_PRODUCT_SIZE_RANGE":            fmt.Sprintf("%d-%d", len(p.Frag.Seq), len(p.Frag.Seq)+1),
	}

	var fileContents string
	for key, val := range settings {
		fileContents += fmt.Sprintf("%s=%s\n", key, val)
	}
	fileContents += "=" // required

	// write to the fs
	inputFile, err := os.Create(p.In)
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
		p.In,
		"-output", p.Out,
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
func (p *p3exec) parse() ([]primer, error) {
	file, err := ioutil.ReadFile(p.Out)
	if err != nil {
		return nil, err
	}
	fileS := string(file)

	// read in results into map, they're all 1:1
	results := map[string]string{}
	for _, line := range strings.Split(fileS, "\n") {
		keyVal := strings.Split(line, "=")
		if len(keyVal) == 2 {
			results[keyVal[0]] = results[keyVal[1]]
		}
	}

	// read in a single primer from the output string file
	// side is either "LEFT" or "RIGHT"
	parsePrimer := func(side string) primer {
		seq := results[fmt.Sprintf("PRIMER_%s_0_SEQUENCE", side)]

		tm := results[fmt.Sprintf("PRIMER_%s_0_TM", side)]
		penalty := results[fmt.Sprintf("PRIMER_%s_0_PENALTY", side)]
		pairPenalty := results["PRIMER_PAIR_0_PENALTY"]

		tmfloat, _ := strconv.ParseFloat(tm, 32)
		penaltyfloat, _ := strconv.ParseFloat(penalty, 32)
		pairfloat, _ := strconv.ParseFloat(pairPenalty, 32)

		fmt.Println(tm, penalty)

		return primer{
			seq:         seq,
			strand:      side == "LEFT",
			Tm:          float32(tmfloat),
			Penalty:     float32(penaltyfloat),
			PairPenalty: float32(pairfloat),
		}
	}

	return []primer{
		parsePrimer("LEFT"),
		parsePrimer("RIGHT"),
	}, nil
}

// SetPrimers creates primers on a PCR fragment and returns an error if
//	1. there are off-targets In the primers
//	2. the primers have an unacceptably high primer3 penalty score
func (p *PCR) SetPrimers() error {
	handle := func(err error) {
		// we should fail totally on any primer3 errors -- shouldn't happen
		if err != nil {
			log.Panic(err)
		}
	}

	exec := p3exec{
		Frag: p,
		In:   path.Join(p3Dir, p.ID+".in"),
		Out:  path.Join(p3Dir, p.ID+".out"),
	}

	// make input file
	err := exec.input()
	handle(err)

	// execute
	err = exec.run()
	handle(err)

	// parse the results into primers for storing on the fragment
	primers, err := exec.parse()
	handle(err)
	p.Primers = primers

	return nil
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
