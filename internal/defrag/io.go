package defrag

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
)

// stderr is for logging to Stderr (without an annoying timestamp)
var stderr = log.New(os.Stderr, "", 0)

// Solution is a single solution to build up the target vector
type Solution struct {
	// Count is the number of fragments in this solution
	Count int `json:"count"`

	// Cost estimated from the primer and sequence lengths
	Cost float64 `json:"cost"`

	// Fragments used to build this solution
	Fragments []*Frag `json:"fragments"`
}

// Output is a struct containing design results for the assembly
type Output struct {
	// Target's name. In >example_CDS FASTA its "example_CDS"
	Target string `json:"target"`

	// Target's sequence
	TargetSeq string `json:"seq"`

	// Time, ex:
	// "2018-01-01 20:41:00"
	Time string `json:"time"`

	// Execution is the number of seconds it took to execute the command
	Execution float64 `json:"execution"`

	// VectorSynthesisCost is the cost of a full gene synthesis within a vector
	VectorSynthesisCost float64 `json:"vectorSynthesisCost"`

	// InsertSynthesisCost is the cost of just synthesizing the insert with homology for a linearized backbone
	InsertSynthesisCost float64 `json:"insertSynthesisCost"`

	// Solutions builds
	Solutions []Solution `json:"solutions"`
}

// Flags conatins parsed cobra Flags like "in", "out", "dbs", etc that
// are used by multiple commands
type Flags struct {
	// the name of the file to write the input to
	in string

	// the name of the file to write the output to
	out string

	// a list of dbs to run BLAST against (their names' on the filesystem)
	dbs []string

	// the backbone (optional) to insert the pieces into
	backbone Frag

	// slice of strings to weed out fragments from BLAST matches
	filters []string

	// percentage identity for finding building fragments in BLAST databases
	identity int
}

// inputParser contains methods for parsing flags from the input &cobra.Command
type inputParser struct{}

// NewFlags makes a new flags object manually. for testing
func NewFlags(in, out, backbone, enzymeName, filter string, dbs []string, addgene, igem bool) (*Flags, *config.Config) {
	c := config.New()

	if addgene {
		dbs = append(dbs, config.AddgeneDB)
	}
	if igem {
		dbs = append(dbs, config.IGEMDB)
	}

	p := inputParser{}
	parsedBB, err := p.parseBackbone(backbone, enzymeName, dbs, c)
	if err != nil {
		stderr.Fatal(err)
	}

	if strings.Contains(in, ",") {
		in = p.parseFeatureInput(strings.Fields(in))
	}

	return &Flags{
		in:       in,
		out:      out,
		dbs:      dbs,
		backbone: parsedBB,
		filters:  p.getFilters(filter),
		identity: 100,
	}, c
}

// parseCmdFlags gathers the in path, out path, etc from a cobra cmd object
// returns Flags and a Config struct for defrag.Vector or defrag.Fragments
func parseCmdFlags(cmd *cobra.Command, args []string) (*Flags, *config.Config) {
	var err error
	fs := &Flags{} // parsed flags
	p := inputParser{}
	c := config.New()

	if fs.in, err = cmd.Flags().GetString("in"); err != nil {
		if strings.ToLower(cmd.Name()) == "features" {
			fs.in = p.parseFeatureInput(args)
		} else if fs.in, err = p.guessInput(); err != nil {
			// check whether an input fail was specified
			stderr.Fatal(err)
		}
	}

	if fs.out, err = cmd.Flags().GetString("out"); err != nil {
		fs.out = p.guessOutput(fs.in) // guess at an output name
	}

	addgene, err := cmd.Flags().GetBool("addgene") // use addgene db?
	if err != nil {
		stderr.Fatalf("failed to parse addgene flag: %v", err)
	}

	igem, err := cmd.Flags().GetBool("igem") // use igem db?
	if err != nil {
		stderr.Fatalf("failed to parse igem flag: %v", err)
	}

	dbString, err := cmd.Flags().GetString("dbs")
	if err != nil && !addgene {
		stderr.Fatalf("failed to parse building fragments: %v", err)
	}

	filters, err := cmd.Flags().GetString("exclude")
	if err != nil {
		stderr.Fatalf("failed to parse filters: %v", err)
	}

	identity, err := cmd.Flags().GetInt("identity")
	if err != nil {
		identity = 100 // might be something other than `defrag vector`
	}

	if dbString == "" && !addgene && !igem {
		fmt.Println("no fragment databases chosen [-adi]: setting Addgene and iGEM to true")
		addgene = true
		igem = true
	}
	// read in the BLAST DB paths
	if fs.dbs, err = p.parseDBs(dbString, addgene, igem); err != nil || len(fs.dbs) == 0 {
		stderr.Fatalf("failed to find any fragment databases: %v", err)
	}

	// check if user asked for a specific backbone, confirm it exists in one of the dbs
	backbone, _ := cmd.Flags().GetString("backbone")

	// check if they also specified an enzyme
	enzymeName, _ := cmd.Flags().GetString("enzyme")

	// try to digest the backbone with the enzyme
	fs.backbone, err = p.parseBackbone(backbone, enzymeName, fs.dbs, c)
	if err != nil {
		stderr.Fatal(err)
	}

	// try to split the filter fields into a list
	fs.filters = p.getFilters(filters)

	// set identity for blastn searching
	fs.identity = identity

	return fs, c
}

// guessInput returns the first fasta file in the current directory. Is used
// if the user hasn't specified an input file
func (p *inputParser) guessInput() (in string, err error) {
	dir, _ := filepath.Abs(".")
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return
	}

	for _, file := range files {
		if file.IsDir() {
			continue
		}

		ext := strings.ToUpper(filepath.Ext(file.Name()))
		if ext == ".fa" || ext == ".fasta" {
			return file.Name(), nil
		}
	}

	return "", fmt.Errorf("failed: no input argument set and no fasta file found in %s", dir)
}

// parseFeatureInput turns the arguments to the features command into a
// CSV list of features in a single string
func (p *inputParser) parseFeatureInput(args []string) (out string) {
	commaSeparated := false
	for _, a := range args {
		if strings.Contains(a, ",") {
			commaSeparated = true
		}
	}

	// if it's a features command, concatenate the arguments in case they're feature names
	// with 'defrag features' the arguments are the feature names to use
	if commaSeparated {
		spacedSeparated := strings.Join(args, " ")
		splitByComma := strings.Split(spacedSeparated, ",")
		trimmedSeqs := []string{}
		for _, entry := range splitByComma {
			trimmedSeqs = append(trimmedSeqs, strings.TrimSpace(entry))
		}
		return strings.Join(trimmedSeqs, ",")
	}

	return strings.Join(args, " ")
}

// guessOutput gets an outpath path from an input path (if no output path is
// specified). It uses the same name as the input path to create an output
func (p *inputParser) guessOutput(in string) (out string) {
	ext := filepath.Ext(in)
	noExt := in[0 : len(in)-len(ext)]
	return noExt + ".defrag.json"
}

// parseDBs returns a list of absolute paths to BLAST databases
func (p *inputParser) parseDBs(dbs string, addgene, igem bool) (paths []string, err error) {
	if addgene {
		dbs += "," + config.AddgeneDB
	}
	if igem {
		dbs += "," + config.IGEMDB
	}

	if paths, err = p.dbPaths(dbs); err != nil {
		return nil, err
	}

	// make sure all the blast databases exist in the user's FS
	for _, path := range paths {
		if _, err := os.Stat(path); os.IsNotExist(err) {
			return nil, fmt.Errorf("failed to find a BLAST database at %s", path)
		}
	}

	return paths, nil
}

// dbPaths turns a single string of comma separated BLAST dbs into a
// slice of absolute paths to the BLAST dbs on the local fs
func (p *inputParser) dbPaths(dbList string) (paths []string, err error) {
	dbRegex := regexp.MustCompile(",")
	for _, db := range dbRegex.Split(dbList, -1) {
		dbTrimmed := strings.Trim(db, " ,")
		if dbTrimmed == "" {
			continue
		}

		absPath, err := filepath.Abs(dbTrimmed)
		if err != nil {
			return nil, fmt.Errorf("failed to create absolute path: %v", err)
		}
		paths = append(paths, absPath)
	}

	return
}

// parseBackbone takes a backbone, referenced by its id, and an enzyme to cleave the
// backbone, and returns the linearized backbone as a Frag
func (p *inputParser) parseBackbone(backbone, enzyme string, dbs []string, c *config.Config) (f Frag, err error) {
	// if no backbone was specified, return an empty Frag
	if backbone == "" {
		return Frag{}, nil
	}

	// try to digest the backbone with the enzyme
	if enzyme == "" {
		return Frag{}, fmt.Errorf("backbone passed, %s, without an enzyme to digest it", backbone)
	}

	// confirm that the backbone exists in one of the dbs (or local fs) gather it as a Frag if it does
	bbFrag, err := queryDatabases(backbone, dbs)
	if err != nil {
		return Frag{}, err
	}
	bbFrag.circular = true // assume it's circular here (hackish, fix)

	// gather the enzyme by name, err if it's unknown
	enz, err := p.getEnzyme(enzyme)
	if err != nil {
		return Frag{}, err
	}

	if f, err = digest(bbFrag, enz); err != nil {
		return Frag{}, err
	}

	return
}

// getEnzymes return the enzyme with the name passed. errors out if there is none
func (p *inputParser) getEnzyme(enzymeName string) (enzyme, error) {
	enzymeDB := NewEnzymeDB()
	if e, exists := enzymeDB.enzymes[enzymeName]; exists {
		return newEnzyme(e), nil
	}

	return enzyme{}, fmt.Errorf(
		`failed to find enzyme with name %s use "defrag enzymes" for a list of recognized enzymes`,
		enzymeName,
	)
}

// getFilters takes an input string and returns a list of strings to run against matches
// when filtering out possible building fragments
func (p *inputParser) getFilters(filterFlag string) []string {
	splitFunc := func(c rune) bool {
		return c == ' ' || c == ',' // space or comma separated
	}

	return strings.FieldsFunc(strings.ToUpper(filterFlag), splitFunc)
}

// read a FASTA file (by its path on local FS) to a slice of Fragments
func read(path string, feature bool) (fragments []Frag, err error) {
	if !filepath.IsAbs(path) {
		path, err = filepath.Abs(path)
		if err != nil {
			return nil, fmt.Errorf("failed to create path to input file: %s", err)
		}
	}

	dat, err := ioutil.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("failed to read input file: %s", err)
	}
	file := string(dat)

	path = strings.ToLower(path)
	if strings.HasSuffix(path, "fa") || strings.HasSuffix(path, "fasta") || file[0] == '>' {
		return readFasta(path, file)
	}

	fmt.Println(path)
	if strings.HasSuffix(path, "gb") || strings.HasSuffix(path, "gbk") || strings.HasSuffix(path, "genbank") {
		return readGenbank(path, file, feature)
	}

	return nil, fmt.Errorf("failed to parse %s: unrecognized file type", path)
}

// readFasta parses the multifasta file to fragments
func readFasta(path, contents string) (fragments []Frag, err error) {
	// split by newlines
	lines := strings.Split(contents, "\n")

	// read in the fragments
	var headerIndices []int
	var ids []string
	for i, line := range lines {
		if strings.HasPrefix(line, ">") {
			headerIndices = append(headerIndices, i)
			ids = append(ids, line[1:])
		}
	}

	// create a regex for cleaning the sequence
	var unwantedChars = regexp.MustCompile(`(?im)[^atgc]|\W`)

	// accumulate the sequences from between the headers
	var seqs []string
	for i, headerIndex := range headerIndices {
		nextLine := len(lines)
		if i < len(headerIndices)-1 {
			nextLine = headerIndices[i+1]
		}
		seqLines := lines[headerIndex+1 : nextLine]
		seqJoined := strings.Join(seqLines, "")
		seq := unwantedChars.ReplaceAllString(seqJoined, "")
		seqs = append(seqs, seq)
	}

	// build and return the new fragments
	for i, id := range ids {
		fragments = append(fragments, Frag{
			ID:  id,
			Seq: seqs[i],
		})
	}

	// opened and parsed file but found nothing
	if len(fragments) < 1 {
		return fragments, fmt.Errorf("failed to parse fragment(s) from %s", path)
	}

	return
}

// readGenbank parses a genbank file to fragments. Returns either fragments or parseFeatures,
// depending on the parseFeatures parameter
func readGenbank(path, contents string, parseFeatures bool) (fragments []Frag, err error) {
	genbankSplit := strings.Split(contents, "ORIGIN")

	if len(genbankSplit) != 2 {
		return nil, fmt.Errorf("failed to parse %s: improperly formatted genbank file", path)
	}

	seq := strings.ToUpper(genbankSplit[1])
	nonBpRegex := regexp.MustCompile("[^ATGC]")
	cleanedSeq := nonBpRegex.ReplaceAllString(seq, "")

	if parseFeatures {
		// parse each feature to a fragment (misnomer)
		splitOnFeatures := strings.Split(genbankSplit[0], "FEATURES")

		if len(splitOnFeatures) < 2 {
			return nil, fmt.Errorf("failed to parse features from %s", path)
		}

		featureSplitRegex := regexp.MustCompile("\\w+\\s+\\w+")
		featureStrings := featureSplitRegex.Split(splitOnFeatures[1], -1)

		features := []Frag{}
		for featureIndex, feature := range featureStrings {
			rangeRegex := regexp.MustCompile("(\\d*)\\.\\.(\\d*)")
			rangeIndexes := rangeRegex.FindStringSubmatch(feature)

			if len(rangeIndexes) < 3 {
				continue
			}

			start, err := strconv.Atoi(rangeIndexes[1])
			if err != nil {
				return nil, err
			}

			end, err := strconv.Atoi(rangeIndexes[2])
			if err != nil {
				return nil, err
			}
			featureSeq := cleanedSeq[start-1 : end] // make 0-indexed

			labelRegex := regexp.MustCompile("\\/label=(.*)")
			labelMatch := labelRegex.FindStringSubmatch(feature)
			label := ""
			if len(labelMatch) > 1 {
				label = labelMatch[1]
			} else {
				label = strconv.Itoa(featureIndex)
			}

			features = append(features, Frag{
				ID:  label,
				Seq: featureSeq,
			})
		}

		return features, nil
	}

	// parse just the file's sequence
	idRegex := regexp.MustCompile("LOCUS[ \\t]*([^ \\t]*)")
	id := idRegex.FindString(genbankSplit[0])

	if id == "" {
		return nil, fmt.Errorf("failed to parse locus from %s", path)
	}

	return []Frag{
		Frag{
			ID:  id,
			Seq: cleanedSeq,
		},
	}, nil
}

// write a slice of assemblies to the fs at the output path
//
// filename is the output file to write to
// target is the vector we tried to assemble
// assemblies are the solutions that can build up the target vector
func write(filename, targetName, targetSeq string, assemblies [][]*Frag, insertSeqLength int, conf *config.Config, seconds float64) (output []byte, err error) {
	// store save time, using same format as log.Println https://golang.org/pkg/log/#Println
	t := time.Now() // https://gobyexample.com/time-formatting-parsing
	time := fmt.Sprintf("%d/%02d/%02d %02d:%02d:%02d", t.Year(), t.Month(), t.Day(), t.Hour(), t.Minute(), t.Second())

	// calculate final cost of the assembly and fragment count
	solutions := []Solution{}
	for _, assembly := range assemblies {
		assemblyCost := 0.0
		assemblyFragmentIDs := make(map[string]bool)
		for _, f := range assembly {
			// freeze fragment type
			f.Type = f.fragType.String()

			// if it's already in the assembly, don't count cost twice
			if _, contained := assemblyFragmentIDs[f.ID]; contained {
				continue
			} else {
				assemblyFragmentIDs[f.ID] = true
			}

			// round to two decimal places
			f.Cost, err = roundCost(f.cost())
			if err != nil {
				return nil, err
			}

			// accumulate assembly cost
			assemblyCost += f.Cost
		}

		solutionCost, err := roundCost(assemblyCost)
		if err != nil {
			return nil, err
		}

		solutions = append(solutions, Solution{
			Count:     len(assembly),
			Cost:      solutionCost,
			Fragments: assembly,
		})
	}

	// sort solutions in increasing fragment count order
	sort.Slice(solutions, func(i, j int) bool {
		return solutions[i].Count < solutions[j].Count
	})

	// get the cost of full synthesis (comes in vector)
	fullSynthCost, err := roundCost(conf.SynthGeneCost(insertSeqLength))
	if err != nil {
		return nil, err
	}

	insertSynthCost, err := roundCost(conf.SynthFragmentCost(insertSeqLength))
	if err != nil {
		return nil, err
	}

	out := Output{
		Time:                time,
		Target:              targetName,
		TargetSeq:           strings.ToUpper(targetSeq),
		Execution:           seconds,
		Solutions:           solutions,
		VectorSynthesisCost: fullSynthCost,
		InsertSynthesisCost: insertSynthCost,
	}

	output, err = json.MarshalIndent(out, "", "  ")
	if err != nil {
		return output, fmt.Errorf("failed to serialize output: %v", err)
	}

	if err = ioutil.WriteFile(filename, output, 0666); err != nil {
		return output, fmt.Errorf("failed to write the output: %v", err)
	}
	return output, nil
}

// roundCost returns a float for cost to 2 decimal places
func roundCost(cost float64) (float64, error) {
	roundedString := fmt.Sprintf("%.2f", cost)

	return strconv.ParseFloat(roundedString, 64)
}
