package repp

import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/jjtimmons/repp/config"
	"github.com/spf13/cobra"
)

var (
	// stderr is for logging to Stderr (without an annoying timestamp)
	stderr = log.New(os.Stderr, "", 0)
)

// Flags contains parsed cobra Flags like "in", "out", "dbs", etc that are used by multiple commands.
type Flags struct {
	// the name of the file to write the input from
	in string

	// the name of the file to write the output to
	out string

	// a list of dbs to run BLAST against (their names' on the filesystem)
	dbs []string

	// the backbone (optional) to insert the pieces into
	backbone *Frag

	// backbone meta (name, enzyme used to cut it, cut index)
	backboneMeta *Backbone

	// slice of strings to weed out fragments from BLAST matches
	filters []string

	// percentage identity for finding building fragments in BLAST databases
	identity int
}

// inputParser contains methods for parsing flags from the input &cobra.Command.
type inputParser struct{}

// NewFlags makes a new flags object manually. for testing.
func NewFlags(
	in, out, backbone, enzymeName, filter string,
	dbs []string, addgene, igem, dnasu bool) (*Flags, *config.Config) {
	c := config.New()

	if addgene {
		dbs = append(dbs, config.AddgeneDB)
	}
	if igem {
		dbs = append(dbs, config.IGEMDB)
	}
	if dnasu {
		dbs = append(dbs, config.DNASUDB)
	}

	p := inputParser{}
	parsedBB, bbMeta, err := p.parseBackbone(backbone, enzymeName, dbs, c)
	if err != nil {
		stderr.Fatal(err)
	}

	if strings.Contains(in, ",") {
		in = p.parseFeatureInput(strings.Fields(in))
	}

	return &Flags{
		in:           in,
		out:          out,
		dbs:          dbs,
		backbone:     parsedBB,
		backboneMeta: bbMeta,
		filters:      p.getFilters(filter),
		identity:     98,
	}, c
}

// parseCmdFlags gathers the in path, out path, etc from a cobra cmd object
// returns Flags and a Config struct for repp.Vector or repp.Fragments.
func parseCmdFlags(cmd *cobra.Command, args []string, strict bool) (*Flags, *config.Config) {
	var err error
	fs := &Flags{} // parsed flags
	p := inputParser{}
	c := config.New()

	if fs.in, err = cmd.Flags().GetString("in"); fs.in == "" || err != nil {
		cmdName := strings.ToLower(cmd.Name())
		if cmdName == "features" {
			fs.in = p.parseFeatureInput(args)
		} else if cmdName == "sequence" && len(args) > 0 {
			fs.in = "input.fa"
			if err = ioutil.WriteFile(fs.in, []byte(fmt.Sprintf(">target_sequence\n%s", args[0])), 0644); err != nil {
				stderr.Fatal(err)
			}
		} else if fs.in, err = p.guessInput(); strict && err != nil {
			// check whether an input fail was specified
			cmd.Help()
			stderr.Fatal(err)
		}
	}

	if fs.out, err = cmd.Flags().GetString("out"); strict && (fs.out == "" || err != nil) {
		fs.out = p.guessOutput(fs.in) // guess at an output name

		if fs.out == "" {
			cmd.Help()
			stderr.Fatal("no output path")
		}
	}

	addgene, err := cmd.Flags().GetBool("addgene") // use addgene db?
	if strict && err != nil {
		cmd.Help()
		stderr.Fatalf("failed to parse addgene flag: %v", err)
	}

	igem, err := cmd.Flags().GetBool("igem") // use igem db?
	if strict && err != nil {
		cmd.Help()
		stderr.Fatalf("failed to parse igem flag: %v", err)
	}

	dnasu, err := cmd.Flags().GetBool("dnasu") // use dnasu db?
	if strict && err != nil {
		cmd.Help()
		stderr.Fatalf("failed to parse dnasu flag: %v", err)
	}

	dbString, err := cmd.Flags().GetString("dbs")
	if strict && err != nil && !addgene {
		cmd.Help()
		stderr.Fatalf("failed to parse building fragments: %v", err)
	}

	filters, err := cmd.Flags().GetString("exclude")
	if strict && err != nil {
		cmd.Help()
		stderr.Fatalf("failed to parse filters: %v", err)
	}
	// try to split the filter fields into a list
	fs.filters = p.getFilters(filters)

	identity, err := cmd.Flags().GetInt("identity")
	if err != nil {
		identity = 100 // might be something other than `repp vector`
	}
	// set identity for blastn searching
	fs.identity = identity

	if dbString == "" && !addgene && !igem && !dnasu {
		fmt.Println("no fragment databases chosen [-agu]: using Addgene, DNASU, and iGEM by default")
		addgene = true
		igem = true
		dnasu = true
	}
	// read in the BLAST DB paths
	if fs.dbs, err = p.parseDBs(dbString, addgene, igem, dnasu); err != nil || len(fs.dbs) == 0 {
		stderr.Fatalf("failed to find any fragment databases: %v", err)
	}

	// check if user asked for a specific backbone, confirm it exists in one of the dbs
	backbone, _ := cmd.Flags().GetString("backbone")

	// check if they also specified an enzyme
	enzymeName, _ := cmd.Flags().GetString("enzyme")

	// try to digest the backbone with the enzyme
	fs.backbone, fs.backboneMeta, err = p.parseBackbone(backbone, enzymeName, fs.dbs, c)
	if strict && err != nil {
		stderr.Fatal(err)
	}

	return fs, c
}

// guessInput returns the first fasta file in the current directory. Is used
// if the user hasn't specified an input file.
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
// CSV list of features in a single string.
func (p *inputParser) parseFeatureInput(args []string) (out string) {
	commaSeparated := false
	for _, a := range args {
		if strings.Contains(a, ",") {
			commaSeparated = true
		}
	}

	// if it's a features command, concatenate the arguments in case they're feature names
	// with 'repp features' the arguments are the feature names to use
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
// specified). It uses the same name as the input path to create an output.
func (p *inputParser) guessOutput(in string) (out string) {
	ext := filepath.Ext(in)
	noExt := in[0 : len(in)-len(ext)]
	return noExt + ".output.json"
}

// parseDBs returns a list of absolute paths to BLAST databases.
func (p *inputParser) parseDBs(dbs string, addgene, igem, dnasu bool) (paths []string, err error) {
	if addgene {
		dbs += "," + config.AddgeneDB
	}
	if igem {
		dbs += "," + config.IGEMDB
	}
	if dnasu {
		dbs += "," + config.DNASUDB
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
// slice of absolute paths to the BLAST dbs on the local fs.
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
// backbone, and returns the linearized backbone as a Frag.
func (p *inputParser) parseBackbone(
	bbName, enzyme string,
	dbs []string,
	c *config.Config,
) (f *Frag, backbone *Backbone, err error) {
	// if no backbone was specified, return an empty Frag
	if bbName == "" {
		return &Frag{}, &Backbone{}, nil
	}

	// confirm that the backbone exists in one of the dbs (or local fs) gather it as a Frag if it does
	bbFrag, err := queryDatabases(bbName, dbs)
	if err != nil {
		return &Frag{}, &Backbone{}, err
	}

	// check if it's an iGEM backbone (don't digest)
	if igemBackbone(bbFrag.ID) {
		return &Frag{
				ID:       bbFrag.ID,
				uniqueID: "backbone",
				Seq:      bbFrag.Seq[:len(bbFrag.Seq)/2],
				fragType: linear,
				db:       bbFrag.db,
			},
			&Backbone{
				URL:              parseURL(bbFrag.ID, bbFrag.db),
				Seq:              bbFrag.Seq,
				Enzyme:           "",
				RecognitionIndex: 0,
				Forward:          true,
			},
			nil
	}

	// try to digest the backbone with the enzyme
	if enzyme == "" {
		return &Frag{},
			&Backbone{},
			fmt.Errorf("backbone passed, %s, without an enzyme to digest it", bbName)
	}

	// gather the enzyme by name, err if it's unknown
	enz, err := p.getEnzyme(enzyme)
	enz.name = enzyme
	if err != nil {
		return &Frag{}, &Backbone{}, err
	}

	if f, backbone, err = digest(bbFrag, enz); err != nil {
		return &Frag{}, &Backbone{}, err
	}

	return
}

// getEnzymes return the enzyme with the name passed. errors out if there is none.
func (p *inputParser) getEnzyme(enzymeName string) (enzyme, error) {
	enzymeDB := NewEnzymeDB()
	if e, exists := enzymeDB.enzymes[enzymeName]; exists {
		return newEnzyme(e), nil
	}

	return enzyme{}, fmt.Errorf(
		`failed to find enzyme with name %s use "repp enzymes" for a list of recognized enzymes`,
		enzymeName,
	)
}

// getFilters takes an input string and returns a list of strings to run against matches
// when filtering out possible building fragments.
func (p *inputParser) getFilters(filterFlag string) []string {
	splitFunc := func(c rune) bool {
		return c == ' ' || c == ',' // space or comma separated
	}

	return strings.FieldsFunc(strings.ToUpper(filterFlag), splitFunc)
}

// read a FASTA file (by its path on local FS) to a slice of Fragments.
func read(path string, feature bool) (fragments []*Frag, err error) {
	if !filepath.IsAbs(path) {
		path, err = filepath.Abs(path)
		if err != nil {
			return nil, fmt.Errorf("failed to create path to input file: %s", err)
		}
	}

	dat, err := ioutil.ReadFile(path)
	if err != nil {
		return nil, err
	}
	file := string(dat)

	path = strings.ToLower(path)
	if strings.HasSuffix(path, "fa") ||
		strings.HasSuffix(path, "fasta") ||
		file[0] == '>' {
		return readFasta(path, file)
	}

	if strings.HasSuffix(path, "gb") ||
		strings.HasSuffix(path, "gbk") ||
		strings.HasSuffix(path, "genbank") {
		return readGenbank(path, file, feature)
	}

	return nil, fmt.Errorf("failed to parse %s: unrecognized file type", path)
}

// readFasta parses the multifasta file to fragments.
func readFasta(path, contents string) (frags []*Frag, err error) {
	// split by newlines
	lines := strings.Split(contents, "\n")

	// read in the frags
	var headerIndices []int
	var ids []string
	var fragTypes []fragType
	for i, line := range lines {
		if strings.HasPrefix(line, ">") {
			headerIndices = append(headerIndices, i)
			ids = append(ids, line[1:])
			if strings.Contains(line, "circular") {
				fragTypes = append(fragTypes, circular)
			} else {
				fragTypes = append(fragTypes, linear)
			}
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
		seq = strings.ToUpper(seq)
		seqs = append(seqs, seq)
	}

	// build and return the new frags
	for i, id := range ids {
		frags = append(frags, &Frag{
			ID:       id,
			Seq:      seqs[i],
			fragType: fragTypes[i],
		})
	}

	// opened and parsed file but found nothing
	if len(frags) < 1 {
		return frags, fmt.Errorf("failed to parse fragment(s) from %s", path)
	}

	return
}

// readGenbank parses a genbank file to fragments. Returns either fragments or parseFeatures,
// depending on the parseFeatures parameter.
func readGenbank(path, contents string, parseFeatures bool) (fragments []*Frag, err error) {
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

		features := []*Frag{}
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
			featureSeq = strings.ToUpper(featureSeq)

			labelRegex := regexp.MustCompile("\\/label=(.*)")
			labelMatch := labelRegex.FindStringSubmatch(feature)
			label := ""
			if len(labelMatch) > 1 {
				label = labelMatch[1]
			} else {
				label = strconv.Itoa(featureIndex)
			}

			features = append(features, &Frag{
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

	return []*Frag{
		&Frag{
			ID:  id,
			Seq: cleanedSeq,
		},
	}, nil
}

// igemBackbone returns a backbone, as it was in the database,
// if it corresponds to an iGEM backbone. They are not digested.
// see: http://parts.igem.org/Help:Prefix-Suffix
func igemBackbone(backbone string) bool {
	for _, bb := range []string{"pSB1A3", "pSB1T3", "pSB1K3", "pSB1C3"} {
		if bb == backbone {
			return true
		}
	}

	return false
}
