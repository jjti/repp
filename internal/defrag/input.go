package defrag

import (
	"fmt"
	"io/ioutil"
	"path"
	"path/filepath"
	"regexp"
	"strings"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
)

// flags conatins parsed cobra flags like "in", "out", "dbs", etc that
// are used by multiple commands
type flags struct {
	in       string
	out      string
	dbs      []string
	backbone string
	enzyme   string
}

type inputParser struct{}

// parseFlags gathers the in path, out path, etc from the cobra cmd object
func parseFlags(cmd *cobra.Command) (parsedFlags *flags, err error) {
	p := inputParser{}

	if parsedFlags.in, err = cmd.Flags().GetString("in"); err != nil {
		// check whether an input fail was specified
		if parsedFlags.in, err = p.guessInput(); err != nil {
			// try to guess at the input file the user wanted to use
			return
		}
	}

	if parsedFlags.out, err = cmd.Flags().GetString("out"); err != nil {
		parsedFlags.out = p.guessOutput(parsedFlags.in)
	}

	addgene, err := cmd.Flags().GetBool("addgene")
	if err != nil {
		return nil, fmt.Errorf("failed to parse addgene flag: %v", err)
	}

	dbString, err := cmd.Flags().GetString("dbs")
	if err != nil && !addgene {
		return nil, fmt.Errorf("failed to parse building fragments: %v", err)
	}

	// read in the BLAST DB paths
	if parsedFlags.dbs, err = p.parseDBs(dbString, config.Root, addgene); err != nil || len(parsedFlags.dbs) == 0 {
		return nil, fmt.Errorf("failed to find any fragment databases: %v", err)
	}

	// check if user asked for a specific backbone
	parsedFlags.backbone, _ = cmd.Flags().GetString("backbone")

	// check if they also specified an enzyme
	parsedFlags.enzyme, _ = cmd.Flags().GetString("enzyme")

	// enzyme should only be specified if backbone is
	if parsedFlags.enzyme != "" && parsedFlags.backbone == "" {
		return nil, fmt.Errorf("enzymes can only be chosen if a backbone is as well")
	}

	return
}

// guessInput returns the first fasta file in the current directory. Is used
// if the user hasn't specified an input file
func (parser *inputParser) guessInput() (in string, err error) {
	dir, _ := filepath.Abs(".")
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		return
	}

	for _, file := range files {
		if file.IsDir() {
			continue
		}

		ext := strings.ToLower(filepath.Ext(file.Name()))
		if ext == ".fa" || ext == ".fasta" {
			return file.Name(), nil
		}
	}

	return "", fmt.Errorf("failed: no input argument set and no fasta file found in %s", dir)
}

// guessOutput gets an outpath path from an input path (if no output path is
// specified). It uses the same name as the input path to create an output
func (parser *inputParser) guessOutput(in string) (out string) {
	ext := filepath.Ext(in)
	noExt := in[0 : len(in)-len(ext)]
	return noExt + ".defrag.json"
}

// parseDBs returns a list of absolute paths to BLAST databases
func (parser *inputParser) parseDBs(dbsInput, root string, addgene bool) (paths []string, err error) {
	if addgene {
		addgenePath := path.Join(root, "assets", "addgene", "db", "addgene")
		return parser.dbPaths(dbsInput + "," + addgenePath)
	}

	return parser.dbPaths(dbsInput)
}

// dbPaths turns a single string of comma separated BLAST dbs into a
// slice of absolute paths to the BLAST dbs on the local fs
func (parser *inputParser) dbPaths(dbList string) (paths []string, err error) {
	re := regexp.MustCompile(",")

	for _, db := range re.Split(dbList, -1) {
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

	// TODO: add a test to make sure all the blast databases exist on the FS

	return
}
