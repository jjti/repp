package defrag

import (
	"fmt"
	"io/ioutil"
	"log"
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
	in   string
	out  string
	dbs  []string
	igem bool
}

type inputParser struct{}

// parseFlags gathers the in path, out path, etc from the cobra cmd object
func parseFlags(cmd *cobra.Command) *flags {
	p := inputParser{}

	in, err := cmd.Flags().GetString("in")
	if err != nil {
		in = p.guessInput()
	}

	out, err := cmd.Flags().GetString("out")
	if err != nil {
		out = p.guessOutput(in)
	}

	addgene, err := cmd.Flags().GetBool("addgene")
	if err != nil {
		log.Fatal(err)
	}

	dbs, err := cmd.Flags().GetString("dbs")
	if err != nil && !addgene {
		log.Fatalf("failed to find dbs of building fragments: %v", err)
	}

	// read in the BLAST DB paths
	dbPaths, err := p.parseDBs(dbs, config.Root, addgene)
	if err != nil {
		log.Fatalf("failed to find a fragment database: %v", err)
	}

	return &flags{
		in:  in,
		out: out,
		dbs: dbPaths,
	}
}

// guessInput returns the first fasta file in the current directory. Is used
// if the user hasn't specified an input file
func (parser *inputParser) guessInput() (in string) {
	dir, _ := filepath.Abs(".")
	files, err := ioutil.ReadDir(dir)
	if err != nil {
		log.Fatal(err)
	}
	for _, file := range files {
		if file.IsDir() {
			continue
		}

		ext := strings.ToLower(filepath.Ext(file.Name()))
		if ext == ".fa" || ext == ".fasta" {
			return file.Name()
		}
	}

	log.Fatalf("failed: no input argument set and no fasta file found in %s", dir)
	return
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
	return
}
