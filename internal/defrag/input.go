package defrag

import (
	"fmt"
	"io/ioutil"
	"os"
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
	backbone Frag
}

// inputParser contains methods for parsing flags from the input &cobra.Command
type inputParser struct{}

// newFlags makes a new flags object out of the in, out, dbs passed
//
// just for testing right now, but may also be useful for AWS lambda down the line
func newFlags(in, out string, dbs []string, addgene bool) *flags {
	dbPaths := dbs
	if addgene {
		dbPaths = append(dbPaths, path.Join(config.Root, "assets", "addgene", "db", "addgene"))
	}

	return &flags{
		in:  in,
		out: out,
		dbs: dbPaths,
	}
}

// parseFlags gathers the in path, out path, etc from the cobra cmd object
func parseFlags(cmd *cobra.Command, conf *config.Config) (parsedFlags *flags, err error) {
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

	// check if user asked for a specific backbone, confirm it exists in one of the dbs
	backbone, _ := cmd.Flags().GetString("backbone")

	// check if they also specified an enzyme
	enzymeName, _ := cmd.Flags().GetString("enzyme")

	// try to digest the backbone with the enzyme
	if backbone != "" {
		if enzymeName == "" {
			return nil, fmt.Errorf("backbone passed, %s, without an enzyme to digest it", backbone)
		}

		// confirm that the backbone exists in one of the dbs (or local fs) gather it as a Frag if it does
		backbone, err := p.getBackbone(backbone, parsedFlags.dbs, conf)
		if err != nil {
			return nil, err
		}

		// gather the enzyme by name, err if it's unknown
		enzyme, err := p.getEnzyme(enzymeName)
		if err != nil {
			return nil, err
		}

		if parsedFlags.backbone, err = digest(backbone, enzyme); err != nil {
			return nil, err
		}
	}

	return
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

		ext := strings.ToLower(filepath.Ext(file.Name()))
		if ext == ".fa" || ext == ".fasta" {
			return file.Name(), nil
		}
	}

	return "", fmt.Errorf("failed: no input argument set and no fasta file found in %s", dir)
}

// guessOutput gets an outpath path from an input path (if no output path is
// specified). It uses the same name as the input path to create an output
func (p *inputParser) guessOutput(in string) (out string) {
	ext := filepath.Ext(in)
	noExt := in[0 : len(in)-len(ext)]
	return noExt + ".defrag.json"
}

// parseDBs returns a list of absolute paths to BLAST databases
func (p *inputParser) parseDBs(dbsInput, root string, addgene bool) (paths []string, err error) {
	dbs := dbsInput
	if addgene {
		dbs += path.Join(root, "assets", "addgene", "db", "addgene")
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

// getBackbone is for confirming that the fragment passed actually exists in
// one of the databases. Avoid nonsense backbones
//
// TODO: use goroutine
// TODO: test
func (p *inputParser) getBackbone(backbone string, dbs []string, c *config.Config) (f Frag, err error) {
	// first try to get the backbone out of a local file
	if frags, err := read(backbone); err == nil && len(frags) > 0 {
		return frags[0], nil // it was a local file
	}

	// move through each db and see if it contains the backbone
	for _, db := range dbs {
		// if outFile is defined here we managed to query it from the db
		if outFile, _ := blastdbcmd(backbone, db, c); outFile != "" {
			frags, err := read(outFile)
			frags[0].Type = circular // assume its circular here, used as backbone
			return frags[0], err
		}
	}

	dbMessage := strings.Join(dbs, "\n")
	return Frag{}, fmt.Errorf("failed to find backbone %s in local filesystem or any of:\n%s", backbone, dbMessage)
}

// getEnzymes return the enzyme with the name passed. errors out if there is none
func (p *inputParser) getEnzyme(enzymeName string) (enzyme, error) {
	if e, exists := enzymes[enzymeName]; exists {
		return e, nil
	}

	return enzyme{}, fmt.Errorf(`failed to find enzyme with name %s
use "defrag enzymes" for a list of recognized enzymes`, enzymeName)
}
