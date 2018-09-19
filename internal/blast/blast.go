// Package blast is for finding building fragments that may be able to contribute
// to the target sequence the user is trying to construct
package blast

import (
	"bytes"
	"errors"
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/jjtimmons/decvec/config"
	"github.com/jjtimmons/decvec/internal/frag"
)

// blastExec is a small utility function for executing BLAST
// on a fragment
type blastExec struct {
	// the fragment we're BLASTing
	f *frag.Fragment

	// the path to the input BLAST file
	in string

	// the path the output BLAST file
	out string

	// the path to the blast executable
	blast string

	// path to the database to blast against
	db string
}

// init ensures there's a blast subdirectory in the binary's execution enviornment
// for the BLAST database that this is about to create
func init() {
	blastPath := filepath.Join(".", "blast")
	os.MkdirAll(blastPath, os.ModePerm)
}

// BLAST the passed Fragment against a set from the command line and create
// matches for those that are long enough
func BLAST(f *frag.Fragment) error {
	c := config.NewConfig()

	// get current path
	exPath, err := os.Executable()
	if err != nil {
		return fmt.Errorf("failed to get binary's path: ", err)
	}

	// exPath is to the binary, step up one
	exPath = path.Join(exPath, "..")

	// create the utility struct with paths
	b := blastExec{
		f:     f,
		in:    exPath + ".input.fa",
		out:   exPath + ".output.json",
		blast: path.Join(exPath, "..", "tools", "ncbi-blast-2.7.1+", "bin", "blastn"),
		db:    c.Make.DBPath,
	}

	// create the input file
	if err = input(b); err != nil {
		return fmt.Errorf("failed at creating BLAST input file at %s: %v", b.in, err)
	}

	// execute BLAST on it
	if err = run(b); err != nil {
		return fmt.Errorf("failed at executing BLAST: %v", err)
	}

	// parse the BLAST output to Matches for the Fragment
	matches, err := parse(b)
	if err != nil {
		return fmt.Errorf("failed to parse BLAST output: %v", err)
	}
	f.Matches = matches

	return nil
}

// input is for creating an input file for BLAST
// return the path to the file and an error if there was one
func input(b blastExec) error {
	// create the file contents, double sequence because its circular
	fileContents := ">" + b.f.ID + "\n" + b.f.Seq + b.f.Seq

	// create file
	inputFile, err := os.Create(b.in)

	// close at the end of execution
	defer inputFile.Close()

	// write to it
	_, err = inputFile.WriteString(fileContents)
	if err != nil {
		return err
	}

	return nil
}

// run is for calling BLAST and returning when it finished
// returns the local path to the output fils
func run(b blastExec) error {
	// error output
	var errOut bytes.Buffer

	// create the blast command
	// https://www.ncbi.nlm.nih.gov/books/NBK279682/
	blastCmd := &exec.Cmd{
		Path: b.blast,
		Args: []string{
			"-db " + b.db,
			"-query " + b.in,
			"-out " + b.out,
			"-outfmt \"7 qseqid qstart qend sseqid sstart send\" ",
			"-ungapped ",
			"-perc_identity 100.0",
		},
		Stderr: &errOut,
	}

	// execute BLAST and wait on it to finish
	err := blastCmd.Run()
	if err != nil {
		return err
	}

	// check err results
	if errOut.Len() > 0 {
		return errors.New(errOut.String())
	}

	return nil
}

// parse is for reading the output file into Matches on the Fragment
// returns a slice of Matches for the blasted fragment
func parse(b blastExec) ([]frag.Match, error) {
	// read in the results
	file, err := ioutil.ReadFile(b.out)
	fileS := string(file)
	if err != nil {
		return nil, err
	}

	// read it into Matches
	var ms []frag.Match
	for _, line := range strings.Split(fileS, "\n") {
		// comment lines start with a #
		if strings.HasPrefix(line, "#") {
			continue
		}

		cols := strings.Split(line, " ")

		if len(cols) != 6 {
			continue
		}

		start, _ := strconv.Atoi(cols[1])
		end, _ := strconv.Atoi(cols[2])

		// create and append the new match
		ms = append(ms, frag.Match{
			ID:    cols[3],
			Start: start,
			End:   end,
		})
	}
	return ms, nil
}
