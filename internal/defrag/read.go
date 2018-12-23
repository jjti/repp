package defrag

import (
	"fmt"
	"io/ioutil"
	"path/filepath"
	"regexp"
	"strings"
)

// read a FASTA file (by its path on local FS) to a slice of Fragments
func read(path string) (fragments []Fragment, err error) {
	if !filepath.IsAbs(path) {
		path, err = filepath.Abs(path)
		if err != nil {
			return nil, fmt.Errorf("Failed to create path to FASTA file: %s", err)
		}
	}

	dat, err := ioutil.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("Failed to read input FASTA path: %s", err)
	}
	file := string(dat)

	// split by newlines
	lines := strings.Split(file, "\n")

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
		fragments = append(fragments, Fragment{
			ID:  id,
			Seq: seqs[i],
		})
	}

	// opened and parsed file but found nothing
	if len(fragments) < 1 {
		return fragments, fmt.Errorf("Failed to parse fragment(s) from %s", file)
	}

	return
}
