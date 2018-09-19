package io

import (
	"io/ioutil"
	"log"
	"regexp"
	"strings"

	"github.com/jjtimmons/decvec/internal/frag"
)

// ReadFASTA file to a slice of Fragments
func ReadFASTA(file string) []frag.Fragment {
	// read a file into memory
	dat, err := ioutil.ReadFile(file)
	if err != nil {
		log.Fatalln("failed to open fasta file", err)
	}

	// read it into a string
	file := string(dat)

	// split by newlines
	lines := strings.Split(file, "\n")

	// read in the fragments
	var headerIndices []int
	ids := []string{}
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
	var fragments []frag.Fragment
	for i, id := range ids {
		fragments = append(fragments, frag.Fragment{
			ID:  id,
			Seq: seqs[i],
		})
	}

	// opened and parsed file but found nothing
	if len(fragments) < 1 {
		log.Fatalln("failed to parse building fragments from %s", file)
	}

	return fragments
}
