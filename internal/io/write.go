package io

import (
	"encoding/json"
	"io/ioutil"
	"log"

	"github.com/jjtimmons/decvec/internal/dvec"
)

// Write a slice of possible assemblies to the fs at the output path
func Write(filename string, assemblies [][]dvec.Fragment) {
	output, err := json.MarshalIndent(assemblies, "", "  ")
	if err != nil {
		log.Fatalf("Failed to serialize the output data: %v", err)
	}

	err = ioutil.WriteFile(filename, output, 0666)
	if err != nil {
		log.Fatalf("Failed to write the results to the file system: %v", err)
	}
}
