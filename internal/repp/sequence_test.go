package repp

import (
	"path"
	"testing"
)

// if an input fragment being built is exactly the same as one in a DB, it should be used
// as is and without PCR or any preparation
func Test_sequence(t *testing.T) {
	fs, c := NewFlags(
		path.Join("..", "..", "test", "input", "BBa_K1649003.fa"),
		path.Join("..", "..", "test", "output", "BBa_K1649003.json"),
		"pSB1A3",
		"2015,2016,2017,2018",
		[]string{"PstI"},
		[]string{},
		false,
		true,
		false,
	)

	results := Sequence(fs, c) // use addgene database

	if len(results) < 1 {
		t.Fail()
	}
}
