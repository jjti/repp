package defrag

import (
	"path"
	"strings"
	"testing"
)

// if an input fragment being built is exactly the same as one in a DB, it should be used
// as is and without PCR or any preparation
func Test_vector_single_vector(t *testing.T) {
	fs, c := NewFlags(
		path.Join("..", "..", "test", "input", "109049.addgene.fa"),
		path.Join("..", "..", "test", "output", "109049.output.json"),
		"",
		"",
		"",
		[]string{},
		true,
		false,
		false,
	)

	_, _, assemblies, err := sequence(fs, c) // use addgene database
	if err != nil {
		t.Error(err)
	}

	if !strings.Contains(assemblies[0][0].ID, "109049") && !strings.Contains(assemblies[1][0].ID, "109049") {
		t.Fatal("failed to use 109049 to build the vector")
	}
}
