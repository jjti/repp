package blast

import (
	"path"
	"path/filepath"
	"testing"
)

func TestMismatch(t *testing.T) {
	db, _ = filepath.Abs(path.Join("..", "..", "test", "blast", "db"))
	_, err := Mismatch(
		"CAAATTTTGTAATCCAGAGGTTGATTATCGATAAGCTTGATATCGAATTCATAACTTCGTATAGCTACATTATACGAAGTTATTCT",
		"gnl|addgene|85039.2",
	)

	if err != nil {
		t.Error(err)
	}
}
