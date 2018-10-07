package blast

import (
	"fmt"
	"path"
	"path/filepath"
	"testing"
)

// test that a primer's off-target is caught and flagged
func TestMismatch(t *testing.T) {
	db, _ = filepath.Abs(path.Join("..", "..", "test", "blast", "db"))
	primer := "AGTATAGGATAGGTAGTCATTCTT"
	parent := "gnl|addgene|107006"
	mismatchFound, m, err := Mismatch(primer, parent)

	if err != nil {
		t.Error(err)
		return
	}

	if !mismatchFound {
		t.Error("no mismatch caught")
		return
	}

	fmt.Printf("%+v", m)
}
