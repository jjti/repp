// Package blast is for finding building fragments that may be able to contribute
// to the target sequence the user is trying to construct
package blast

import (
	"os"
	"path/filepath"
)

// init ensures there's a blast subdirectory in the binary's execution enviornment
// for the BLAST database that this is about to create
func init() {
	blastPath := filepath.Join(".", "blast")
	os.MkdirAll(blastPath, os.ModePerm)
}
