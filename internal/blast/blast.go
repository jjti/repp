package blast

import (
	"os"
	"path/filepath"
)

// init ensures there's a blast subdirectory in the binary's execution enviornment
// for the BLAST database this is about to create
func init() {
	blastPath := filepath.Join(".", "blast")
	os.MkdirAll(blastPath, os.ModePerm)
}
