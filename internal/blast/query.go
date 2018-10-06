package blast

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"os/exec"
	"path"
)

// Query the BLAST db to get the template sequence (needed to avoid
// off-target primers later on)
func Query(entry string, db string) (string, error) {
	// path to the blastcmd binary
	blastcmd := path.Join("..", "..", "vendor", "ncbi-blast-2.7.1+", "bin", "blastcmd")
	// path to the output sequence file
	out := path.Join(blastDir, entry+".query")

	// make a blast command
	queryCmd := exec.Command(
		blastcmd,
		"-task", "blastn",
		"-db", db,
		"-dbtype", "nucl",
		"-entry", entry,
		"-out", out,
		// %s means sequence data (without defline)
		"-outfmt", "%s",
	)

	var stderr bytes.Buffer
	queryCmd.Stderr = &stderr

	// query the blast db
	err := queryCmd.Run()
	if err != nil {
		return "", fmt.Errorf("failed to execute blastcmd for entry %s: %v: %s", entry, err, stderr.String())
	}

	// read the file into a sequence
	file, err := ioutil.ReadFile(out)
	if err != nil {
		return "", fmt.Errorf("failed to read blastcmd output file for entry %s: %v", entry, err)
	}
	fileS := string(file)

	// return the sequence associated with the fragment in the DB
	return fileS, nil
}
