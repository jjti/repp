package main

import (
	"os"

	"github.com/jjtimmons/repp/cmd"
)

func main() {
	if len(os.Args) > 1 && os.Args[1] == "docs" {
		makeDocs() // create the documentation pages. see ./docs.go
	} else {
		cmd.Execute() // initialize cobra commands
	}
}
