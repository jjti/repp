package main

import (
	"fmt"
	"os"

	"github.com/jjtimmons/repp/cmd"
	"github.com/spf13/cobra/doc"
)

func main() {
	if len(os.Args) > 1 && os.Args[1] == "docs" {
		if err := doc.GenMarkdownTree(cmd.RootCmd, "./docs"); err != nil {
			fmt.Println(err.Error())
		} else {
			// need to update the default repp entry to be the doc entry page
			os.Rename("./docs/repp.md", "./docs/index.md")
		}
	} else {
		cmd.Execute() // initialize cobra commands
	}
}
