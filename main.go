package main

import (
	"fmt"

	"github.com/jjtimmons/plade/cmd"
	"github.com/spf13/cobra/doc"
)

func main() {
	// cmd.Execute() // initialize cobra commands
	// log.Println(http.ListenAndServe("localhost:6060", nil)) // for profiling

	if err := doc.GenMarkdownTree(cmd.RootCmd, "./docs"); err != nil {
		fmt.Println(err.Error())
	}

	title := &doc.GenManHeader{
		Title:   "plade",
		Section: "1",
	}
	if err := doc.GenManTree(cmd.RootCmd, title, "./docs"); err != nil {
		fmt.Println(err.Error())
	}
}
