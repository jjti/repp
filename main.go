package main

import (
	"fmt"

	"github.com/jjtimmons/plade/cmd"
	"github.com/spf13/cobra/doc"
)

func main() {
	// cmd.Execute() // initialize cobra commands
	// log.Println(http.ListenAndServe("localhost:6060", nil)) // for profiling

	err := doc.GenMarkdownTree(cmd.RootCmd, "./docs")
	if err != nil {
		fmt.Println(err.Error())
	}

	err = doc.GenManTree(cmd.RootCmd, &doc.GenManHeader{
		Title:   "plade",
		Section: "1",
	}, "./docs")

	if err != nil {
		fmt.Println(err.Error())
	}
}
