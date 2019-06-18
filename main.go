package main

import (
	"github.com/jjtimmons/plade/cmd"
)

func main() {
	cmd.Execute() // initialize cobra commands
	// log.Println(http.ListenAndServe("localhost:6060", nil)) // for profiling

	// if err := doc.GenMarkdownTree(cmd.RootCmd, "./docs"); err != nil {
	// 	fmt.Println(err.Error())
	// } else {
	// 	// need to update the default plade entry to be the doc entry page
	// 	os.Rename("./docs/plade.md", "./docs/index.md")
	// }
}
