package main

import (
	"github.com/jjtimmons/repp/cmd"
)

func main() {
	cmd.Execute() // initialize cobra commands
	// log.Println(http.ListenAndServe("localhost:6060", nil)) // for profiling

	// if err := doc.GenMarkdownTree(cmd.RootCmd, "./docs"); err != nil {
	// 	fmt.Println(err.Error())
	// } else {
	// 	// need to update the default repp entry to be the doc entry page
	// 	os.Rename("./docs/repp.md", "./docs/index.md")
	// }
}
