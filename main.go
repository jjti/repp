package main

import (
	"github.com/jjtimmons/rvec/cmd"
)

func main() {
	cmd.Execute() // initialize cobra commands

	// log.Println(http.ListenAndServe("localhost:6060", nil))
}
