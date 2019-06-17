package main

import (
	_ "net/http/pprof"

	"github.com/jjtimmons/plade/cmd"
)

func main() {
	cmd.Execute() // initialize cobra commands
	// log.Println(http.ListenAndServe("localhost:6060", nil)) // for profiling
}
