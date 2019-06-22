// Package cmd is for command line interactions with the repp application
package cmd

import (
	"log"

	"github.com/jjtimmons/repp/internal/repp"
	"github.com/spf13/cobra"
)

var (
	featureDB = repp.NewFeatureDB()

	enzymeDB = repp.NewEnzymeDB()
)

// RootCmd represents the base command when called without any subcommands.
var RootCmd = &cobra.Command{
	Use: "repp",
	Short: `PLAsmid DEfragger
	
Repository-based plasmid design. Specify and build plasmids using
their sequence, features, or fragments`,
	Version: "0.1.0",
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := RootCmd.Execute(); err != nil {
		log.Fatalf("%v", err)
	}
}
