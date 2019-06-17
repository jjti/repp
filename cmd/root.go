// Package cmd is for command line interactions with the plade application
package cmd

import (
	"log"

	"github.com/jjtimmons/plade/internal/plade"
	"github.com/spf13/cobra"
)

var (
	featureDB = plade.NewFeatureDB()

	enzymeDB = plade.NewEnzymeDB()
)

// rootCmd represents the base command when called without any subcommands.
var rootCmd = &cobra.Command{
	Use: "plade",
	Short: `Build plasmids using DNA sequences available in public repositories.
Specify plasmids using their sequence, features, or fragments`,
	Version: "0.1.0",
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatalf("%v", err)
	}
}
