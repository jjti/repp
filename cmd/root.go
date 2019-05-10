// Package cmd is for command line interactions with the rvec application
package cmd

import (
	"log"

	"github.com/jjtimmons/rvec/internal/rvec"
	"github.com/spf13/cobra"
)

var (
	featureDB = rvec.NewFeatureDB()

	enzymeDB = rvec.NewEnzymeDB()
)

// rootCmd represents the base command when called without any subcommands.
var rootCmd = &cobra.Command{
	Use: "rvec",
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
