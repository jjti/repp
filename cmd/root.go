package cmd

import (
	"log"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "defrag",
	Short: "Build vectors from any of their vector sequence, features or fragments",
	Long:  ``,
}

// backboneHelp is the help message for the backbone CLI argument
var backboneHelp = `backbone to insert the fragments into. Can either be an entry 
in one of the dbs or a file on the local filesystem.`

var enzymeHelp = `enzyme to linearize the backbone with (backbone must be specified).
'defrag enzymes' prints a list of recognized enzymes.`

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatalf("Failed to execute defrag: %v", err)
	}
}

func init() {
	// settings is an optional parameter for a settings file (that overrides the fields in BaseSettingsFile)
	rootCmd.PersistentFlags().StringP("settings", "s", config.BaseSettingsFile, "defrag settings file")
	viper.BindPFlag("settings", rootCmd.PersistentFlags().Lookup("settings"))
}
