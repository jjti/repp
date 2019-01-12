package cmd

import (
	"log"
	"path"
	"path/filepath"

	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "defrag",
	Short: "Build vectors from any of their vector sequence, features or fragments",
	Long: `Entry name of a backbone to insert the fragments into. Must match an entry 
	in one of the dbs (either those passed manually or in AddGene, iGEM). If an
	enzyme chosen, the backbone will be linearized with that enzyme and
	the largest resulting fragment will be used as the backbone. If no enzyme
	is specified, defrag will chose one nearest the first bp of the backbone with a
	single recognition site`,
}

// backboneHelp is the help message for the backbone CLI argument
var backboneHelp = `Backbone to insert the fragments into. Can either be an entry 
	in one of the dbs or a FASTA file on the local filesystem. An enzyme must also
	be chosen to linearize the backbone. Pre-linearized backbones are not yet
	supported`

var enzymeHelp = `Enzyme to linearize the backbone with (backbone must be specified).
	The enzyme's name must be recognizable by defrag. Use 'defrag enzymes' for a list
	of recognized enzyme names and their recognition sequences`

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatalf("Failed to execute defrag: %v", err)
	}
}

func init() {
	settings, _ := filepath.Abs(path.Join("..", "config"))

	// Here you will define your flags and configuration settings.
	// Cobra supports persistent flags, which, if defined here,
	// will be global for your application.
	rootCmd.PersistentFlags().StringP("config", "c", settings, "config file")
	viper.BindPFlag("config", rootCmd.PersistentFlags().Lookup("config"))
}
