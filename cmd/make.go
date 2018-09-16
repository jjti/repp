package cmd

import (
	"encoding/json"
	"fmt"

	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var targetPath string
var inputFastaPath string

// makeCmd represents the make command
var makeCmd = &cobra.Command{
	Use:   "make",
	Short: "Make a vector from its target sequence and a database of building fragments",
	Long: `Make a vector from its target sequence and a database of building fragments
	
"decvec make" is for assembling a vector, using Gibson Assembly, in the most efficient 
way possible. It's a declaractive approach to vector design. This means that, rather than 
telling make which DNA fragments go together, and in what order, make simply builds the 
vector it's told to. It does this by:

1. Screening the fragments in a database of existing fragments, "db-fasta," to find
   fragments that might be useful for creating the vector
2. Creating a list of possible assemblies using fragments from 1, and ranking them by
   estimated assembly cost
3. Finding the cheapest assembly-plan with the best Gibson Assembly fragments possible`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("make called")
		b, err := json.MarshalIndent(viper.AllSettings(), "", "  ")
		if err == nil {
			fmt.Print(string(b))
		}
	},
}

func init() {
	rootCmd.AddCommand(makeCmd)

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// makeCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Flags for specifying the paths to the input file, input fragment files, and output file
	makeCmd.Flags().StringVarP(&targetPath, "target", "t", "", "path to a FASTA file with a target sequence to make")
	makeCmd.Flags().StringVarP(&inputFastaPath, "db-fasta", "f", "", "path to a multi-FASTA file with building sequences")

	// Mark required flags
	makeCmd.MarkFlagRequired("target")
	makeCmd.MarkFlagRequired("db-fasta")

	// Bind the paramters to viper
	viper.BindPFlag("target", makeCmd.Flags().Lookup("target"))
	viper.BindPFlag("db-fasta", makeCmd.Flags().Lookup("db-fasta"))
}
