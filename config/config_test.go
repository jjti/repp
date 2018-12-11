// Package config is for app wide settings that are unmarshalled
// from Viper (see: /cmd)
package config

import (
	"fmt"
	"path/filepath"
	"reflect"
	"testing"
)

func TestConfig_SynthCost(t *testing.T) {
	type fields struct {
		DBs       string
		Fragments FragmentConfig
		PCR       PCRConfig
		Synthesis SynthesisConfig
		filled    bool
	}
	type args struct {
		fragLength int
	}

	configFields := fields{
		DBs:       "",
		Fragments: FragmentConfig{},
		PCR:       PCRConfig{},
		Synthesis: SynthesisConfig{
			Cost: map[int]SynthCost{
				50: SynthCost{
					Fixed:   true,
					Dollars: 5,
				},
				200: SynthCost{
					Fixed:   false,
					Dollars: 0.1,
				},
				10000: SynthCost{
					Fixed:   false,
					Dollars: 0.50,
				},
			},
			MaxLength: 1000,
			MinLength: 20,
		},
		filled: true,
	}

	tests := []struct {
		name   string
		fields fields
		args   args
		want   float32
	}{
		{
			"fixed cost synthesis",
			configFields,
			args{
				fragLength: 20,
			},
			5.0,
		},
		{
			"variable cost synthesis (small)",
			configFields,
			args{
				fragLength: 100,
			},
			10.0,
		},
		{
			"variable cost synthesis (large)",
			configFields,
			args{
				fragLength: 1000,
			},
			500.0,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			c := Config{
				DBs:       tt.fields.DBs,
				Fragments: tt.fields.Fragments,
				PCR:       tt.fields.PCR,
				Synthesis: tt.fields.Synthesis,
				filled:    tt.fields.filled,
			}
			if got := c.SynthCost(tt.args.fragLength); got != tt.want {
				t.Errorf("Config.SynthCost() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_parseDBs(t *testing.T) {
	db1 := "../exampleDir/exampleDB.fa"
	dbAbs1, _ := filepath.Abs(db1)

	db2 := "otherBLASTDir/otherDB.fa"
	dbAbs2, _ := filepath.Abs(db2)

	type args struct {
		dbList string
	}
	tests := []struct {
		name      string
		args      args
		wantPaths []string
		wantError error
	}{
		{
			"single blast path",
			args{
				dbList: db1,
			},
			[]string{dbAbs1},
			nil,
		},
		{
			"multi fasta db paths",
			args{
				dbList: fmt.Sprintf("%s, %s", db1, db2),
			},
			[]string{dbAbs1, dbAbs2},
			nil,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotPaths, _ := parseDBs(tt.args.dbList); !reflect.DeepEqual(gotPaths, tt.wantPaths) {
				t.Errorf("parseDBs() = %v, want %v", gotPaths, tt.wantPaths)
			}
		})
	}
}
