// Package config is for app wide settings that are unmarshalled
// from Viper (see: /cmd)
package config

import (
	"reflect"
	"testing"
)

func TestNew(t *testing.T) {
	tests := []struct {
		name  string
		wantC Config
	}{
		{
			"create a map from fragment size upper limit to cost",
			Config{},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotC := New(); !reflect.DeepEqual(gotC, tt.wantC) {
				t.Errorf("New() = %+v, want %+v", gotC, tt.wantC)
			}
		})
	}
}
