package assemble

import (
	"reflect"
	"testing"
)

func Test_build(t *testing.T) {
	type args struct {
		nodes []node
	}
	tests := []struct {
		name         string
		args         args
		wantCountMap map[int][]assembly
	}{
		// TODO: write table tests
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotCountMap := build(tt.args.nodes); !reflect.DeepEqual(gotCountMap, tt.wantCountMap) {
				t.Errorf("build() = %v, want %v", gotCountMap, tt.wantCountMap)
			}
		})
	}
}
