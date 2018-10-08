package assemble

import (
	"reflect"
	"testing"
)

func Test_node_distTo(t *testing.T) {
	type fields struct {
		id         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		other node
	}
	tests := []struct {
		name       string
		fields     fields
		args       args
		wantBpDist int
	}{
		{
			"negative distance with overlap",
			fields{
				id:         "1",
				uniqueID:   "1",
				start:      0,
				end:        40,
				assemblies: []assembly{},
			},
			args{
				other: node{
					id:         "2",
					uniqueID:   "2",
					start:      20,
					end:        60,
					assemblies: []assembly{},
				},
			},
			-20,
		},
		{
			"positive distance without overlap",
			fields{
				start: 0,
				end:   40,
			},
			args{
				other: node{
					start: 60,
					end:   100,
				},
			},
			20,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &node{
				id:         tt.fields.id,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
			}
			if gotBpDist := n.distTo(tt.args.other); gotBpDist != tt.wantBpDist {
				t.Errorf("node.distTo() = %v, want %v", gotBpDist, tt.wantBpDist)
			}
		})
	}
}

func Test_node_synthDist(t *testing.T) {
	conf.Synthesis.MaxLength = 100

	type fields struct {
		id         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		other node
	}
	tests := []struct {
		name           string
		fields         fields
		args           args
		wantSynthCount int
	}{
		{
			"synth dist is 0 with overlap",
			fields{
				start: 0,
				end:   40,
			},
			args{
				other: node{
					start: 20,
					end:   60,
				},
			},
			0,
		},
		{
			"synth dist is 1 without overlap",
			fields{
				start: 0,
				end:   40,
			},
			args{
				other: node{
					start: 60,
					end:   80,
				},
			},
			1,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &node{
				id:         tt.fields.id,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
			}
			if gotSynthCount := n.synthDist(tt.args.other); gotSynthCount != tt.wantSynthCount {
				t.Errorf("node.synthDist() = %v, want %v", gotSynthCount, tt.wantSynthCount)
			}
		})
	}
}

func Test_node_costTo(t *testing.T) {
	conf.Fragments.MinHomology = 20
	conf.PCR.BPCost = 0.03
	conf.Synthesis.BPCost = 0.05

	type fields struct {
		id         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		other node
	}
	tests := []struct {
		name     string
		fields   fields
		args     args
		wantCost float32
	}{
		{
			"just cost of PCR of new node if they overlap",
			fields{
				start: 0,
				end:   50,
			},
			args{
				other: node{
					start: 20,
					end:   100,
				},
			},
			40 * conf.PCR.BPCost,
		},
		{
			"cost of synthesis if they don't overlap",
			fields{
				start: 0,
				end:   50,
			},
			args{
				other: node{
					start: 80,
					end:   120,
				},
			},
			((float32(conf.Fragments.MinHomology) + 30.0) * conf.Synthesis.BPCost),
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &node{
				id:         tt.fields.id,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
			}
			if gotCost := n.costTo(tt.args.other); gotCost != tt.wantCost {
				t.Errorf("node.costTo() = %v, want %v", gotCost, tt.wantCost)
			}
		})
	}
}

func Test_node_reach(t *testing.T) {
	type fields struct {
		id         string
		uniqueID   string
		start      int
		end        int
		assemblies []assembly
	}
	type args struct {
		nodes      []node
		i          int
		synthCount int
	}
	tests := []struct {
		name          string
		fields        fields
		args          args
		wantReachable []node
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := &node{
				id:         tt.fields.id,
				uniqueID:   tt.fields.uniqueID,
				start:      tt.fields.start,
				end:        tt.fields.end,
				assemblies: tt.fields.assemblies,
			}
			if gotReachable := n.reach(tt.args.nodes, tt.args.i, tt.args.synthCount); !reflect.DeepEqual(gotReachable, tt.wantReachable) {
				t.Errorf("node.reach() = %v, want %v", gotReachable, tt.wantReachable)
			}
		})
	}
}
