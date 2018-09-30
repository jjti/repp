# decvec

declarative vector design

## Steps

1. Filter out fragments that aren't going to be in assemblies beneath the upper-fragment count limit
2. Run a DP algo on the nodes to find each node's cost from end
3. Sort all "starting-nodes" by their estimated total assembly cost (low -> high)
4. Traverse the cheapest assembly and "fill in the node"
