# GUESS

Great Until Every Sequence is Synthesized

## Steps

FRAG_COUNT_LIMIT = maximum number of fragments to include in the assembly

ASSEMBLY_LISTS = map from fragment count to the list of assemblies, as a list of nodes, that will make an assembly with that number of nodes (in the final assembly). Sorted by estimated assembly cost

BLACKLISTED_NODES = list of nodes that are deemed unusable (ex: inverted repeat in a junction region or primers that have off-targets in the parent vector sequence)

1.  Build up ASSEMBLY_LISTS

    1.1 Using reverse induction on a directed acyclic graph, traverse each node and build, on the node, a list of possible assemblies that could be used to get from that node to the "end of" the vector

        1.1.1 Add the current node to the start of all the assemblies on the "reachable" nodes and add to self's list of assemblies

        1.1.2 Only build up and add assemblies to self if the total number of fragments in the assembly are less than FRAG_COUNT_LIMIT

        1.1.3 Also estimate cost at this step, storing the estimated cost to get from this assembly to each of the children in the list (should be the cost from this node to the first in the assembly list + the cost of the first child in the assembly list)

    1.2 After building up all assembly lists on all nodes, gather all the assembly lists from the "entry nodes." (ie, nodes that are a seq-length away from the end of the target vector sequence)

2.  Find all the pareto optimal solutions

    2.1 Rank assemblies in lists for those involving a given number of solutions. So if the FRAG_COUNT_LIMIT is 5, there will be a cheapest solution with 5 fragments, a cheapest solution with 4 fragments, 3, etc

    2.2 Find the pareto optimal solutions: the assemblies that have the fewest number of fragments and the cheapest estimated assembly cost. If the assembly with 3 fragments has a cheaper estimated assembly cost than the assembly with 4 fragments, the 4-fragment assembly is not pareto optimal and should not be included

    2.3 Remove assemblies that include a node from BLACKLISTED_NODES

3.  Traverse each pareto optimal solution from #3:

    3.1 Create primers if it's going to be PCR'ed, create a synthetic fragment otherwise

    3.2 Fail out (removing the assembly from ASSEMBLY_LISTS) and repeat #2 if:

         3.2.1 The primers have a very high primer3 pair penalty score OR
               The primers have off-target's in their source vector OR
               The node has an inverted repeat in its junction OR
               The synthetic fragment will be dificult to synthesize

               3.2.1.1 If its a PCR Fragment, add it to BLACKLISTED_NODES

        3.2.2 The node has a duplicate end region with another fragment in the assembly

## List of things that need to be considered when creating building fragments

1.  Don't create primers for a fragment if they have off-targets in the parent fragment

2.  Don't create primers that have an excessive primer3 pair penalty

3.  Don't create synthetic fragments that have high synthesis complexities

4.  Don't create fragments that have off-target end-homology

## Reachable fragments

Connect each node to all those that are "REACHABLE"

"REACHABLE" =
number_to_consider = Math.max(10, 0.05 \* number of fragments)
all fragments the current fragment overlaps with + number_to_consider

## Caveats

- Addgene verified sequence information doesn't include all the of the vector sequence for some of the vectors. Therefore there might be off-targets that are invisible

- Pareto optimal solutions right now are only so for cost and fragment count. If there are multiple assemblies with close the same cost/assembly-count, we might progress into creating primers as well and then comparing on primer3 pair penalty score
