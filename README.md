# DEFRAG

Declarative Engineering using Fragments, Repositories, and Gison
of
Declarative Engineering with FRAGments

## Steps

FRAG_COUNT_LIMIT = maximum number of fragments to include in the assembly

ASSEMBLY_LIST = list of assemblies that can be made circular (first node wraps around on itself across the zero-index)

BLACKLISTED = list of nodes that are deemed unusable (ex: inverted repeat in a junction region or primers with off-targets in the parent vector sequence)

1.  Build up ASSEMBLY_LIST

    1.1 Using reverse induction on a directed acyclic graph, traverse each node and build, on the node, a list of possible assemblies that could be used to get from that node to the "end of" the vector

        1.1.1 Add the current node to the start of all the assemblies on the "reachable" nodes and add to self's list of assemblies

        1.1.2 Only build up and add assemblies to self if the total number of fragments in the assembly are less than FRAG_COUNT_LIMIT

        1.1.3 Also estimate cost at this step, storing the estimated cost to get from this assembly to each of the children in the list (should be the cost from this node to the first in the assembly list + the cost of the first child in the assembly list)

2.  Find all the pareto optimal solutions

    2.1 Rank assemblies in lists for those involving a given number of solutions. So if the FRAG_COUNT_LIMIT is 5, there will be a cheapest solution with 5 fragments, a cheapest solution with 4 fragments, 3, etc

    2.2 Find the pareto optimal solutions: the assemblies with the fewest number of fragments and the cheapest estimated assembly cost. If the assembly with 3 fragments has a cheaper estimated assembly cost than the assembly with 4 fragments, the 4-fragment assembly is not pareto optimal and should not be included

3.  Traverse each pareto optimal solution from #3:

    3.1 Create primers if it's going to be PCR'ed, create a synthetic fragment otherwise

        3.1.1 Try and minimize the number of fragments with primers adding homology to adjacent fragments (they're going to be inferior primers). Could either do something where every other fragment gets homology-holding primers (spreading it out) or have adjacent fragments share in the homology (each gets half the homology)

                3.1.1.1 As an example: Fragment A and Fragment B coming together. We could add 20 additonal basepairs to the end of Fragment A and put nothing at the start of Fragment B... or we could add 10 basepairs from A to match B and 10 basepairs from A to match B (spreading out the shit primers). These primers would also be small enough so that we could pass a larger range to primer3 and let it figure out the best primers (as opposed to concatenating the sequence to its ends)

    3.2 Fail out (removing the assembly from ASSEMBLY_LIST) and repeat #2 if:

         3.2.1 The primers have a very high primer3 pair penalty score OR
               The primers have off-target's in their source vector OR
               The node has an inverted repeat in its junction OR
               The synthetic fragment will be dificult to synthesize

        3.2.2 The node has a duplicate end region with another fragment in the assembly

## List of red flags for building fragments

1.  Don't create primers for a fragment with off-targets

2.  Don't create primers with an excessive primer3 pair penalty

3.  Don't create synthetic fragments with synthesis complexities

4.  Don't create fragments with off-target end-homology

## Reachable fragments

Connect each node to all those that are "REACHABLE"

"REACHABLE" =
number_to_consider = Math.max(10, 0.05 \* number of fragments)
all fragments the current fragment overlaps with + number_to_consider

## Caveats

- Addgene verified sequence information doesn't include all the of the vector sequence for some of the vectors. Therefore there might be off-targets that are invisible

- Pareto optimal solutions right now are only so for cost and fragment count. If there are multiple assemblies with close the same cost/assembly-count, we might progress into creating primers as well and then compare primer3 pair penalty scores

## TODO

- Allow passing paths to multiple BLAST databases at once (like other vendors)
- Allow users to select databases, FASTA files with sequences they own (making cost zero)
