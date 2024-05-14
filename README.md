# shapemapper2-mut-filtering
A script that takes --mutation-counts file from shapemapper2 and performs custom mutation signature filtering

*what it does*
This python script takes a modified and untreated mutation-counts file output from shapemapper2 along with a
mutation_type.csv filled out by the user and performs custom mutation signature filtering to calculate 2-8%
normalized reactivities.

*motivation*
Existing computational programs for analyzing mutational profiling data limit mutation signatures limit the data
to 1nt mismatches. With shapemapper2, the --dms flag will enact filtering that retains all 1nt mismatches for A, C,
and T, as well as G>C and G>T. RNAFramework allows for any combination of 1nt mismatches, but if filtering is
enacted for one nucleotide, all other nucleotides will also only count 1 nt mismatches. This script allows for custom
filtering and automatic reactivity calculations.

**Arguments***
-ut               untreated input file
-t                treated input file
output            the output file identifier (it is not necessary to specify file type)
-muttype          (optional) csv file that encodes for the mutation signature filtering to be applied
-primer           (optional) RT primer length for experiments, default is 10
--raw_react       outputs raw reactivity file (treated - untreated for all mutation types)
--raw_react_filt  outputs the raw reactivity file post filtering (treated - untreated for filtered mutation types)
--nind            flag indicating Normalizing the reactivity INDependently for each nucleotide.


*example data from Douds, Catherine A., Paul Babitzke, and Philip Bevilacqua. "A new reagent for in vivo structure probing of RNA G and U residues that improves RNA structure prediction alone and combined with DMS." RNA (2024): rna-079974.
