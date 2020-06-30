# MTMAS: Multi-trait MArker Selection

MTMAS is a procedure and R script to optimize explained variance across multiple traits with a smallest possible set of markers. It uses a genetic algorithm and a penalty parameter to tune how small the marker set should be.

MTMAS is not an R package, but just a script that contains the mtmas() function. You can 'source' it directly from github, or download and source it from a local drive (see example). 
The mtmas() function needs as input files GWAS summary statistics for the traits to be optimized, the mtmas() function was designed based on GAPIT output, but it is quite easy to format other outputs to match (see below). The second impportant input file is a file with LD information between the markers - also see below how you need to compute it in the form that mtmas() expects.

If you have improvements on the script, you are welcome to make a pull request and I will integrate it in the version stored here.

# Installation

