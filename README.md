## MTMAS: Multi-trait MArker Selection

MTMAS is a procedure and R script to optimize explained variance across multiple traits with a smallest possible set of markers. It uses a genetic algorithm and a penalty parameter to tune how small the marker set should be.

MTMAS is not an R package, but just a script that contains the mtmas() function. You can 'source' it directly from github, or download and source it from a local drive (see example). 
The mtmas() function needs as input files GWAS summary statistics for the traits to be optimized, the mtmas() function was designed based on GAPIT output, but it is quite easy to format other outputs to match (see below). The second important input file is a file with LD information between the markers - also see below how you need to compute it in the form that mtmas() expects.

There is example data provided with GWAS summary statistics from an analysis of 5 traits.

If you have improvements on the script, you are welcome to make a pull request and I will integrate it in the version stored here.

### Installation

To use mtmas, you only have to 'source' the mtmas.R script like:

```R
source("https://github.com/ljanss/mtmas.R")
```

### Input data

The GWAS summary statistics are provided as one file per trait. We have developed the mtmas() function based on GAPIT output, but any other GWAS output can be used when it has at least the following three columns with matching names:
- SNP: name for the marker (for matching across GWAS output files and marker LD file)
- P.value: P-value from the GWAS analysis, used to make marker-selections
- effect: the effect ('beta') for the marker. computed in a GWAS as regression on 0,1,2 genotype values

The LD-information needed is a matrix with the variance-covariance between all the markers, where the original marker data was coded as 0,1,2. To allow for missing genotypes, this covariance-matrix is computed as:

```R
LDmatrix = var(geno,na.rm=TRUE,use="pair")
```

where 'geno' is a matrix (or data frame?) with 0,1,2 coded genotypes. 

### mtmas() function

Here needs a list and explanation of all arguments for the mtmas() function.

### Example

An example is provided from a 5-trait mapping study.
