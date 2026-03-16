
<!-- README.md is generated from README.Rmd. Please edit that file -->

# txcutr

<!-- badges: start -->

[![R build
status](https://github.com/mfansler/txcutr/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/mfansler/txcutr/actions)
[![codecov](https://codecov.io/gh/mfansler/txcutr/branch/bioc-check/graph/badge.svg?token=CGGZP68G67)](https://codecov.io/gh/mfansler/txcutr)
[![Anaconda-Server
Badge](https://anaconda.org/bioconda/bioconductor-txcutr/badges/version.svg)](https://anaconda.org/bioconda/bioconductor-txcutr)
<!-- badges: end -->

## Overview

Various mRNA sequencing library preparation methods generate sequencing
reads from the transcript ends. Quantification of isoform usage can be
improved by using truncated versions of transcriptome annotations when
assigning such reads to isoforms. The `txcutr` package implements some
convenience methods for readily generating such truncated annotations
from either their 5’ or 3’ transcript ends and their corresponding
sequences.

## Installation instructions

### Bioconductor

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `txcutr` using from
[Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("txcutr")
```

or the development version with:

``` r
# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("txcutr")
```

### Conda

Users managing R environments with Conda can install the package with:

**Conda**

``` bash
conda install -c conda-forge -c bioconda bioconductor-txcutr
```

We strongly encourage users to create dedicated R environments. **Do not
install this in your *base* environment!**

## Example

A typical workflow for `txcutr` involves

- loading an existing annotation as `TxDb` object
- truncating the annotation from the 3’ or 5’ end
- exporting the truncated annotation (GTF)
- exporting supporting files (FASTA, merge TSV)

``` r
library(rtracklayer)
library(txcutr)
library(BSgenome.Hsapiens.UCSC.hg38)

## load human genome
hg38 <- BSgenome.Hsapiens.UCSC.hg38

## load human GENCODE annotation
txdb <- makeTxDbFromGFF("gencode.v38.annotaton.gtf.gz", organism="Homo sapiens")

## truncate to maximum of 500 nts
txdb_w500 <- truncate3primeTxome(txdb, maxTxLength=500) # use `truncate5primeTxome` for 5' truncation

## export annotation
exportGTF(txdb_w500, file="gencode.v38.txcutr_w500.gtf.gz")

## export FASTA
exportFASTA(txdb_w500, genome=hg38, file="gencode.v38.txcutr_w500.fa.gz")

## export merge-table
exportMergeTable(txdb_w500, minDistance=200,
                 file="gencode.v38.txcutr_w500.merge.tsv.gz")
```

## Citation

Below is the citation output from using `citation('txcutr')` in R.
Please run this yourself to check for any updates on how to cite
**txcutr**.

``` r
print(citation('txcutr'), bibtex = TRUE)
#> To cite package 'txcutr' in publications use:
#> 
#>   Fansler M (2025). _txcutr: Transcriptome CUTteR_.
#>   doi:10.18129/B9.bioc.txcutr
#>   <https://doi.org/10.18129/B9.bioc.txcutr>, R package version 1.16.0,
#>   <https://bioconductor.org/packages/txcutr>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {txcutr: Transcriptome CUTteR},
#>     author = {Mervin Fansler},
#>     year = {2025},
#>     note = {R package version 1.16.0},
#>     url = {https://bioconductor.org/packages/txcutr},
#>     doi = {10.18129/B9.bioc.txcutr},
#>   }
```

Please note that the `txcutr` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `txcutr` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.22/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://mfansler.github.io/txcutr) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.22/biocthis)*.
