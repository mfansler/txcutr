
<!-- README.md is generated from README.Rmd. Please edit that file -->

# txcutr

<!-- badges: start -->

[![R build
status](https://github.com/mfansler/txcutr/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/mfansler/txcutr/actions)
[![Anaconda-Server
Badge](https://anaconda.org/merv/r-txcutr/badges/installer/conda.svg)](https://conda.anaconda.org/merv)
[![Anaconda-Server
Badge](https://anaconda.org/merv/r-txcutr/badges/version.svg)](https://anaconda.org/merv/r-txcutr)
<!-- badges: end -->

## Overview

Various mRNA sequencing library preparation methods generate sequencing
reads from the transcript ends. Quantification of isoform usage can be
improved by using truncated versions of transcriptome annotations when
assigning such reads to isoforms. The `txcutr` package implements some
convenience methods for readily generating such truncated annotations
and their corresponding sequences.

## Installation instructions

**Interim Note:** This package is intended for release on Bioconductor,
but is currently only available through GitHub or Anaconda Cloud.

### Bioconductor

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `txcutr` using from
[Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

## not currently available
#BiocManager::install("txcutr")
```

And the development version from
[GitHub](https://github.com/mfansler/txcutr) with:

``` r
BiocManager::install("mfansler/txcutr")
```

### Conda/Mamba

Users managing R environments with Conda/Mamba can install the package
with:

**Conda**

``` bash
conda install -c conda-forge -c bioconda merv::r-txcutr
```

**Mamba**

``` bash
mamba install -c conda-forge -c bioconda merv::r-txcutr
```

We strongly encourage users to create dedicated R environments. **Do not
install this in your *base* environment!**

## Example

A typical workflow for `txcutr` involves

-   loading an existing annotation as `TxDb` object
-   truncating the annotation
-   exporting the truncated annotation (GTF)
-   exporting supporting files (FASTA, merge TSV)

``` r
library(rtracklayer)
library(txcutr)
library(BSgenome.Hsapiens.UCSC.hg38)

## load human genome
hg38 <- BSgenome.Hsapiens.UCSC.hg38

## load human GENCODE annotation
txdb <- makeTxDbFromGFF("gencode.v38.annotaton.gtf.gz", organism="Homo sapiens")

## truncate to maximum of 500 nts
txdb_w500 <- truncateTxome(txdb, maxTxLength=500)

## export annotation
exportGTF(txdb_w500, file="gencode.v38.txcutr_w500.gtf.gz")

## export FASTA and merge-table
exportFASTA(txdb_w500, genome=BSgenome.Hsapiens.UCSC.hg38, 
            file="gencode.v38.txcutr_w500.fa.gz")
exportMergeTable(txdb_w500, minDistance=200,
                 file="gencode.v38.txcutr_w500.merge.tsv.gz")
```

## Citation

Below is the citation output from using `citation('txcutr')` in R.
Please run this yourself to check for any updates on how to cite
**txcutr**.

``` r
print(citation('txcutr'), bibtex = TRUE)
#> 
#> To cite package 'txcutr' in publications use:
#> 
#>   Mervin Fansler (2021). txcutr: Transcriptome CUTteR. R package
#>   version 0.3.1.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {txcutr: Transcriptome CUTteR},
#>     author = {Mervin Fansler},
#>     year = {2021},
#>     note = {R package version 0.3.1},
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

-   Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.13/BiocCheck)*.
-   Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation website](http://mfansler.github.io/txcutr) is
    automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
-   The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.13/biocthis)*.
