## miMediation

miMediation is an R package for performing various mediation tests for microbiome data. It currently implements two methods CAMRA and PhyloMed listed below

- CAMRA reference: Wang Q, Li Y, Peng Y, Tang, ZZ (2026). Error control in microbiome mediator discovery: benchmark and remedy. Submitted.
- PhyloMed reference: Hong Q, Chen G, Tang ZZ (2023). [PhyloMed: a phylogeny-based test of mediation effect in microbiome.](https://doi.org/10.1186/s13059-023-02902-3) Genome Biology, 24(1), 1-21.


See the following for comprehensive and up-to-date documentation:

- The [miMediation R package manual](https://github.com/tangzheng1/miMediation/blob/main/miMediation_1.0.pdf).
- The [tutorial walkthrough of the proposed PhyloMed](https://github.com/tangzheng1/miMediation/blob/main/vignettes/miMediation.pdf).
- The [tutorial walkthrough of the proposed CAMRA](https://github.com/tangzheng1/miMediation/blob/main/vignettes/CAMRA.pdf)

## Author

Yunfei Peng, Qiyu Wang, Qilin Hong

Department of Biostatistics and Medical Informatics, University of Wisconsin-Madison

## Installation

You can install the package from github with:

``` r
devtools::install_github("tangzheng1/miMediation")
```
You can download the [package source](https://github.com/tangzheng1/miMediation/blob/main/miMediation_1.0.tar.gz) and install it manually with:

``` r
install.packages("miMediation_1.0.tar.gz", repos = NULL, type ="source", dependencies = c("Depends", "Imports")) 
```

You can force installation if you already have old version with:

``` r
devtools::install_github("tangzheng1/miMediation", force = TRUE)
```
## Troubleshoot Dependencies

At this point, there may be complaints about missing dependencies. To install missing dependencies on either [CRAN](https://cran.r-project.org/) or [Bioconductor](http://bioconductor.org/install/), start a fresh R session and enter the following:

``` r
# For CRAN
install.packages("missing_package")
# For Biocondocutor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("missing_package") 
# ... and so on
```

## Re-attempt miMediation Installation

Install again, after dependencies have been installed. You should now be done if the installation was successful.

## Load Package, Explore Function Documentation and Vignette

``` r
library(miMediation)
help(package = "miMediation")
?CAMRA
?data.camra
?data.cecal
?data.zeeviD
?phyloMed
?prepareTree
```

If you install the package manually from source, a vignette describing the use of the package is available from within R. Load the package and then use the vignette function.

``` r
vignette("miMediation", package = "miMediation")
```

Otherwise, it will not build vignette by default if you install the package from github because they are time consuming and may require additional packages (here, require `prettydoc` R package, install it before building vignette). You can force building (take ~7 mins) with:

``` r
devtools::install_github("tangzheng1/miMediation", build_vignettes = TRUE)
```
Then, the vignette would be available from within R.

## Getting help

Please use the [issue tracker](https://github.com/tangzheng1/miMediation/issues) to post any bugs, suggestions, or the installation problem.

## News

- Version 0.1 (06/12/2023): Initial version released for publication.
- Version 0.2 (12/27/2023): Minor update. 
  - Implemented simple permutation when no confounding variable is present.
  - Refactored the code for Smith's permutation.
- Version 0.3 (03/25/2024): Minor enhancements. 
  - Address scenarios where the mediator matrix is rank-deficient in certain subcompositions. 
  - Added an intermediate argument outcome_type to handle cases with only two unique values in the continuous outcome.
- Version 1.0 (02/26/2026): Incorporate CAMRA method. 
  
## License

This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3
