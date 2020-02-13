# R Package - Capybara
Capybara is a tool to measure cell identity and fate transitions. This approach is designed to measure cell identity as a continuum, at a single-cell resolution. Capybara enables classification of discrete entities as well as cells with multiple identities. This package have a dependency on R version (R >= 3.5.0). For details regarding the methods, usage and application, please refer to the following papaer - *\<Fill this in\>*.

## Installation
### Dependencies
Most dependency package can be installed along with Capybara through CRAN. The following dependency may need to be installed manually through BioConductor (Instructions can also be found here: https://bioconductor.org/).

Install BiocManager
```r
install.packages("BiocManager")
```
Install dependency packages
```r
BiocManager::install("limma")
```

### Install the package
Install devtools
```r
install.packages("devtools")
```
Install the package from GitHub.
```r
library("devtools")
devtools::install_github("morris-lab/Capybara")
```
Load the package
```r
library("Capybara")
```



