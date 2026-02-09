## R CMD check results

0 errors | 0 warnings | 1 note

## Test environments

* Local Linux (AlmaLinux 9.6, R 4.4.0)
* GitHub Actions (ubuntu-latest, macOS-latest, windows-latest; R-release, R-oldrel)

## Notes

This is the first CRAN submission of the `gedi` package.

### Conflicting package name

The incoming feasibility check flags a conflicting name with Bioconductor's
`GeDi` package. Our package `gedi` (Gene Expression Decomposition and
Integration) is unrelated to `GeDi` (Gene Distance based on Enrichment and
Depletion for Inference). The names differ in capitalisation and the packages
serve entirely different purposes.

### C++14 specification

The package requires C++14 for the RcppEigen linear algebra backend.
OpenMP is optional and detected at compile time via Makevars.
