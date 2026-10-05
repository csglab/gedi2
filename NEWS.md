# gedi 2.3.6

## Maintenance

* Compatibility with Eigen 5.0 (upcoming RcppEigen 0.4.x). Eigen 5.0 removes
  `Eigen::MappedSparseMatrix`, so the sparse-matrix helpers `compute_Yp()`,
  `compute_Mp()`, `compute_s_0()`, `compute_o_0()` and
  `eigenSparseMatVecProduct()` now take `Eigen::Map<Eigen::SparseMatrix<double>>`
  instead. The package still builds against the current CRAN RcppEigen
  (0.3.4.0.2) and the RcppEigen release candidate (0.4.9.9-2). No user-visible
  changes. Contributed by Dirk Eddelbuettel (@eddelbuettel) in #28; see
  RcppCore/RcppEigen#151 for the wider Eigen 5.0 migration.

## Bug fixes

* Count matrices are now coerced to `dgCMatrix` before reaching the C++
  backend. Matrix >= 1.8 stores integer counts as `igCMatrix`, which caused
  `CreateGEDIObject()` to fail with "Need S4 class dgCMatrix for a mapped
  sparse matrix".

## Documentation

* DESCRIPTION now cites the GEDI 2.0 paper
  (<doi:10.1093/bioinformatics/btag334>), and `LICENSE` is restored to the
  CRAN `YEAR`/`COPYRIGHT HOLDER` stub (full MIT text in `LICENSE.md`).

# gedi 2.3.5

## New features

* `plot_features()` gains two differential projection types (#26):
  * `projection = "diffexp"` plots the per-cell differential expression for
    selected genes, given a `contrast` (optionally adding the global offset via
    `include_O = TRUE`). The full J x N matrix is never materialised.
  * `projection = "diffadb"` plots the per-cell differential pathway activity
    for selected pathways (requires a gene-level prior `C`).
* New method `model$diffADB(contrast)` returns the differential pathway
  activity (num_pathways x N). It is the exact differential of `ADB`
  (i.e. `ADB(Z + dQ) - ADB(Z)`), applying the same `solve_A` shrinkage so it
  stays on the same scale as `model$projections$ADB`.

# gedi 2.3.1

## CRAN Compliance

* Replace all `cat()` / `print()` calls with `message()` or `warning()` for
  suppressible console output.
* Add `verbose` parameters to `seurat_to_gedi()`, `gedi_to_seurat()`,
  `check_optional_dependencies()`, and `install_optional_dependencies()`.
* Replace `cat()`-based progress bars with `txtProgressBar()` across imputation
  and training routines.
* Replace `installed.packages()` with `requireNamespace()` for dependency
  checking.
* Move `hdf5r` from Imports to Suggests (optional dependency for H5AD I/O).
* Add `ggplot2` and `scales` to Imports; add `uwot` and `digest` to Suggests.
* Use CRAN-required two-line LICENSE format.
* Add `@return` documentation tags to all exported and documented functions.
* Replace non-ASCII characters in C++ source files with ASCII equivalents.
* Remove redundant Maintainer field from DESCRIPTION.

# gedi 2.3.0

* Initial public release with C++ backend and R6 interface.
* Support for multiple data modalities (count matrices, paired data, binary
  indicators).
* Latent variable model with block coordinate descent optimization.
* Dimensionality reduction, batch correction, and imputation.
* Differential expression and pathway association analysis.
* H5AD file I/O for Python interoperability.
* Seurat and SingleCellExperiment integration.
