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
