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
