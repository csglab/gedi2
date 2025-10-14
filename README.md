# GEDI v2: Gene Expression Data Integration

[![R-CMD-check](https://github.com/yourusername/gedi/workflows/R-CMD-check/badge.svg)](https://github.com/yourusername/gedi/actions)
[![test-coverage](https://github.com/yourusername/gedi/workflows/test-coverage/badge.svg)](https://github.com/yourusername/gedi/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A high-performance, memory-efficient R package for integrating gene expression data from single-cell RNA sequencing experiments. GEDI v2 implements a unified generative model for interpretable latent embedding of multi-sample, multi-condition single-cell data.

## Overview

GEDI (Gene Expression Decomposition and Integration) enables:

- **Cross-sample integration** on par with state-of-the-art methods
- **Cluster-free differential expression** analysis along the continuum of cell states
- **Pathway and regulatory network** activity inference in single cells
- **Machine learning-based prediction** of sample characteristics from single-cell data

### Key Features

- **Memory-efficient architecture**: All data lives in C++ backend (~10GB data with ~1KB R objects)
- **High performance**: OpenMP parallelization with C++14 optimization
- **Sparse matrix support**: Efficiently handles sparse single-cell data
- **Multiple data modalities**: Count matrices (M), paired data (CITE-seq), binary indicators (X), or pre-processed expression (Y)
- **Flexible modeling**: Dimensionality reduction with batch effect correction

## Installation

### System Requirements

- **R** >= 4.0.0
- **C++ Compiler** with C++14 support
- **Eigen** >= 3.3.0 (linear algebra library)
- **OpenMP** (optional, for parallelization)

### Installing from GitHub

```r
# Install devtools if not already installed
install.packages("devtools")

# Install gedi from GitHub
devtools::install_github("yourusername/gedi")
```

### Installing from Source

```bash
# Download the source package
R CMD INSTALL gedi_2.0.0.tar.gz
```

### System Dependencies

**macOS:**
```bash
brew install eigen
```

**Ubuntu/Debian:**
```bash
sudo apt-get install libeigen3-dev
```

**CentOS/RHEL:**
```bash
sudo yum install eigen3-devel
```

## Quick Start

```r
library(gedi)

# Load example data
data("pbmc_small", package = "Seurat")

# Create GEDI model
model <- CreateGEDIObject(
  Samples = pbmc_small@meta.data$orig.ident,
  M = pbmc_small@assays$RNA@counts,
  K = 10,
  verbose = 1
)

# Train the model
model$train(iterations = 50, track_interval = 5)

# Access latent representation
Z <- model$Z          # Shared metagenes
B <- model$B          # Cell embeddings
params <- model$params  # Model parameters
```

## Usage Examples

### Basic Workflow

```r
# 1. Setup: Create GEDI object from count matrix
model <- CreateGEDIObject(
  Samples = sample_labels,  # Factor indicating sample of origin
  M = count_matrix,         # Raw count matrix (genes x cells)
  K = 15,                   # Number of latent variables
  verbose = 1
)

# 2. Initialize: Initialize latent variables via randomized SVD
model$initialize_lvs()

# 3. Optimize: Fit model via block coordinate descent
model$optimize(iterations = 100)

# 4. Extract: Access integrated embedding and parameters
embedding <- model$Z %*% model$D %*% model$B
```

### Advanced Features

#### Working with Pre-processed Data

```r
# Use log-transformed expression matrix instead of counts
model <- CreateGEDIObject(
  Samples = sample_labels,
  Y = log_expression_matrix,  # Pre-processed expression
  K = 10
)
```

#### CITE-seq Paired Data

```r
# Analyze paired RNA and protein measurements
model <- CreateGEDIObject(
  Samples = sample_labels,
  M = list(RNA_counts, protein_counts),  # Paired count matrices
  K = 10
)
```

#### Incorporating Prior Knowledge

```r
# Include gene-level biological prior
model <- CreateGEDIObject(
  Samples = sample_labels,
  M = count_matrix,
  C = gene_prior_matrix,  # Gene-level prior (e.g., pathways)
  K = 10
)
```

## Architecture

GEDI v2 implements a three-layer architecture:

1. **C++ Core**: Stateful GEDI class with full optimization using Eigen and OpenMP
2. **R6 Wrapper**: Thin R6 class exposing methods and active bindings (~1KB memory)
3. **Factory Function**: `CreateGEDIObject()` for user-friendly model creation

This design enables analysis of datasets with millions of cells while maintaining minimal memory footprint in R.

## Model Components

The GEDI model decomposes gene expression as:

```
Y ~ o + oi + (Z + Qi) × D × Bi + si
```

Where:
- **Y**: Log-transformed expression matrix
- **o**: Global gene-specific offsets
- **oi**: Sample-specific gene offsets
- **Z**: Shared metagenes (latent factors)
- **Qi**: Sample-specific metagene variations
- **D**: Diagonal scaling matrix
- **Bi**: Cell embeddings for sample i
- **si**: Cell-specific library size factors

## Performance

GEDI v2 is optimized for large-scale single-cell datasets:

- **Memory**: ~1KB per R object (data stored in C++)
- **Speed**: Parallelized with OpenMP
- **Scalability**: Handles millions of cells efficiently
- **Sparsity**: Native sparse matrix support

## Documentation

Full documentation is available within R:

```r
# Package documentation
?gedi

# Main function
?CreateGEDIObject

# View methods
?GEDI
```

## Citation

If you use GEDI in your research, please cite:

> Madrigal, A., Lu, T., Soto, L. M., & Najafabadi, H. S. (2024). A unified model for interpretable latent embedding of multi-sample, multi-condition single-cell data. *Nature Communications*, 15(1), 6573.

## Development

GEDI v2 is developed by the Computational and Statistical Genomics Laboratory at McGill University, Montreal.

### Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

### Reporting Issues

Report bugs and request features at: [https://github.com/yourusername/gedi/issues](https://github.com/yourusername/gedi/issues)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Computational and Statistical Genomics Laboratory, McGill University
- Original GEDI implementation: [https://github.com/csglab/GEDI](https://github.com/csglab/GEDI)
- Built with [Rcpp](https://www.rcpp.org/) and [RcppEigen](http://dirk.eddelbuettel.com/code/rcpp.eigen.html)
