# GEDI
<img src="/Logo.svg" align="right" height="240"/>

[![R-CMD-check](https://github.com/Arshammik/gedi/workflows/R-CMD-check/badge.svg)](https://github.com/Arshammik/gedi/actions)
[![test-coverage](https://github.com/Arshammik/gedi/workflows/test-coverage/badge.svg)](https://github.com/Arshammik/gedi/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A high-performance, memory-efficient R package for integrating gene expression data from single-cell RNA sequencing experiments. GEDI v2 implements a unified generative model for interpretable latent embedding of multi-sample, multi-condition single-cell data.

```r
devtools::install_github("Arshammik/gedi", auth_token = "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")
```

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
devtools::install_github("Arshammik/gedi")
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
### Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

### Reporting Issues

Report bugs and request features at: [https://github.com/Arshammik/gedi/issues](https://github.com/Arshammik/gedi/issues)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


