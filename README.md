# GEDI-2
<img src="/Logo.svg" align="right" height="240"/>

[![R-CMD-check](https://github.com/Arshammik/gedi/workflows/R-CMD-check/badge.svg)](https://github.com/caglab/gedi2/actions)
[![test-coverage](https://github.com/Arshammik/gedi/workflows/test-coverage/badge.svg)](https://github.com/caglab/gedi2/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/Docs-Learn%20More-blue.svg)](https://github.com/caglab/gedi2/wiki)

A high-performance, memory-efficient R package for integrating gene expression data from single-cell RNA sequencing experiments. GEDI 2.0 implements a unified generative model for interpretable latent embedding of multi-sample, multi-condition single-cell data.


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
devtools::install_github("csglab/gedi2")
```
### Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

### Reporting Issues
Report bugs and request features at: [https://github.com/caglab/gedi2/issues](https://github.com/Arshammik/gedi/issues)
## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


