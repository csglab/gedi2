# GEDI-2
<img src="/Logo.svg" align="right" height="240"/>

[![R-CMD-check](https://github.com/csglab/gedi2/workflows/R-CMD-check/badge.svg)](https://github.com/csglab/gedi2/actions)
[![test-coverage](https://github.com/csglab/gedi2/workflows/test-coverage/badge.svg)](https://github.com/csglab/gedi2/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/gedi)](https://CRAN.R-project.org/package=gedi)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/MIT)
[![Documentation](https://img.shields.io/badge/Docs-Learn%20More-blue.svg)](https://github.com/csglab/gedi2/wiki)

A high-performance, memory-efficient R package for integrating gene expression data from single-cell RNA sequencing experiments. GEDI 2.0 implements a unified generative model for interpretable latent embedding of multi-sample, multi-condition single-cell data.

See the **full Documentation** in the [wiki page.](https://github.com/csglab/gedi2/wiki) <br/>

**Python implementation** of GEDI 2 is available at:  
  https://github.com/csglab/gedi2py

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


### Reproducible Code
All reproducible code, scripts, and resources used in this project are available at: https://github.com/csglab/gedi2_manuscript


### Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

### Reporting Issues
Report bugs and request features at: [https://github.com/csglab/gedi2/issues](https://github.com/csglab/gedi2/issues)
## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


