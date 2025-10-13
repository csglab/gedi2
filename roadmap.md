# GEDI Multi-Language Package Development Progress

**Project:** <br/>
GEDI v2 - Gene Expression Data Integration  <br/>
**Goal:** <br/>
Create memory-efficient R and Python packages wrapping C++ core  <br/>
**Status:** <br/>
**R Package** <br/>
IN PROGRESS <br/>  
**Python Package** <br/>
IN PROGRESS  <br/>

**Last Updated:** 2025-10-12

---

## Project Overview

### Core Architecture (Completed)
- [x] **C++ Core:** `gedi_v2_standalone.cpp` - Stateful GEDI class with full optimization
- [x] **C++ Utilities:** `matrix_operations.cpp` - Sparse/dense matrix operations
- [x] **Design:** Single source of truth in C++, thin wrappers in R/Python

### Key Design Decisions
1. **Memory Strategy:** All data lives in C++ backend (~10 GB), wrapper objects are tiny (~1 KB)
2. **API Pattern:** `setup()` -> `initialize_lvs()` -> `optimize()` or `train()` (all-in-one)
3. **Results Caching:** Wrappers cache last result from C++ in `.lastResult` / `_last_result`
4. **No Duplication:** Data created in R/Python, moved to C++, then deleted from wrapper language

---

## Phase 1: R Package (COMPLETE)

### 1.1 Core Implementation
- [x] **File:** `gedi_v2_r6_class.R`
- [x] **Class:** R6 class with proper encapsulation
- [x] **Methods:**
  - [x] `CreateGEDIObject()` - Factory function
  - [x] `$setup()` - Data preparation and C++ object creation
  - [x] `$initialize_lvs()` - Initialize latent variables via rSVD
  - [x] `$optimize()` - Run block coordinate descent
  - [x] `$train()` - Combined initialize + optimize
  - [x] `$print()` - Display model summary
- [x] **Active Bindings:** `$params`, `$aux`, `$hyperparams`, `$target`, `$tracking`, `$Z`, `$geneIDs`, `$cellIDs`, `$sampleNames`

### 1.2 Data Handling
- [x] Sparse matrix support (preserves sparsity during setup)
- [x] Converts to dense for C++ (required by Eigen::MatrixXd)
- [x] Handles three input types:
  - [x] Single M matrix (counts)
  - [x] Paired M matrices (CITE-seq style)
  - [x] Binary X matrix (probit model)
  - [x] Pre-processed Y matrix

### 1.3 Helper Functions
- [x] `newLogger()` - Verbosity control
- [x] `generateRsvdO()` - Random matrix for rSVD initialization
- [x] C++ functions: `compute_s_0()`, `compute_Mp()`, `compute_o_0()` (sparse)
- [x] C++ functions: `compute_s_0_dense()`, `compute_Yp_dense()`, `compute_o_0_dense()` (dense)
- [x] C++ functions: `VecVecProduct()`, `MatVecProduct()` (general utilities)

### 1.4 Issues Resolved
- [x] Fixed: R6 `initialize()` name conflict (renamed to `initialize_lvs()`)
- [x] Fixed: "Not a matrix" error (convert sparse to dense before slicing by sample)
- [x] Fixed: Memory duplication (data deleted after C++ object creation)
- [x] Fixed: Logger initialization moved to `$setup()` from constructor

### 1.5 Testing
- [x] Tested with real data: 20,000 cells, sparse matrix
- [x] Confirmed memory efficiency: R object size ~1 KB
- [x] Confirmed methods work: `$initialize_lvs()` -> `$optimize()` -> access `$Z`
- [x] Confirmed `$train()` one-shot method works

### 1.6 R Package Usage Pattern
```r
# Load
sourceCpp("matrix_operations.cpp")
sourceCpp("gedi_v2_standalone.cpp")
source("gedi_v2_r6_class.R")

# Create model
set.seed(42)
model <- CreateGEDIObject(
  Samples = sample_labels,
  M = count_matrix,
  K = 10,
  verbose = 1,
  num_threads = 4
)

# Option 1: Step-by-step
model$initialize_lvs()
model$optimize(iterations = 50, track_interval = 5)

# Option 2: One-shot
model$train(iterations = 50, track_interval = 5)
```

---

## Phase 2: Python Package (IN PROGRESS)

### 2.1 Project Structure (To Do)
```
python-package/
├── gedi/
│   ├── __init__.py          [ ] Not started
│   ├── setup.py             [ ] Not started - Data preparation functions
│   ├── core.py              [ ] Not started - Main GEDI class
│   └── utils.py             [ ] Not started - Helper functions
├── src/
│   └── pybind11_wrapper.cpp [ ] Not started - C++ bindings
├── setup.py                 [ ] Not started - Package installation
├── pyproject.toml           [ ] Not started - Build configuration
└── README.md                [ ] Not started
```

### 2.2 Core Components (To Do)

#### A. Data Preparation (`gedi/setup.py`)
- [ ] `prepare_gedi_data()` - Main preparation function
- [ ] Sparse matrix handling (scipy.sparse)
- [ ] Three input types: M (single/paired), Y, X
- [ ] Hyperparameter initialization
- [ ] Sample organization
- [ ] Prior matrix processing (C, H with SVD)
- [ ] Random matrix generation for rSVD

#### B. pybind11 Wrapper (`src/pybind11_wrapper.cpp`)
- [ ] Module definition: `PYBIND11_MODULE(gedi_cpp, m)`
- [ ] Class binding: `py::class_<GEDI>(m, "GEDICore")`
- [ ] Constructor binding with dict arguments
- [ ] Method bindings:
  - [ ] `initialize(multimodal=False)`
  - [ ] `optimize(iterations, track_interval=5)`
  - [ ] `train(iterations, track_interval=5, multimodal=False)`
- [ ] Automatic type conversion setup (NumPy <-> Eigen)

#### C. Python Class Wrapper (`gedi/core.py`)
- [ ] `class GEDI` definition
- [ ] `__init__()` method
  - [ ] Call `prepare_gedi_data()`
  - [ ] Store metadata (gene_ids, cell_ids, sample_names)
  - [ ] Create `GEDICore` C++ object
  - [ ] Delete Python-side data
- [ ] `initialize_lvs(multimodal=False)` method
- [ ] `optimize(iterations, track_interval)` method
- [ ] `train(iterations, track_interval, multimodal)` method
- [ ] Properties (using `@property` decorator):
  - [ ] `params_`
  - [ ] `aux_`
  - [ ] `hyperparams_`
  - [ ] `target_`
  - [ ] `tracking_`
  - [ ] `Z_`
  - [ ] `embedding_`
  - [ ] `gene_ids_`
  - [ ] `cell_ids_`
  - [ ] `sample_names_`
- [ ] `__repr__()` method
- [ ] `__sizeof__()` method (should return ~1 KB)

#### D. Utility Functions (`gedi/utils.py`)
- [ ] Logger class (similar to R's `newLogger`)
- [ ] Random matrix generator
- [ ] Helper functions for common operations

#### E. Installation Setup (`setup.py`)
- [ ] Package metadata
- [ ] Extension module configuration
- [ ] Include directories (pybind11, Eigen3)
- [ ] Compilation flags (C++14, OpenMP, optimization)
- [ ] Dependencies: numpy, scipy, pybind11

### 2.3 Installation Requirements
- [ ] pybind11 >= 2.6
- [ ] numpy >= 1.19
- [ ] scipy >= 1.5
- [ ] Eigen3 (system library)
- [ ] C++14 compiler
- [ ] OpenMP support (optional, for parallelization)

### 2.4 Python Implementation Notes

**Key Differences from R:**
1. **Type hints:** Use proper Python type annotations
2. **Property naming:** Use trailing underscore (scikit-learn convention): `Z_` not `Z`
3. **Method naming:** Use snake_case: `initialize_lvs()` not `initializeLVs()`
4. **Error handling:** Use Python exceptions, not R's `stop()`
5. **Logging:** Use Python logging module or custom logger
6. **Memory management:** pybind11 handles C++ object lifetime automatically

**pybind11 Type Conversions (Automatic):**
- `np.ndarray` <-> `Eigen::MatrixXd` / `Eigen::VectorXd`
- `scipy.sparse.csr_matrix` <-> `Eigen::SparseMatrix<double>`
- `dict` <-> `std::map` or C++ struct
- `list` <-> `std::vector`

### 2.5 Python Package Usage Pattern (Target)
```python
import numpy as np
from scipy import sparse
from gedi import GEDI

# Create model
model = GEDI(
    samples=sample_labels,
    M=count_matrix,
    K=10,
    verbose=1,
    num_threads=4
)

# Check memory
import sys
sys.getsizeof(model)  # Should be ~1 KB

# Option 1: Step-by-step
model.initialize_lvs()
model.optimize(iterations=50, track_interval=5)

# Option 2: One-shot
model.train(iterations=50, track_interval=5)

# Access results
Z = model.Z_
params = model.params_
```

---

## Phase 3: Cross-Language Validation (NOT STARTED)

### 3.1 Validation Tests
- [ ] Create shared test dataset
- [ ] Run same data through R and Python
- [ ] Compare results numerically:
  - [ ] Z matrices (allow for sign flips)
  - [ ] params$sigma2 / params_['sigma2']
  - [ ] Bi matrices
  - [ ] Convergence tracking
- [ ] Assert equivalence within floating-point precision

### 3.2 Performance Benchmarks
- [ ] Memory usage comparison
- [ ] Speed comparison (single-threaded)
- [ ] Speed comparison (multi-threaded)
- [ ] Comparison with GEDI v1

---

## Phase 4: Documentation & Packaging (NOT STARTED)

### 4.1 R Package
- [ ] Write vignettes
- [ ] Add roxygen2 documentation
- [ ] Create DESCRIPTION file
- [ ] Add unit tests (testthat)
- [ ] Prepare for CRAN submission
- [ ] Create pkgdown website

### 4.2 Python Package
- [ ] Write Sphinx documentation
- [ ] Add docstrings (NumPy style)
- [ ] Create example notebooks
- [ ] Add unit tests (pytest)
- [ ] Upload to PyPI
- [ ] Create Read the Docs site

### 4.3 Integration Examples
- [ ] R: Seurat integration
- [ ] R: SingleCellExperiment integration
- [ ] Python: scanpy integration
- [ ] Python: anndata integration

---

## Known Issues & Limitations

### Current Limitations
1. **Dense conversion required:** Sparse M matrices converted to dense for C++ (memory overhead during setup)
2. **No streaming:** All data must fit in memory
3. **Limited error messages:** C++ errors not always informative in R/Python

### Potential Improvements (Future)
- [ ] Keep Yi sparse in C++ (requires Eigen::SparseMatrix support)
- [ ] Add progress bars for long optimizations
- [ ] Add early stopping based on convergence criteria
- [ ] Add model serialization (save/load fitted models)
- [ ] Add GPU support for large datasets

---

## Development Environment

### Required Tools
- **R:** >= 4.0, Rcpp, RcppEigen, Matrix, R6
- **Python:** >= 3.7, pybind11, numpy, scipy
- **C++:** C++14 compiler, Eigen3, OpenMP (optional)
- **Build:** CMake (for Python), Rcpp (for R)

### Compilation Commands

**R Package:**
```r
Rcpp::sourceCpp("matrix_operations.cpp")
Rcpp::sourceCpp("gedi_v2_standalone.cpp")
source("gedi_v2_r6_class.R")
```

**Python Package (when ready):**
```bash
python -m pip install pybind11
python setup.py build_ext --inplace
python -m pip install -e .
```

---

## Next Steps (Priority Order)

1. [ ] **IMMEDIATE:** Create `gedi/setup.py` with data preparation functions
2. [ ] **IMMEDIATE:** Create `src/pybind11_wrapper.cpp` with basic bindings
3. [ ] **IMMEDIATE:** Create `gedi/core.py` with GEDI class
4. [ ] **HIGH:** Test basic Python workflow (create -> initialize -> optimize)
5. [ ] **HIGH:** Cross-validate R vs Python results
6. [ ] **MEDIUM:** Add comprehensive error handling
7. [ ] **MEDIUM:** Write documentation and examples
8. [ ] **LOW:** Package for PyPI/CRAN

---

## Questions & Decisions Log

**Q:** Why convert sparse to dense for C++?  
**A:** Eigen::MatrixXd (dense) required for current C++ implementation. Future: use Eigen::SparseMatrix.

**Q:** Why cache results in wrapper instead of calling getters?  
**A:** C++ already returns everything after each operation. Caching avoids repeated calls.

**Q:** Why `initialize_lvs()` instead of `initialize()`?  
**A:** R6 reserves `initialize()` for constructor. Python doesn't have this issue but we keep names consistent.

**Q:** Should wrappers expose all C++ parameters?  
**A:** Yes, for maximum flexibility. Use sensible defaults.

---

## Contact & Resources

- **C++ Core Author:** [Your name]
- **Repository:** [Add when available]
- **Documentation:** [Add when available]
- **Issues:** [Add when available]



3. Paste and save

This gives you a complete progress tracker with checkboxes you can mark as tasks complete!
