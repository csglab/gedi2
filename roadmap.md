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

## Phase 2: Auxiliary Functions and Advanced Features (IN PROGRESS)

### Overview

This phase focuses on implementing auxiliary functions and advanced analytical capabilities from the legacy GEDI package. These functions extend the core GEDI v2 functionality to provide comprehensive tools for downstream analysis, visualization, and biological interpretation of single-cell integration results.

### 2.1 Data Retrieval and Basic Accessors

#### 2.1.1 colData.gedi
**Purpose:** Extract cell-level metadata and sample information from GEDI model.

**Implementation Priority:** MUST HAVE

**Functionality:**
- Returns a data frame containing cell identifiers, sample assignments, and associated metadata
- Provides standardized interface compatible with SingleCellExperiment and Seurat objects
- Enables integration with downstream analysis pipelines

**Technical Specifications:**
```r
colData.gedi(object, ...)
```

**Return Value:** data.frame with columns:
- `cell_id`: Unique cell identifier
- `sample`: Sample of origin
- Additional metadata columns if provided during model creation

**Implementation Notes:**
- Should leverage existing `$cellIDs` and `$sampleNames` active bindings
- Maintain compatibility with Bioconductor's colData generic function
- Consider memory efficiency for large datasets (millions of cells)

**Usability Assessment:**
- Essential for downstream workflows
- Standard interface expected by users familiar with scRNA-seq tools
- Low complexity implementation

---

#### 2.1.2 getY.gedi
**Purpose:** Extract normalized/transformed gene expression matrix from GEDI model.

**Implementation Priority:** MUST HAVE

**Functionality:**
- Returns the Y matrix (log-transformed or normalized expression)
- Provides access to the data used for model fitting
- Supports both dense and sparse matrix output formats

**Technical Specifications:**
```r
getY.gedi(object, sparse = TRUE, ...)
```

**Parameters:**
- `object`: GEDI model object
- `sparse`: Logical, whether to return sparse matrix format (default: TRUE)

**Return Value:** Matrix (dgCMatrix or matrix) with dimensions genes x cells

**Implementation Notes:**
- Should check if Y exists in target list
- For count-based models, reconstruct Y from M if not stored
- Implement sparse conversion for memory efficiency
- Consider lazy evaluation for large datasets

**Performance Considerations:**
- Sparse format critical for scalability (10x memory reduction typical)
- May require conversion from C++ backend storage format
- Cache result if called repeatedly

---

#### 2.1.3 getY.var.gedi
**Purpose:** Calculate gene-wise variance of normalized expression across all cells or within samples.

**Implementation Priority:** SHOULD HAVE

**Functionality:**
- Computes variance of Y matrix gene-wise
- Supports per-sample or global variance calculation
- Used for feature selection and quality control

**Technical Specifications:**
```r
getY.var.gedi(object, per.sample = FALSE, ...)
```

**Parameters:**
- `object`: GEDI model object
- `per.sample`: Logical, compute variance per sample vs. global

**Return Value:**
- If `per.sample = FALSE`: Numeric vector of length n_genes
- If `per.sample = TRUE`: Matrix (n_genes x n_samples)

**Implementation Notes:**
- Leverage efficient matrix operations for variance calculation
- For sparse Y, use specialized sparse variance algorithms
- Consider Bessel's correction (n-1 denominator)

**Use Cases:**
- Identifying highly variable genes
- Quality control and outlier detection
- Feature selection for downstream analysis

---

### 2.2 Model Component Extractors

#### 2.2.1 getDB.gedi
**Purpose:** Extract the product of scaling matrix D and cell embeddings B (D %*% B).

**Implementation Priority:** MUST HAVE

**Functionality:**
- Returns scaled cell embeddings D × B for all samples
- Fundamental component for computing reconstructed expression
- Used in visualization and downstream analysis

**Technical Specifications:**
```r
getDB.gedi(object, sample = NULL, ...)
```

**Parameters:**
- `object`: GEDI model object
- `sample`: Optional, specific sample name or index (returns all if NULL)

**Return Value:**
- Matrix (K x N_total) if sample = NULL
- Matrix (K x N_i) for specific sample i

**Implementation Notes:**
- D is K x K diagonal matrix (stored as vector)
- B is list of K x N_i matrices (one per sample)
- Efficient computation: diag(D) * B for each sample
- Consider caching if called repeatedly

**Performance Optimization:**
- Diagonal scaling is O(K * N) operation
- Parallelize across samples if multiple samples requested
- Return sparse format if B is sparse

---

#### 2.2.2 getZDB.gedi
**Purpose:** Compute full embedding Z × D × B (shared metagenes × scaled cell embeddings).

**Implementation Priority:** MUST HAVE

**Functionality:**
- Returns the complete low-dimensional representation of cells
- Primary output for visualization (UMAP, t-SNE input)
- Captures both shared biology (Z) and cell-specific variation (DB)

**Technical Specifications:**
```r
getZDB.gedi(object, sample = NULL, ...)
```

**Parameters:**
- `object`: GEDI model object
- `sample`: Optional, specific sample name or index

**Return Value:** Matrix (n_genes x N) representing integrated cell embeddings

**Implementation Notes:**
- Sequential matrix multiplication: (Z %*% D) %*% B
- Optimize order: compute Z %*% D once (J x K), then multiply by each B
- For large J, consider dimension reduction first

**Use Cases:**
- Primary input for UMAP/t-SNE visualization
- Batch-corrected expression for downstream clustering
- Cross-sample cell type identification

**Performance Considerations:**
- Most computationally expensive accessor function
- Cache result if used multiple times
- Consider chunked computation for very large datasets

---

#### 2.2.3 getADB.gedi
**Purpose:** Extract sample-specific adjusted embeddings (Z + Q_i) × D × B_i.

**Implementation Priority:** SHOULD HAVE

**Functionality:**
- Returns sample-specific embedding including sample variation Q_i
- Captures batch effects and sample-specific biology
- Used for understanding sample heterogeneity

**Technical Specifications:**
```r
getADB.gedi(object, sample, ...)
```

**Parameters:**
- `object`: GEDI model object
- `sample`: Required, sample name or index

**Return Value:** Matrix (n_genes x N_i) for sample i

**Implementation Notes:**
- A_i = Z + Q_i (sample-specific metagenes)
- Returns A_i × D × B_i
- Q may be NULL if not using sample-specific factors

**Use Cases:**
- Sample-level differential expression
- Identifying sample-specific cell states
- Batch effect visualization

---

#### 2.2.4 getA.gedi
**Purpose:** Extract sample-specific metagenes A_i = Z + Q_i without cell embeddings.

**Implementation Priority:** NICE TO HAVE

**Functionality:**
- Returns metagene matrix adjusted for sample-specific effects
- Useful for understanding sample-level biology
- Smaller memory footprint than full embedding

**Technical Specifications:**
```r
getA.gedi(object, sample, ...)
```

**Parameters:**
- `object`: GEDI model object
- `sample`: Required, sample name or index

**Return Value:** Matrix (n_genes x K) for sample i

**Implementation Notes:**
- Simple addition: Z + Q_i
- Q_i may be zero matrix if not estimated
- Returns reference to Z if Q not used in model

---

#### 2.2.5 svd.gedi
**Purpose:** Perform Singular Value Decomposition on model components for dimensionality reduction visualization.

**Implementation Priority:** SHOULD HAVE (but requires refactoring)

**Functionality:**
- Computes SVD of combined embedding for 2D/3D visualization
- Reduces high-dimensional GEDI output to interpretable dimensions
- Returns both cell coordinates and gene loadings

**Technical Specifications:**
```r
svd.gedi(object, rank = 2, ...)
```

**Parameters:**
- `object`: GEDI model object
- `rank`: Number of SVD components to compute (default: 2 for visualization)

**Return Value:** List containing:
- `u`: Left singular vectors (gene loadings, J x rank)
- `d`: Singular values (rank-length vector)
- `v`: Right singular vectors (cell coordinates, N x rank)

**Implementation Notes - Legacy Issues:**
The legacy implementation has significant performance issues:

**Problem 1: Redundant Computations**
```r
# Legacy code performs SVD 3 separate times:
svd1 <- svd(Z %*% D, nu = rank, nv = 0)          # SVD on ZD
svd2 <- svd(t(do.call(cbind, Bi)), nu = rank, nv = 0)  # SVD on B
svd3 <- svd(svd1$u %*% diag(svd1$d) %*% t(svd2$u))    # Final SVD
```
**Solution:** Single SVD on ZDB is mathematically equivalent and far more efficient:
```r
ZDB <- getZDB.gedi(object)
svd_result <- svd(ZDB, nu = rank, nv = rank)
```

**Problem 2: Memory Inefficiency**
- Legacy code materializes full ZDB matrix multiple times
- Unnecessary intermediate matrix allocations

**Recommended Refactor:**
```r
svd.gedi <- function(object, rank = 2, method = c("full", "randomized")) {
  method <- match.arg(method)

  # Compute ZDB once
  ZDB <- getZDB.gedi(object)

  # Use randomized SVD for large matrices
  if (method == "randomized" || ncol(ZDB) > 10000) {
    result <- rsvd::rsvd(ZDB, k = rank)
  } else {
    result <- svd(ZDB, nu = rank, nv = rank)
  }

  return(list(
    u = result$u,     # Gene loadings (J x rank)
    d = result$d[1:rank],  # Singular values
    v = result$v      # Cell coordinates (N x rank)
  ))
}
```

**Use Cases:**
- 2D scatter plots of cells (v[,1] vs v[,2])
- Gene contribution analysis (u matrix)
- Variance explained by each component (d^2)

---

### 2.3 Differential Analysis Functions

#### 2.3.1 getDiffQ.gedi
**Purpose:** Identify genes with differential sample-specific metagene contributions.

**Implementation Priority:** MUST HAVE

**Functionality:**
- Tests for differences in Q_i matrices across samples
- Identifies genes with sample-specific expression patterns
- Accounts for batch effects vs. biological differences

**Technical Specifications:**
```r
getDiffQ.gedi(object, sample1, sample2, genes = NULL, ...)
```

**Parameters:**
- `object`: GEDI model object
- `sample1`, `sample2`: Sample identifiers for comparison
- `genes`: Optional, restrict to specific genes

**Return Value:** data.frame with columns:
- `gene_id`: Gene identifier
- `log_fold_change`: Log2 fold change in Q contribution
- `p_value`: Statistical significance
- `adj_p_value`: FDR-adjusted p-value

**Statistical Method:**
- Compare Q_1[g,] vs Q_2[g,] for each gene g
- Use appropriate test (t-test, Wilcoxon, or permutation)
- Multiple testing correction (Benjamini-Hochberg)

**Implementation Notes:**
- Requires Q matrices estimated during model fitting
- May need to account for gene-level variance (sigma2)
- Consider baseline expression level in interpretation

---

#### 2.3.2 getDiffO.gedi
**Purpose:** Identify genes with differential sample-specific offsets.

**Implementation Priority:** SHOULD HAVE

**Functionality:**
- Tests for differences in o_i (sample-specific offset terms)
- Detects global expression shifts between samples
- Complements getDiffQ by focusing on additive effects

**Technical Specifications:**
```r
getDiffO.gedi(object, sample1, sample2, genes = NULL, ...)
```

**Parameters:**
- Similar to getDiffQ.gedi

**Return Value:** Similar to getDiffQ.gedi but for offset differences

**Implementation Notes:**
- Offsets represent sample-level expression shifts
- May correlate with library size or technical factors
- Should be interpreted alongside Q differences

---

#### 2.3.3 getDiffExp.gedi
**Purpose:** Perform cluster-free differential expression along continuous cell state trajectories.

**Implementation Priority:** MUST HAVE (core GEDI feature)

**Functionality:**
- Identifies genes with expression changes along latent dimensions
- No discrete clustering required
- Tests for continuous expression gradients

**Technical Specifications:**
```r
getDiffExp.gedi(object, dimension = 1, genes = NULL, method = c("correlation", "regression"), ...)
```

**Parameters:**
- `object`: GEDI model object
- `dimension`: Which latent dimension to test (default: 1)
- `genes`: Optional, restrict to specific genes
- `method`: Statistical approach (correlation or regression)

**Return Value:** data.frame with columns:
- `gene_id`: Gene identifier
- `correlation`: Spearman or Pearson correlation with dimension
- `slope`: Regression slope (if method = "regression")
- `p_value`: Statistical significance
- `adj_p_value`: FDR-adjusted p-value
- `direction`: "increasing" or "decreasing"

**Implementation Notes:**
- Project cells onto latent dimension (use B[dimension,])
- Compute correlation between cell coordinates and gene expression
- Rank genes by absolute correlation or regression coefficient

**Use Cases:**
- Trajectory analysis without pseudotime
- Identifying markers of continuous cell state transitions
- Understanding biological meaning of latent dimensions

**Performance Considerations:**
- Can parallelize across genes
- Consider sparse expression matrices for efficiency
- May need to subset cells for very large datasets

---

### 2.4 Advanced Visualization Functions

#### 2.4.1 svd.vectorField.gedi
**Purpose:** Generate vector field representation of cell state dynamics in latent space.

**Implementation Priority:** NICE TO HAVE

**Functionality:**
- Computes directional trends in cell state transitions
- Visualizes potential developmental trajectories
- Based on local density gradients in latent space

**Technical Specifications:**
```r
svd.vectorField.gedi(object, dims = c(1, 2), resolution = 20, ...)
```

**Parameters:**
- `object`: GEDI model object
- `dims`: Which dimensions to visualize (e.g., c(1,2) for first two PCs)
- `resolution`: Grid resolution for vector field

**Return Value:** data.frame with columns:
- `x`, `y`: Grid coordinates
- `dx`, `dy`: Vector field directions

**Implementation Notes - Legacy Issues:**
- Legacy code has significant code duplication with svd.activityGradient.gedi
- Both functions share 80% identical code
- Opportunity for refactoring into shared utilities

**Recommended Refactor:**
Create shared helper functions:
```r
.compute_cell_projections <- function(object, dims) { ... }
.create_grid <- function(coords, resolution) { ... }
.interpolate_field <- function(grid, values) { ... }
```

**Use Cases:**
- Visualizing differentiation trajectories
- Understanding cell state transitions
- Identifying attractor states in development

---

#### 2.4.2 svd.activityGradient.gedi
**Purpose:** Compute activity gradients for gene sets or pathways along latent dimensions.

**Implementation Priority:** NICE TO HAVE

**Functionality:**
- Projects pathway activity onto low-dimensional representation
- Identifies biological processes varying along trajectories
- Complements vector field visualization

**Technical Specifications:**
```r
svd.activityGradient.gedi(object, gene_sets, dims = c(1, 2), ...)
```

**Parameters:**
- `object`: GEDI model object
- `gene_sets`: List of gene sets or pathways
- `dims`: Dimensions for projection

**Return Value:** List containing gradient information for each gene set

**Implementation Notes:**
- Shares substantial code with svd.vectorField.gedi
- Should be refactored to use common utilities
- Consider gene set enrichment score algorithms (ssGSEA, AUCell)

---

#### 2.4.3 svd.joint_vectorField_gradient.gedi
**Purpose:** Combined visualization of vector field and activity gradients.

**Implementation Priority:** CAN SKIP (composite function)

**Functionality:**
- Wrapper combining vectorField and activityGradient
- Creates integrated visualization
- Primarily convenience function

**Implementation Notes:**
- Can be user-level function calling both sub-functions
- Low priority given individual functions provide same capability
- May be better as vignette example than exported function

---

#### 2.4.4 getActivityGradients.gedi
**Purpose:** Compute pathway activity scores for cells along latent dimensions.

**Implementation Priority:** SHOULD HAVE

**Functionality:**
- Calculates gene set enrichment scores per cell
- Projects onto latent dimensions
- Links GEDI embedding to biological pathways

**Technical Specifications:**
```r
getActivityGradients.gedi(object, gene_sets, method = c("mean", "gsva", "aucell"), ...)
```

**Parameters:**
- `object`: GEDI model object
- `gene_sets`: Named list of gene sets
- `method`: Enrichment score calculation method

**Return Value:** Matrix (n_gene_sets x n_cells) of activity scores

**Implementation Notes:**
- Support multiple enrichment algorithms
- Integrate with existing tools (GSVA, AUCell packages)
- Consider computational cost for large gene set databases

**Use Cases:**
- Pathway activity visualization
- Functional annotation of cell states
- Hypothesis generation for biological processes

---

### 2.5 Visualization Helper Functions

#### 2.5.1 plot_vectorField
**Purpose:** Plotting function for vector field visualization.

**Implementation Priority:** SHOULD HAVE

**Functionality:**
- Creates publication-quality vector field plots
- Overlays cell projections with directional arrows
- Customizable aesthetics

**Technical Specifications:**
```r
plot_vectorField(object, vector_field, color_by = NULL, ...)
```

**Parameters:**
- `object`: GEDI model object
- `vector_field`: Output from svd.vectorField.gedi
- `color_by`: Optional, metadata column for cell coloring

**Return Value:** ggplot2 object

**Implementation Notes:**
- Use ggplot2 for modern, customizable graphics
- Support multiple color schemes
- Allow overlay of multiple features (cell density, metadata, gene expression)

**Design Considerations:**
- Separate plotting logic from computation
- Enable faceting by sample or condition
- Export vector formats (PDF, SVG) for publications

---

#### 2.5.2 plot_embedding
**Purpose:** General-purpose embedding visualization function.

**Implementation Priority:** SHOULD HAVE

**Functionality:**
- Scatter plot of cells in 2D latent space
- Color by gene expression, metadata, or clusters
- Supports multiple layout algorithms (SVD, UMAP, t-SNE)

**Technical Specifications:**
```r
plot_embedding(object, dims = c(1, 2), color_by = NULL, embedding = c("svd", "umap", "tsne"), ...)
```

**Parameters:**
- `object`: GEDI model object
- `dims`: Dimensions to plot
- `color_by`: Feature for coloring (gene name, metadata column, etc.)
- `embedding`: Type of embedding (default: "svd")

**Return Value:** ggplot2 object

**Implementation Notes:**
- Integrate with Seurat and SingleCellExperiment visualizations
- Support both continuous and discrete color scales
- Enable interactive plotting (plotly) as option

**Use Cases:**
- Exploratory data analysis
- Quality control visualization
- Figure generation for publications

---

### 2.6 Quality Control Functions

#### 2.6.1 dispersion.gedi
**Purpose:** Calculate and analyze gene dispersion patterns for quality control.

**Implementation Priority:** MUST HAVE

**Functionality:**
- Computes mean-variance relationship across genes
- Identifies highly variable genes
- Detects technical artifacts and outliers

**Technical Specifications:**
```r
dispersion.gedi(object, method = c("deviance", "pearson"), ...)
```

**Parameters:**
- `object`: GEDI model object
- `method`: Dispersion calculation method

**Return Value:** data.frame with columns:
- `gene_id`: Gene identifier
- `mean_expression`: Mean log expression
- `dispersion`: Dispersion estimate
- `residual_dispersion`: Residual after trend fitting

**Implementation Notes:**
- Model mean-variance relationship (e.g., loess regression)
- Standardize dispersion by expected variance
- Flag genes with extreme dispersion values

**Use Cases:**
- Feature selection for downstream analysis
- Technical quality assessment
- Identifying problematic genes or samples

---

#### 2.6.2 plot_dispersion
**Purpose:** Visualization of gene dispersion patterns.

**Implementation Priority:** SHOULD HAVE

**Functionality:**
- Mean-variance scatter plot with fitted trend
- Highlights highly variable genes
- Flags outliers and potential artifacts

**Technical Specifications:**
```r
plot_dispersion(object, highlight_genes = NULL, ...)
```

**Parameters:**
- `object`: GEDI model object
- `highlight_genes`: Optional, specific genes to label

**Return Value:** ggplot2 object

**Implementation Notes:**
- Plot mean expression vs. dispersion
- Overlay fitted trend line
- Color code genes by selection status
- Interactive version for gene identification

---

### 2.7 Implementation Priorities and Timeline

#### Phase 2A: Essential Functions (Weeks 1-3)
**Priority: MUST HAVE**
- [ ] colData.gedi - Cell metadata accessor
- [ ] getY.gedi - Expression matrix accessor
- [ ] getDB.gedi - Scaled embedding extractor
- [ ] getZDB.gedi - Full embedding computation
- [ ] getDiffExp.gedi - Continuous differential expression
- [ ] dispersion.gedi - Quality control metrics

#### Phase 2B: Extended Analysis (Weeks 4-6)
**Priority: SHOULD HAVE**
- [ ] getY.var.gedi - Variance calculation
- [ ] getADB.gedi - Sample-specific embeddings
- [ ] getDiffQ.gedi - Sample-specific differential metagenes
- [ ] getDiffO.gedi - Sample-specific offset differences
- [ ] svd.gedi (refactored) - Efficient dimensionality reduction
- [ ] getActivityGradients.gedi - Pathway activity scores
- [ ] plot_embedding - General visualization
- [ ] plot_vectorField - Trajectory visualization
- [ ] plot_dispersion - QC visualization

#### Phase 2C: Advanced Features (Weeks 7-8)
**Priority: NICE TO HAVE**
- [ ] getA.gedi - Sample metagenes accessor
- [ ] svd.vectorField.gedi - Vector field analysis
- [ ] svd.activityGradient.gedi - Pathway gradients

#### Phase 2D: Optional Extensions (Future)
**Priority: CAN SKIP**
- [ ] svd.joint_vectorField_gradient.gedi - Composite visualization

---

### 2.8 Code Quality and Maintainability Goals

#### Refactoring Priorities
1. **Eliminate Redundancy:**
   - svd.gedi: Replace three SVDs with single efficient computation
   - Visualization functions: Extract shared grid computation and interpolation logic

2. **Dependency Management:**
   - Remove hard-coded library() calls from function bodies
   - Use @importFrom in roxygen2 documentation
   - Declare all dependencies in DESCRIPTION

3. **Performance Optimization:**
   - Implement caching for expensive computations (getZDB, svd results)
   - Use sparse matrix operations where applicable
   - Parallelize gene-wise operations in differential analysis

4. **Documentation Standards:**
   - Comprehensive roxygen2 documentation for all exported functions
   - Usage examples in function documentation
   - Vignettes demonstrating workflows

5. **Testing Strategy:**
   - Unit tests for each accessor function
   - Integration tests for analysis workflows
   - Performance benchmarks for large datasets

---

### 2.9 Integration with Ecosystem

#### Seurat Integration
- [ ] S3 methods for Seurat objects
- [ ] Export GEDI results to Seurat assays and reductions
- [ ] Compatible with Seurat visualization tools

#### SingleCellExperiment Integration
- [ ] S4 methods for SCE objects
- [ ] Store GEDI results in reducedDims and colData
- [ ] Leverage existing Bioconductor infrastructure

#### Pathway Databases
- [ ] Support MSigDB gene sets
- [ ] Integration with clusterProfiler for enrichment
- [ ] Custom gene set input formats

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
