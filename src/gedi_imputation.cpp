// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>

using namespace Rcpp;
using namespace Eigen;

// ============================================================================
// Imputation Helper Functions
// ============================================================================

//' Compute Residual After Removing Sample Effects (Internal)
//'
//' Removes sample-specific effects (QiDBi, si, oi) and global offset (o) from Yi.
//' Used to extract the shared biological signal (ZDBi component).
//'
//' @param Yi Dense matrix (J x Ni) - fitted log-expression for sample i
//' @param QiDBi Dense matrix (J x Ni) - sample-specific metagene projections
//' @param si Vector (Ni) - cell-specific library size offsets for sample i
//' @param o Vector (J) - global gene-specific offsets
//' @param oi Vector (J) - sample-specific gene offsets
//'
//' @return Dense matrix (J x Ni) - residual Yi after removing sample effects
//'
//' @details
//' Computes: Yi - QiDBi - (si (x) 1^T) - (o + oi) (x) 1^T
//' This leaves only the ZDBi component (shared biological signal).
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP Yi_resZ(
    Eigen::MatrixXd Yi,
    Eigen::MatrixXd QiDBi,
    Eigen::VectorXd si,
    Eigen::VectorXd o,
    Eigen::VectorXd oi) {
  
  // Yi - QiDBi
  Eigen::MatrixXd result = Yi - QiDBi;
  
  // Subtract si from each row: rowwise() - si.transpose()
  result.rowwise() -= si.transpose();
  
  // Subtract (o + oi) from each column: colwise() - (o + oi)
  Eigen::VectorXd gene_offset = o + oi;
  result.colwise() -= gene_offset;
  
  return Rcpp::wrap(result);
}

//' Predict Yi from Model Parameters (Internal)
//'
//' Reconstructs the model's prediction of Yi (log-expression) from all components.
//'
//' @param ZDBi Dense matrix (J x Ni) - shared metagene projections
//' @param QiDBi Dense matrix (J x Ni) - sample-specific metagene projections
//' @param si Vector (Ni) - cell-specific library size offsets
//' @param o Vector (J) - global gene-specific offsets
//' @param oi Vector (J) - sample-specific gene offsets
//'
//' @return Dense matrix (J x Ni) - predicted Yi
//'
//' @details
//' Computes: Y_hati = ZDBi + QiDBi + (si (x) 1^T) + (o + oi) (x) 1^T
//' This is the full model prediction for log-expression.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP predict_Yhat(
    Eigen::MatrixXd ZDBi,
    Eigen::MatrixXd QiDBi,
    Eigen::VectorXd si,
    Eigen::VectorXd o,
    Eigen::VectorXd oi) {
  
  // ZDBi + QiDBi
  Eigen::MatrixXd result = ZDBi + QiDBi;
  
  // Add si to each row: rowwise() + si.transpose()
  result.rowwise() += si.transpose();
  
  // Add (o + oi) to each column: colwise() + (o + oi)
  Eigen::VectorXd gene_offset = o + oi;
  result.colwise() += gene_offset;
  
  return Rcpp::wrap(result);
}

// ============================================================================
// Variance Functions
// ============================================================================

//' Variance of Yi for Single Count Matrix (Internal)
//'
//' Computes posterior variance of Yi given the fitted values and model variance.
//' For observation type "M" (single count matrix).
//'
//' @param Yi Dense matrix (J x Ni) - fitted log-expression
//' @param sigma2 Scalar - model variance parameter
//'
//' @return Dense matrix (J x Ni) - posterior variance at each position
//'
//' @details
//' Variance formula: Var(Yi | Mi, model) = 1 / (exp(Yi) + 1/sigma2)
//'
//' This comes from the Poisson-lognormal model where:
//' - Mi ~ Poisson(exp(Yi))
//' - Yi ~ N(Y_hati, sigma2)
//'
//' The posterior variance decreases where:
//' - Counts are high (exp(Yi) large)
//' - Model uncertainty is low (sigma2 small)
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP Yi_var_M(
    Eigen::MatrixXd Yi,
    double sigma2) {
  
  // Convert to array for element-wise operations
  Eigen::ArrayXXd Yi_array = Yi.array();
  
  // Variance = 1 / (exp(Yi) + 1/sigma2)
  Eigen::ArrayXXd variance = 1.0 / (Yi_array.exp() + 1.0 / sigma2);
  
  return Rcpp::wrap(variance.matrix());
}

//' Variance of Yi for Paired Count Matrices (Internal)
//'
//' Computes posterior variance of Yi for paired count data (e.g., CITE-seq, ATAC+RNA).
//' For observation type "M_paired".
//'
//' @param Yi Dense matrix (J x Ni) - fitted log-ratio: log((M1+1)/(M2+1))
//' @param M1i Sparse matrix (J x Ni) - first count matrix (e.g., RNA)
//' @param M2i Sparse matrix (J x Ni) - second count matrix (e.g., protein)
//' @param sigma2 Scalar - model variance parameter
//'
//' @return Dense matrix (J x Ni) - posterior variance at each position
//'
//' @details
//' Variance formula: Var(Yi | M1i, M2i, model) = 1 / (M * exp(-|Yi|) / (1+exp(-|Yi|))^2 + 1/sigma2)
//' where M = M1i + M2i
//'
//' This comes from the binomial-logistic-normal model where:
//' - M1i ~ Binomial(M1i + M2i, p)
//' - logit(p) = Yi ~ N(Y_hati, sigma2)
//'
//' The variance depends on:
//' - Total counts M (more counts = less variance)
//' - Yi magnitude (variance minimized at Yi = 0, i.e., p = 0.5)
//' - Model uncertainty sigma2
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP Yi_var_M_paired(
    Eigen::MatrixXd Yi,
    Eigen::SparseMatrix<double> M1i,
    Eigen::SparseMatrix<double> M2i,
    double sigma2) {
  
// Compute exp(-|Yi|) efficiently
  Eigen::ArrayXXd expYi = (-Yi.array().abs()).exp();
  
  // Sparse addition (stays sparse until converted)
  Eigen::SparseMatrix<double> M_sum_sparse = M1i + M2i;
  
  // Convert sparse to dense, then to array (two separate steps)
  Eigen::MatrixXd M_sum_dense = Eigen::MatrixXd(M_sum_sparse);
  Eigen::ArrayXXd M_sum = M_sum_dense.array();
  
  // Compute (1 + exp(-|Yi|))^2
  Eigen::ArrayXXd one_plus_expYi = 1.0 + expYi;
  Eigen::ArrayXXd one_plus_expYi_sq = one_plus_expYi.square();
  
  // Variance denominator: M * exp(-|Yi|) / (1 + exp(-|Yi|))^2 + 1/sigma2
  Eigen::ArrayXXd denominator = M_sum * expYi / one_plus_expYi_sq + (1.0 / sigma2);
  
  // Variance = 1 / denominator
  Eigen::ArrayXXd variance = 1.0 / denominator;
  
  return Rcpp::wrap(variance.matrix());
}

// ============================================================================
// Dispersion Analysis
// ============================================================================

//' Compute Dispersion Statistics (Sparse-Optimized) (Internal)
//'
//' Analyzes the relationship between predicted and observed variance for count data.
//' Only samples nonzero positions from the sparse count matrix for memory efficiency.
//'
//' @param Yi_fitted Dense matrix (J x Ni) - model prediction of log-expression
//' @param Mi Sparse matrix (J x Ni) - observed count matrix
//' @param subsample Integer - maximum number of nonzero positions to sample
//'
//' @return List containing:
//'   - predicted: Vector of predicted counts (lambda = exp(Yi_fitted))
//'   - observed: Vector of observed counts from Mi
//'
//' @details
//' For Poisson model: Mi ~ Poisson(lambda), expected variance = lambda.
//' This function samples nonzero positions and returns predicted vs observed values
//' for downstream dispersion analysis (binning and aggregation done in R).
//'
//' Memory optimization: Works only on sparse nonzero positions (typically 5-10% of matrix).
//' For 30K genes x 10K cells with 5% nonzero: samples from ~15M positions instead of 300M.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List compute_dispersion_sparse(
    Eigen::MatrixXd Yi_fitted,
    Eigen::SparseMatrix<double> Mi,
    int subsample) {
  
  // Extract all nonzero positions from sparse Mi
  std::vector<int> nz_rows, nz_cols;
  std::vector<double> observed_vals;
  
  // Reserve space for efficiency
  nz_rows.reserve(Mi.nonZeros());
  nz_cols.reserve(Mi.nonZeros());
  observed_vals.reserve(Mi.nonZeros());
  
  // Iterate through nonzero elements (sparse iteration - very fast)
  for (int k = 0; k < Mi.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Mi, k); it; ++it) {
      nz_rows.push_back(it.row());
      nz_cols.push_back(it.col());
      observed_vals.push_back(it.value());
    }
  }
  
  // Subsample if we have more nonzeros than requested
  int n_total = nz_rows.size();
  int n_sample = std::min(subsample, n_total);
  
  // Create random sample indices
  std::vector<int> sample_idx(n_total);
  std::iota(sample_idx.begin(), sample_idx.end(), 0);
  
  // Shuffle and take first n_sample elements
  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(sample_idx.begin(), sample_idx.end(), gen);
  sample_idx.resize(n_sample);
  
  // Extract predicted and observed values at sampled positions
  Eigen::VectorXd predicted(n_sample);
  Eigen::VectorXd observed(n_sample);
  
  for (int i = 0; i < n_sample; ++i) {
    int idx = sample_idx[i];
    int row = nz_rows[idx];
    int col = nz_cols[idx];
    
    // Predicted count: lambda = exp(Yi_fitted)
    predicted(i) = std::exp(Yi_fitted(row, col));
    
    // Observed count from Mi
    observed(i) = observed_vals[idx];
  }
  
  return List::create(
    Named("predicted") = predicted,
    Named("observed") = observed,
    Named("n_total_nonzero") = n_total,
    Named("n_sampled") = n_sample
  );
}

// ============================================================================
// Additional Utility (for reference - used in sigma2 calculation)
// ============================================================================

//' Sum of Squared Errors for Paired Data (Internal)
//'
//' Computes the SSE term used in sigma2 optimization for M_paired observation type.
//' This is included for completeness but is called during optimization, not imputation.
//'
//' @param Yi Dense matrix (J x Ni) - fitted log-ratio
//' @param M1i Sparse matrix (J x Ni) - first count matrix
//' @param M2i Sparse matrix (J x Ni) - second count matrix
//' @param ZDBi Dense matrix (J x Ni) - shared metagene projections
//' @param QiDBi Dense matrix (J x Ni) - sample-specific metagene projections
//' @param si Vector (Ni) - cell-specific offsets
//' @param o Vector (J) - global gene offsets
//' @param oi Vector (J) - sample-specific gene offsets
//' @param sigma2 Scalar - current variance estimate
//'
//' @return Scalar - SSE contribution from sample i
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double Yi_SSE_M_paired(
    Eigen::MatrixXd Yi,
    Eigen::SparseMatrix<double> M1i,
    Eigen::SparseMatrix<double> M2i,
    Eigen::MatrixXd ZDBi,
    Eigen::MatrixXd QiDBi,
    Eigen::VectorXd si,
    Eigen::VectorXd o,
    Eigen::VectorXd oi,
    double sigma2) {
  
  // Compute residual: Yi - Y_hati
  Eigen::MatrixXd residual = Yi - ZDBi - QiDBi;
  residual.rowwise() -= si.transpose();
  residual.colwise() -= (o + oi);
  
  // SSE from residual
  double sse = residual.squaredNorm();
  
  // Add contribution from binomial variance term
  Eigen::ArrayXXd expYi = (-Yi.array().abs()).exp();
  Eigen::SparseMatrix<double> M_sum = M1i + M2i;
  Eigen::MatrixXd M_dense = Eigen::MatrixXd(M_sum);  // Sparse to Dense conversion
  Eigen::ArrayXXd M = M_dense.array();               // Dense to Array conversion
  Eigen::ArrayXXd one_plus_expYi = 1.0 + expYi;
  
  // Sum of 1 / (M * exp(-|Yi|) / (1+exp(-|Yi|))^2 + 1/sigma2)
  Eigen::ArrayXXd variance_term = M * expYi / one_plus_expYi.square() + (1.0 / sigma2);
  sse += (1.0 / variance_term).sum();
  
  return sse;
}



