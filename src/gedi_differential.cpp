// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;
using namespace Eigen;

// ============================================================================
// GEDI Differential Expression Functions
// ============================================================================


/**
 * Compute Differential O (Gene Offset) Effect (Internal C++)
 *
 * Computes the differential effect of sample-level variables on gene-specific
 * offsets: diffO = Ro * H.rotation * contrast.
 * This is the C++ implementation of the old getDiffO.gedi()
 *
 * @param Ro Matrix (J x L) representing effect of sample variables on gene offsets
 * @param H_rotation Rotation matrix (L x num_covariates)
 * @param contrast Vector of length L specifying the contrast
 * @param verbose Integer verbosity level
 *
 * @return Vector of length J representing the differential offset effect
 * for each gene under the specified contrast.
 *
 * @keywords internal
 * @noRd
 */
// [[Rcpp::export]]
Eigen::VectorXd getDiffO_cpp(
    const Eigen::Map<Eigen::MatrixXd>& Ro,
    const Eigen::Map<Eigen::MatrixXd>& H_rotation,
    const Eigen::Map<Eigen::VectorXd>& contrast,
    int verbose = 0
) {

  // Dimension validation
  int J = Ro.rows();
  int L = Ro.cols();
  
  if (H_rotation.rows() != L) {
    stop("Dimension mismatch: H_rotation must have L rows");
  }
  
  if (contrast.size() != L) {
    stop("Dimension mismatch: contrast must have length L");
  }
  
  if (L == 0) {
    stop("Cannot compute diffO: no sample-level prior (H) was provided");
  }
  
  if (verbose >= 1) {
    Rcout << "Computing diffO: " << J << " genes" << std::endl;
  }
  
  // Compute: diffO = Ro * H.rotation * contrast
  VectorXd H_contrast = H_rotation * contrast;
  VectorXd diffO = Ro * H_contrast;
  
  if (verbose >= 1) {
    double mean_val = diffO.mean();
    double std_val = std::sqrt((diffO.array() - mean_val).square().mean());
    Rcout << "diffO computed: mean = " << mean_val 
          << ", std = " << std_val << std::endl;
  }
  
  return diffO;
}


/**
 * Compute Differential Q in Z-space (Internal C++)
 *
 * Computes sample-variable effects on Qi, returning a J x K matrix.
 * This is the C++ implementation of the old getDiffQ.gedi()
 *
 * @param Rk_list List of K matrices (each J x L), representing the effect of
 * sample-level variables on each latent factor k
 * @param H_rotation Rotation matrix (L x num_covariates)
 * @param contrast Vector of length L specifying the contrast
 * @param verbose Integer verbosity level
 *
 * @return Dense matrix diffQ of dimensions J x K representing the predicted
 * differential effect in Z-space.
 *
 * @keywords internal
 * @noRd
 */
// [[Rcpp::export]]
Eigen::MatrixXd getDiffQ_cpp(
    const Rcpp::List& Rk_list,
    const Eigen::Map<Eigen::MatrixXd>& H_rotation,
    const Eigen::Map<Eigen::VectorXd>& contrast,
    int verbose = 0
) {

  // Dimension validation
  int K = Rk_list.size();
  int L = H_rotation.rows();
  int J = -1; // Will be set from first Rk
  
  if (K == 0 || L == 0) {
    stop("Cannot compute diffQ: no sample-level prior (H) was provided (K=0 or L=0)");
  }
  
  // Validate Rk dimensions and get J
  Eigen::MatrixXd Rk_first = as<Eigen::MatrixXd>(Rk_list[0]);
  J = Rk_first.rows();
  if (Rk_first.cols() != L) {
     stop("Dimension mismatch: Rk[1] must have L columns");
  }

  if (contrast.size() != L) {
    stop("Dimension mismatch: contrast must have length L");
  }

  if (verbose >= 1) {
    Rcout << "Computing diffQ (Z-space): " << J << " genes x " << K << " factors" << std::endl;
  }
  
  // Pre-compute H.rotation * contrast
  VectorXd H_contrast = H_rotation * contrast; // L x 1
  
  // Allocate result matrix
  MatrixXd diffQ_Z_space(J, K);
  
  // Compute effect for k=1 (already fetched)
  diffQ_Z_space.col(0) = Rk_first * H_contrast;
  
  // Loop for k = 2 to K
  for (int k = 1; k < K; ++k) {
    Eigen::MatrixXd Rk = as<Eigen::MatrixXd>(Rk_list[k]); // J x L
    
    // Validate dimensions for subsequent Rk
    if (Rk.rows() != J || Rk.cols() != L) {
      stop("Dimension mismatch: all Rk matrices must have the same JxL dimensions");
    }
    
    // Compute effect for this factor: Rk * H_contrast (J x 1)
    diffQ_Z_space.col(k) = Rk * H_contrast;
  }

  if (verbose >= 1) {
    Rcout << "[OK] diffQ (Z-space) computed" << std::endl;
  }
  
  return diffQ_Z_space;
}


/**
 * Compute Differential Expression (Internal C++)
 *
 * Computes the cell-specific differential expression effect (J x N).
 * This is the C++ implementation of the old getDiffExp.gedi()
 *
 * @param Rk_list List of K matrices (each J x L)
 * @param H_rotation Rotation matrix (L x num_covariates)
 * @param contrast Vector of length L specifying the contrast
 * @param D Scaling vector (length K)
 * @param Bi_list List of sample-specific cell projection matrices (K x Ni each)
 * @param verbose Integer verbosity level
 *
 * @return Dense matrix diffExp of dimensions J x N representing the predicted
 * differential expression effect for each gene in each cell.
 *
 * @keywords internal
 * @noRd
 */
// [[Rcpp::export]]
Eigen::MatrixXd getDiffExp_cpp(
    const Rcpp::List& Rk_list,
    const Eigen::Map<Eigen::MatrixXd>& H_rotation,
    const Eigen::Map<Eigen::VectorXd>& contrast,
    const Eigen::Map<Eigen::VectorXd>& D,
    const Rcpp::List& Bi_list,
    int verbose = 0
) {

  // Dimension validation
  int K = Rk_list.size();
  int L = H_rotation.rows();
  int J = -1; // Will be determined from first Rk
  int numSamples = Bi_list.size();
  
  if (contrast.size() != L) {
    stop("Dimension mismatch: contrast must have length L");
  }
  
  if (D.size() != K) {
    stop("Dimension mismatch: D must have length K");
  }
  
  if (K == 0 || L == 0) {
    stop("Cannot compute diffExp: no sample-level prior (H) was provided");
  }
  
  // Validate Rk dimensions and get J
  Eigen::MatrixXd Rk_first = as<Eigen::MatrixXd>(Rk_list[0]);
  J = Rk_first.rows();
  if (Rk_first.cols() != L) {
     stop("Dimension mismatch: Rk[1] must have L columns");
  }
  
  // Count total cells and validate Bi
  int N = 0;
  for (int i = 0; i < numSamples; ++i) {
    Eigen::MatrixXd Bi = as<Eigen::MatrixXd>(Bi_list[i]);
    if (Bi.rows() != K) {
      stop("Dimension mismatch: Bi[%d] must have K rows", i + 1);
    }
    N += Bi.cols();
  }
  
  if (verbose >= 1) {
    Rcout << "Computing diffExp (Cell-space): " << J << " genes x " << N << " cells" << std::endl;
    if (verbose >= 2) {
      Rcout << "  Contrast dimension: " << L << ", Factors: " << K << std::endl;
    }
  }
  
  // === Concatenate B and compute DB ===
  if (verbose >= 2) Rcout << "  Concatenating B and computing DB..." << std::endl;
  MatrixXd B(K, N);
  int col_offset = 0;
  
  for (int i = 0; i < numSamples; ++i) {
    Eigen::MatrixXd Bi = as<Eigen::MatrixXd>(Bi_list[i]);
    int Ni_current = Bi.cols();
    B.block(0, col_offset, K, Ni_current) = Bi;
    col_offset += Ni_current;
  }
  
  MatrixXd DB = D.asDiagonal() * B; // K x N
  
  // === Compute: diffExp = sum_k (Rk * H.rotation * contrast) * DB[k, :] ===
  if (verbose >= 2) Rcout << "  Computing sum of outer products..." << std::endl;
  
  // Pre-compute H.rotation * contrast (L x 1 vector)
  VectorXd H_contrast = H_rotation * contrast;
  
  MatrixXd diffExp = MatrixXd::Zero(J, N);
  
  // Add effect for k=1 (already fetched)
  VectorXd effect_1 = Rk_first * H_contrast;
  diffExp += effect_1 * DB.row(0);
  
  // Loop for k = 2 to K
  for (int k = 1; k < K; ++k) {
    if (verbose >= 2 && K <= 20) {
      Rcout << "    Factor " << (k + 1) << "/" << K << std::endl;
    }
    
    Eigen::MatrixXd Rk = as<Eigen::MatrixXd>(Rk_list[k]); // J x L
    
    // Compute effect for this factor: Rk * H_contrast (J x 1)
    VectorXd effect_k = Rk * H_contrast;
    
    // Add outer product: effect_k (Jx1) with DB.row(k) (1xN)
    diffExp += effect_k * DB.row(k);
  }
  
  if (verbose >= 1) {
    double mean_val = diffExp.mean();
    double std_val = std::sqrt((diffExp.array() - mean_val).square().mean());
    Rcout << "[OK] diffExp computed: mean = " << mean_val 
          << ", std = " << std_val << std::endl;
  }
  
  return diffExp;
}