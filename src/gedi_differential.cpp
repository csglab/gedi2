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
 * @param contrast Vector of length num_covariates (= ncol(H_rotation))
 *   specifying the contrast in the user-facing covariate space
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

  // contrast lives in the original num_covariates space; H_rotation projects
  // num_covariates -> L. So contrast length must equal H_rotation.cols().
  if (contrast.size() != H_rotation.cols()) {
    stop("Dimension mismatch: contrast must have length ncol(H_rotation) "
         "(= number of original H covariates)");
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
 * @param contrast Vector of length num_covariates (= ncol(H_rotation))
 *   specifying the contrast in the user-facing covariate space
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

  // contrast lives in the original num_covariates space; H_rotation projects
  // num_covariates -> L. So contrast length must equal H_rotation.cols().
  if (contrast.size() != H_rotation.cols()) {
    stop("Dimension mismatch: contrast must have length ncol(H_rotation) "
         "(= number of original H covariates)");
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
 * @param contrast Vector of length num_covariates (= ncol(H_rotation))
 *   specifying the contrast in the user-facing covariate space
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

  // contrast lives in the original num_covariates space; H_rotation projects
  // num_covariates -> L. So contrast length must equal H_rotation.cols().
  if (contrast.size() != H_rotation.cols()) {
    stop("Dimension mismatch: contrast must have length ncol(H_rotation) "
         "(= number of original H covariates)");
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


/**
 * Compute Differential Pathway Activity (Internal C++)
 *
 * Computes the cell-specific differential pathway activity effect
 * (num_pathways x N). This is the differential analogue of compute_ADB_cpp,
 * exactly as getDiffExp_cpp is the differential analogue of compute_ZDB_cpp.
 *
 * Mathematical derivation:
 *   ADB = C.rotation * A * diag(D) * B, where solve_A() fits A linearly from Z:
 *     A = (Cmat^T Cmat + lambda I)^{-1} Cmat^T Z,  Cmat = inputC * C_rotation,
 *     lambda = 1 / S_A   (Cmat is the orthonormal SVD basis aux_C).
 *   Because solve_A is linear, the differential of ADB under the contrast is
 *     diffADB = ADB(Z + diffQ_Z) - ADB(Z)
 *             = C.rotation * solve_A(diffQ_Z) * diag(D) * B,
 *   i.e. apply the SAME shrinkage operator to the differential gene loadings
 *   diffQ_Z. Omitting the shrinkage would mis-scale diffADB relative to ADB by a
 *   factor of (1 + 1/lambda); we therefore replicate solve_A exactly.
 *
 * Computation order:
 *   1. H_contrast = H_rotation * contrast                         (L x 1)
 *   2. diffQ_Z    = [Rk * H_contrast for each k]                  (J x K)
 *   3. Cmat       = inputC * C_rotation                           (J x P, = aux_C)
 *   4. A_diff     = (Cmat^T Cmat + (1/S_A) I)^{-1} Cmat^T diffQ_Z (P x K)
 *   5. CA_diff    = C_rotation * A_diff                           (num_pathways x K)
 *   6. DB         = diag(D) * B                                   (K x N)
 *   7. diffADB    = CA_diff * DB                                  (num_pathways x N)
 *
 * The largest intermediates are J x K (diffQ_Z) and J x P (Cmat, bounded by the
 * size of the prior inputC). No J x N intermediate is ever materialised.
 *
 * @param Rk_list List of K matrices (each J x L), effect of sample variables
 *   on each latent factor k
 * @param H_rotation Rotation matrix (L x num_covariates)
 * @param contrast Vector of length num_covariates (= ncol(H_rotation))
 *   specifying the contrast in the user-facing covariate space
 * @param C_rotation Rotation matrix (num_pathways x P) that maps from the
 *   P-dimensional reduced SVD space to the original pathway space
 * @param inputC Gene-by-pathway input prior matrix (J x num_pathways)
 * @param D Scaling vector (length K) representing factor importance
 * @param Bi_list List of sample-specific cell projection matrices (K x Ni each)
 * @param S_A Pathway-loading shrinkage hyperparameter (= model$hyperparams$S_A);
 *   defines lambda = 1 / S_A in the solve_A operator
 * @param verbose Integer verbosity level (0 = silent, 1 = summary, 2 = detailed)
 *
 * @return Dense matrix diffADB of dimensions num_pathways x N representing the
 *   predicted differential pathway activity for each pathway in each cell.
 *
 * @keywords internal
 * @noRd
 */
// [[Rcpp::export]]
Eigen::MatrixXd getDiffADB_cpp(
    const Rcpp::List& Rk_list,
    const Eigen::Map<Eigen::MatrixXd>& H_rotation,
    const Eigen::Map<Eigen::VectorXd>& contrast,
    const Eigen::Map<Eigen::MatrixXd>& C_rotation,
    const Eigen::Map<Eigen::MatrixXd>& inputC,
    const Eigen::Map<Eigen::VectorXd>& D,
    const Rcpp::List& Bi_list,
    double S_A,
    int verbose = 0
) {

  // ---- Dimension validation ----
  int K = Rk_list.size();
  int L = H_rotation.rows();
  int num_pathways = C_rotation.rows();
  int P = C_rotation.cols();
  int numSamples = Bi_list.size();

  if (K == 0 || L == 0) {
    stop("Cannot compute diffADB: no sample-level prior (H) was provided");
  }

  if (P == 0) {
    stop("Cannot compute diffADB: no gene-level prior (C) was provided");
  }

  // contrast lives in the original num_covariates space; H_rotation projects
  // num_covariates -> L. So contrast length must equal H_rotation.cols().
  if (contrast.size() != H_rotation.cols()) {
    stop("Dimension mismatch: contrast must have length ncol(H_rotation) "
         "(= number of original H covariates)");
  }

  if (D.size() != K) {
    stop("Dimension mismatch: D must have length K");
  }

  if (S_A <= 0) {
    stop("Invalid hyperparameter: S_A must be positive");
  }

  // Validate first Rk and extract J
  Eigen::MatrixXd Rk_first = as<Eigen::MatrixXd>(Rk_list[0]);
  int J = Rk_first.rows();
  if (Rk_first.cols() != L) {
    stop("Dimension mismatch: Rk[1] must have L columns");
  }

  // inputC must be J x num_pathways
  if (inputC.rows() != J) {
    stop("Dimension mismatch: inputC must have J rows (got %d, expected %d)",
         inputC.rows(), J);
  }
  if (inputC.cols() != num_pathways) {
    stop("Dimension mismatch: inputC must have ncol(C_rotation) = %d columns "
         "(got %d)", num_pathways, inputC.cols());
  }

  if (numSamples == 0) {
    stop("Bi_list cannot be empty");
  }

  // Count total cells and validate Bi dimensions
  int N = 0;
  for (int i = 0; i < numSamples; ++i) {
    Eigen::MatrixXd Bi = as<Eigen::MatrixXd>(Bi_list[i]);
    if (Bi.rows() != K) {
      stop("Dimension mismatch: Bi[%d] must have K rows", i + 1);
    }
    N += Bi.cols();
  }

  if (verbose >= 1) {
    Rcout << "Computing diffADB (pathway-space): "
          << num_pathways << " pathways x " << N << " cells" << std::endl;
    if (verbose >= 2) {
      Rcout << "  Genes: " << J << ", Factors: " << K
            << ", Contrast dim: " << L
            << ", SVD components: " << P << std::endl;
    }
  }

  // ---- Step 1: H_contrast = H_rotation * contrast  (L x 1) ----
  VectorXd H_contrast = H_rotation * contrast;

  // ---- Step 2: diffQ_Z = [Rk * H_contrast for k=1..K]  (J x K) ----
  if (verbose >= 2) Rcout << "  Building diffQ_Z (J x K)..." << std::endl;

  MatrixXd diffQ_Z(J, K);
  diffQ_Z.col(0) = Rk_first * H_contrast;

  for (int k = 1; k < K; ++k) {
    Eigen::MatrixXd Rk = as<Eigen::MatrixXd>(Rk_list[k]); // J x L
    if (Rk.rows() != J || Rk.cols() != L) {
      stop("Dimension mismatch: all Rk matrices must have the same JxL dimensions");
    }
    diffQ_Z.col(k) = Rk * H_contrast;
  }

  // ---- Steps 3-5: differential pathway loadings via the EXACT solve_A operator ----
  // diffADB must be the differential of ADB, i.e. ADB(Z + diffQ_Z) - ADB(Z).
  // Since ADB = C.rotation * A * DB and solve_A() is linear in Z with
  //   A = (Cmat^T Cmat + lambda I)^{-1} Cmat^T Z,   Cmat = inputC * C_rotation,
  //   lambda = 1 / S_A   (see gedi_core.cpp::solve_A / workspace_CtC_inv),
  // the differential pathway loadings are A_diff = solve_A(diffQ_Z), applying
  // the SAME shrinkage so diffADB lands on the same scale as ADB. (Cmat is the
  // orthonormal SVD basis aux_C, so Cmat^T Cmat = I and the shrinkage reduces to
  // the scalar gamma = S_A/(S_A+1); we form CtC explicitly to mirror solve_A
  // exactly and stay robust to any deviation from orthonormality.)
  if (verbose >= 2) Rcout << "  Projecting through pathway space (solve_A operator)..." << std::endl;

  MatrixXd Cmat   = inputC * C_rotation;          // J x P  (= aux_C, orthonormal)
  MatrixXd CtC    = Cmat.transpose() * Cmat;      // P x P  (~ I_P)
  CtC.diagonal().array() += 1.0 / S_A;            // + lambda I,  lambda = 1 / S_A
  MatrixXd CtdQ   = Cmat.transpose() * diffQ_Z;   // P x K  (= Cmat^T * diffQ_Z)
  MatrixXd A_diff = CtC.ldlt().solve(CtdQ);       // P x K  (= workspace_CtC_inv * CtdQ)
  MatrixXd CA_diff = C_rotation * A_diff;         // num_pathways x K

  // ---- Steps 6-7: concatenate B, compute DB, then diffADB ----
  if (verbose >= 2) Rcout << "  Concatenating B and computing DB..." << std::endl;

  MatrixXd B(K, N);
  int col_offset = 0;
  for (int i = 0; i < numSamples; ++i) {
    Eigen::MatrixXd Bi = as<Eigen::MatrixXd>(Bi_list[i]);
    int Ni_current = Bi.cols();
    B.block(0, col_offset, K, Ni_current) = Bi;
    col_offset += Ni_current;
  }

  MatrixXd DB      = D.asDiagonal() * B; // K x N
  MatrixXd diffADB = CA_diff * DB;       // num_pathways x N

  if (verbose >= 1) {
    double mean_val = diffADB.mean();
    double std_val = std::sqrt((diffADB.array() - mean_val).square().mean());
    Rcout << "[OK] diffADB computed: mean = " << mean_val
          << ", std = " << std_val << std::endl;
  }

  return diffADB;
}
