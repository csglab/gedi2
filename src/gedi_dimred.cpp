// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

/**
 * Compute Factorized SVD (Internal C++)
 *
 * Performs factorized SVD preserving GEDI structure: SVD(Z) x SVD(middle) x SVD(DB).
 * This function computes DB internally from the model parameters (D, Bi_list).
 * This is used for the standard, cached SVD/PCA embeddings.
 *
 * @param Z Shared metagene matrix (J x K)
 * @param D Scaling vector (length K)
 * @param Bi_list List of sample-specific cell projection matrices (K x Ni each)
 * @param verbose Integer verbosity level (0 = silent, 1 = progress, 2 = detailed)
 *
 * @return List with three components:
 * \itemize{
 * \item d: Singular values (length K vector)
 * \item u: Left singular vectors (J x K matrix, genes x factors)
 * \item v: Right singular vectors (N x K matrix, cells x factors)
 * }
 *
 * @keywords internal
 * @noRd
 */
// [[Rcpp::export]]
Rcpp::List compute_svd_factorized_cpp(
    const Eigen::Map<Eigen::MatrixXd>& Z,
    const Eigen::Map<Eigen::VectorXd>& D,
    const Rcpp::List& Bi_list,
    int verbose = 0
) {
  
  // ========================================================================
  // Dimension validation and setup
  // ========================================================================

  // int J = Z.rows();  // Unused local variable
  int K = Z.cols();
  int numSamples = Bi_list.size();
  
  if (D.size() != K) {
    stop("Dimension mismatch: D must have length K");
  }
  
  if (numSamples == 0) {
    stop("Bi_list cannot be empty");
  }
  
  // Count total cells
  int N = 0;
  std::vector<int> Ni(numSamples);
  
  for (int i = 0; i < numSamples; ++i) {
    Eigen::MatrixXd Bi = as<Eigen::MatrixXd>(Bi_list[i]);
    
    if (Bi.rows() != K) {
      stop("Dimension mismatch: Bi[%d] must have K rows", i + 1);
    }
    
    Ni[i] = Bi.cols();
    N += Ni[i];
  }
  
  // ========================================================================
  // Step 1: SVD of Z
  // ========================================================================
  
  JacobiSVD<MatrixXd> svd_Z(Z, ComputeThinU | ComputeThinV);
  MatrixXd U_Z = svd_Z.matrixU();        // J x K
  VectorXd S_Z = svd_Z.singularValues(); // K
  MatrixXd V_Z = svd_Z.matrixV();        // K x K
  
  // ========================================================================
  // Step 2: Concatenate B and compute DB
  // ========================================================================
  
  MatrixXd B(K, N);
  int col_offset = 0;
  
  for (int i = 0; i < numSamples; ++i) {
    Eigen::MatrixXd Bi = as<Eigen::MatrixXd>(Bi_list[i]);
    int Ni_current = Bi.cols();
    B.block(0, col_offset, K, Ni_current) = Bi;
    col_offset += Ni_current;
  }
  
  // Apply diagonal scaling: DB = diag(D) * B
  MatrixXd DB = D.asDiagonal() * B;  // K x N
  
  // ========================================================================
  // Step 3: SVD of DB
  // ========================================================================
  
  JacobiSVD<MatrixXd> svd_DB(DB, ComputeThinU | ComputeThinV);
  MatrixXd U_DB = svd_DB.matrixU();        // K x K
  VectorXd S_DB = svd_DB.singularValues(); // K
  MatrixXd V_DB = svd_DB.matrixV();        // N x K
  
  // ========================================================================
  // Step 4: Form and decompose middle matrix
  // ========================================================================
  
  // Middle = diag(S_Z) * V_Z^T * U_DB * diag(S_DB)
  MatrixXd middle = S_Z.asDiagonal() * V_Z.transpose() * U_DB * S_DB.asDiagonal();
  
  // ========================================================================
  // Step 5: SVD of middle matrix
  // ========================================================================

  JacobiSVD<MatrixXd> svd_middle(middle, ComputeThinU | ComputeThinV);
  MatrixXd U_middle = svd_middle.matrixU();        // K x K
  VectorXd S_middle = svd_middle.singularValues(); // K
  MatrixXd V_middle = svd_middle.matrixV();        // K x K

  // ========================================================================
  // Reconstruct final SVD and return
  // ========================================================================

  MatrixXd u = U_Z * U_middle;      // J x K
  MatrixXd v = V_DB * V_middle;     // N x K
  VectorXd d = S_middle;            // K
  
  return List::create(
    Named("d") = d,
    Named("u") = u,
    Named("v") = v
  );
}


/**
 * Run Factorized SVD on a Custom Projection (Internal C++)
 *
 * Performs the GEDI 3-stage factorized SVD: SVD(Z) x SVD(middle) x SVD(projDB).
 * This general-purpose version takes a *pre-computed* projDB matrix as input.
 * This is used by all dynamics functions (vector field, gradient) to analyze
 * hypothetical or transition states.
 *
 * @param Z Shared metagene matrix (J x K)
 * @param projDB Custom projected cell matrix (K x N_total), where N_total
 * could be N, 2N, 3N, etc. depending on the dynamic analysis.
 * @param verbose Integer verbosity level (0 = silent, 1 = progress, 2 = detailed)
 *
 * @return List with three components:
 * \itemize{
 * \item d: Singular values (length K vector)
 * \item u: Left singular vectors (J x K matrix, genes x factors)
 * \item v: Right singular vectors (N_total x K matrix, custom_cells x factors)
 * }
 *
 * @keywords internal
 * @noRd
 */
// [[Rcpp::export]]
Rcpp::List run_factorized_svd_cpp(
    const Eigen::Map<Eigen::MatrixXd>& Z,
    const Eigen::Map<Eigen::MatrixXd>& projDB,
    int verbose = 0
) {
  
  // ========================================================================
  // Dimension validation
  // ========================================================================

  // int J = Z.rows();  // Unused local variable
  int K = Z.cols();
  int K_proj = projDB.rows();
  // int N_total = projDB.cols();  // Unused local variable

  if (K != K_proj) {
    stop("Dimension mismatch: Z columns (%d) must equal projDB rows (%d)", K, K_proj);
  }

  // ========================================================================
  // Step 1: SVD of Z
  // ========================================================================
  
  JacobiSVD<MatrixXd> svd_Z(Z, ComputeThinU | ComputeThinV);
  MatrixXd U_Z = svd_Z.matrixU();        // J x K
  VectorXd S_Z = svd_Z.singularValues(); // K
  MatrixXd V_Z = svd_Z.matrixV();        // K x K
  
  // ========================================================================
  // Step 2: SVD of projDB (Input matrix)
  // ========================================================================
  
  JacobiSVD<MatrixXd> svd_projDB(projDB, ComputeThinU | ComputeThinV);
  MatrixXd U_DB = svd_projDB.matrixU();        // K x K
  VectorXd S_DB = svd_projDB.singularValues(); // K
  MatrixXd V_DB = svd_projDB.matrixV();        // N_total x K
  
  // ========================================================================
  // Step 3: Form and decompose middle matrix
  // ========================================================================

  // Middle = diag(S_Z) * V_Z^T * U_DB * diag(S_DB)
  MatrixXd middle = S_Z.asDiagonal() * V_Z.transpose() * U_DB * S_DB.asDiagonal();

  // ========================================================================
  // Step 4: SVD of middle matrix
  // ========================================================================

  JacobiSVD<MatrixXd> svd_middle(middle, ComputeThinU | ComputeThinV);
  MatrixXd U_middle = svd_middle.matrixU();        // K x K
  VectorXd S_middle = svd_middle.singularValues(); // K
  MatrixXd V_middle = svd_middle.matrixV();        // K x K

  // ========================================================================
  // Reconstruct final SVD and return
  // ========================================================================

  MatrixXd u = U_Z * U_middle;      // J x K
  MatrixXd v = V_DB * V_middle;     // N_total x K
  VectorXd d = S_middle;            // K
  
  return List::create(
    Named("d") = d,
    Named("u") = u,
    Named("v") = v
  );
}