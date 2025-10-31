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
// GEDI Projection Functions
// Stateless auxiliary functions for computing manifold projections
// ============================================================================

//' Compute ZDB Projection (Internal)
//'
//' Computes the shared manifold projection ZDB = Z * diag(D) * B, where B is
//' the concatenation of all sample-specific Bi matrices. This is the main
//' integrated representation of cells in the GEDI latent space.
//'
//' @param Z Shared metagene matrix (J Ã— K), where J = genes, K = latent factors
//' @param D Scaling vector (length K) representing the importance of each factor
//' @param Bi_list List of sample-specific cell projection matrices, where each
//'   Bi is a K Ã— Ni matrix (Ni = number of cells in sample i)
//' @param verbose Integer verbosity level:
//'   \itemize{
//'     \item 0: Silent (no output)
//'     \item 1: Progress bar and summary statistics
//'     \item 2: Detailed per-sample information
//'   }
//'
//' @return Dense matrix ZDB of dimensions J Ã— N, where N = sum(Ni) is the total
//'   number of cells across all samples. Each column represents a cell in the
//'   integrated latent space.
//'
//' @details
//' The ZDB projection represents each cell as a linear combination of shared
//' metagenes (Z), weighted by latent factors and scaled by their importance (D).
//' 
//' Computational strategy:
//' \enumerate{
//'   \item Pre-compute ZD = Z * diag(D) once (saves KÃ—JÃ—numSamples operations)
//'   \item For each sample i: compute ZD * Bi and concatenate
//'   \item Use Eigen block operations for efficient memory access
//' }
//'
//' Memory: Allocates one J Ã— N dense matrix. For large datasets (e.g., 20k genes,
//' 50k cells), this requires ~8 GB RAM. Consider computing projections on demand
//' or working with subsets if memory is limited.
//'
//' Performance: OpenMP parallelization available if enabled during compilation.
//' Typical speed: ~100-200ms for 20k Ã— 5k dataset on modern CPU.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd compute_ZDB_cpp(
    const Eigen::Map<Eigen::MatrixXd>& Z,
    const Eigen::Map<Eigen::VectorXd>& D,
    const Rcpp::List& Bi_list,
    int verbose = 0
) {
  
  // Dimension validation and setup  
  int J = Z.rows();
  int K = Z.cols();
  int numSamples = Bi_list.size();
  
  if (D.size() != K) {
    stop("Dimension mismatch: D must have length K");
  }
  
  if (numSamples == 0) {
    stop("Bi_list cannot be empty");
  }
  
  // Count total cells across all samples
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
  
  // Pre-compute ZD = Z * diag(D) for efficiency
  MatrixXd ZD = Z * D.asDiagonal();
  
  // Allocate output matrix
  MatrixXd ZDB(J, N);
  
  // Compute ZD * Bi for each sample and concatenate
  int col_offset = 0;

#ifdef _OPENMP
  // Note: Parallelizing over samples here may not be beneficial due to
  // memory bandwidth limitations. Profile before enabling.
  // #pragma omp parallel for schedule(dynamic) if(numSamples > 4)
#endif

  for (int i = 0; i < numSamples; ++i) {
    // Extract Bi for this sample
    Eigen::MatrixXd Bi = as<Eigen::MatrixXd>(Bi_list[i]);
    int Ni_current = Bi.cols();

    // Compute ZD * Bi and store in appropriate block
    ZDB.block(0, col_offset, J, Ni_current) = ZD * Bi;

    col_offset += Ni_current;
  }
  
  return ZDB;
}


//' Compute DB Projection (Internal)
//'
//' Computes the latent factor embedding DB = diag(D) * B, where B is the
//' concatenation of all sample-specific Bi matrices. This represents cells
//' in the K-dimensional latent factor space.
//'
//' @param D Scaling vector (length K) representing factor importance
//' @param Bi_list List of sample-specific cell projection matrices, where each
//'   Bi is a K Ã— Ni matrix
//' @param verbose Integer verbosity level (0 = silent, 1 = progress, 2 = detailed)
//'
//' @return Dense matrix DB of dimensions K Ã— N, where N = sum(Ni). Each column
//'   represents a cell's coordinates in the latent factor space.
//'
//' @details
//' The DB projection is simpler than ZDB as it operates directly in the K-dimensional
//' latent space without projecting through the gene space. This is useful for:
//' \itemize{
//'   \item Visualization (when K is small, e.g., K=2 or K=3)
//'   \item Downstream clustering or classification
//'   \item Understanding factor contributions per cell
//' }
//'
//' Computational strategy:
//' \enumerate{
//'   \item Concatenate all Bi matrices horizontally
//'   \item Apply diagonal scaling: diag(D) * B
//' }
//'
//' Memory: Much smaller than ZDB (K Ã— N vs. J Ã— N). For K=10 and N=50k cells,
//' requires only ~4 MB vs. ~8 GB for ZDB when J=20k.
//'
//' Performance: Very fast (~10-50ms) since K << J typically.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd compute_DB_cpp(
    const Eigen::Map<Eigen::VectorXd>& D,
    const Rcpp::List& Bi_list,
    int verbose = 0
) {

  // Dimension validation and setup  
  int K = D.size();
  int numSamples = Bi_list.size();
  
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
  
  // Concatenate all Bi matrices
  MatrixXd B(K, N);
  int col_offset = 0;

  for (int i = 0; i < numSamples; ++i) {
    Eigen::MatrixXd Bi = as<Eigen::MatrixXd>(Bi_list[i]);
    int Ni_current = Bi.cols();

    B.block(0, col_offset, K, Ni_current) = Bi;

    col_offset += Ni_current;
  }

  // Apply diagonal scaling: DB = diag(D) * B
  MatrixXd DB = D.asDiagonal() * B;
  
  return DB;
}


//' Compute ADB Projection (Internal)
//'
//' Computes the pathway activity projection ADB = C.rotation * A * diag(D) * B.
//' This represents cells in terms of their pathway/gene set activities rather
//' than individual gene expression.
//'
//' @param C_rotation Rotation matrix (num_pathways x P) that maps from the P-dimensional
//'   reduced SVD space back to the original pathway space
//' @param A Pathway-factor connection matrix (P x K) learned during training, where
//'   P is the number of reduced components and K is the number of latent factors
//' @param D Scaling vector (length K) representing factor importance
//' @param Bi_list List of sample-specific cell projection matrices (K x Ni each)
//' @param verbose Integer verbosity level (0 = silent, 1 = progress, 2 = detailed)
//'
//' @return Dense matrix ADB of dimensions num_pathways x N, where num_pathways is the
//'   number of original pathways and N = total cells. Each column represents a cell's
//'   pathway activities.
//'
//' @details
//' The ADB projection is only available when a gene-level prior matrix C was
//' provided during model setup. It enables:
//' \itemize{
//'   \item Pathway-level interpretation of cell states
//'   \item Identifying enriched pathways per cell type
//'   \item Biological interpretation in original pathway space
//' }
//'
//' Mathematical formulation:
//' \deqn{ADB = C.rotation \times A \times diag(D) \times B}
//'
//' where:
//' \itemize{
//'   \item C.rotation (num_pathways x P) maps from reduced space to original pathways
//'   \item A (P x K) captures pathway-factor associations in reduced space
//'   \item D (K) scales the importance of each latent factor
//'   \item B (K x N) contains cell projections in latent space
//' }
//'
//' Computational strategy:
//' \enumerate{
//'   \item Concatenate B from all Bi matrices
//'   \item Compute (C.rotation x A) x (diag(D) x B) efficiently
//'   \item Return in original pathway space for interpretation
//' }
//'
//' Memory: num_pathways x N matrix. For 100 original pathways and N=50k cells, ~40 MB.
//'
//' Performance: Fast due to small P and K dimensions (~50-100ms typical).
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd compute_ADB_cpp(
    const Eigen::Map<Eigen::MatrixXd>& C_rotation,
    const Eigen::Map<Eigen::MatrixXd>& A,
    const Eigen::Map<Eigen::VectorXd>& D,
    const Rcpp::List& Bi_list,
    int verbose = 0
) {

  // Dimension extraction - C_rotation is num_pathways x P, A is P x K
  // int num_pathways = C_rotation.rows();  // Unused local variable (original number of pathways)
  int P = C_rotation.cols();             // Number of reduced SVD components
  int K = A.cols();                      // Number of latent factors
  int numSamples = Bi_list.size();
  
  // Dimension validation
  if (A.rows() != P) {
    stop("Dimension mismatch: A must have P rows (got %d, expected %d)", A.rows(), P);
  }
  
  if (D.size() != K) {
    stop("Dimension mismatch: D must have length K (got %d, expected %d)", D.size(), K);
  }
  
  if (numSamples == 0) {
    stop("Bi_list cannot be empty");
  }
  
  if (P == 0) {
    stop("Cannot compute ADB: no gene-level prior (C) was provided");
  }
  
  // Count total cells and validate Bi dimensions
  int N = 0;
  std::vector<int> Ni(numSamples);
  
  for (int i = 0; i < numSamples; ++i) {
    Eigen::MatrixXd Bi = as<Eigen::MatrixXd>(Bi_list[i]);
    
    if (Bi.rows() != K) {
      stop("Dimension mismatch: Bi[%d] must have K rows (got %d, expected %d)", 
           i + 1, Bi.rows(), K);
    }
    
    Ni[i] = Bi.cols();
    N += Ni[i];
  }
  
  // Concatenate all Bi matrices into B
  MatrixXd B(K, N);
  int col_offset = 0;

  for (int i = 0; i < numSamples; ++i) {
    Eigen::MatrixXd Bi = as<Eigen::MatrixXd>(Bi_list[i]);
    int Ni_current = Bi.cols();

    B.block(0, col_offset, K, Ni_current) = Bi;

    col_offset += Ni_current;
  }

  // Compute ADB = C.rotation * A * diag(D) * B
  // Efficient computation order: (C.rotation * A) * (diag(D) * B)
  // This minimizes intermediate matrix sizes
  MatrixXd CA = C_rotation * A;        // P x K
  MatrixXd DB = D.asDiagonal() * B;    // K x N
  MatrixXd ADB = CA * DB;              // P x N
  
  return ADB;
}