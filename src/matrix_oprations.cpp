#include <RcppEigen.h>
#include <Eigen/Sparse>
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMat;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXd Mat;


//' Compute Cell-Specific Scaling Factors (Dense)
//'
//' Calculates column-wise means of a dense expression matrix for hyperparameter initialization.
//' This function computes s_0 = (J_vec^T * Y) / J, where each element represents the mean
//' expression across all genes for each cell. Used internally during GEDI model setup.
//'
//' @param J_vec Vector of ones with length equal to number of genes (J)
//' @param Y Dense log-transformed expression matrix (genes x cells)
//' @param J Number of genes (rows in Y)
//'
//' @return Vector of length N (number of cells) containing cell-specific scaling factors
//'
//' @details
//' This is the dense matrix version used when Y matrix is provided directly.
//' For count matrices (M), use the sparse version \code{compute_s_0()} instead.
//' The function uses BLAS-optimized dense matrix-vector multiplication for efficiency.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd compute_s_0_dense(const Eigen::Map<Eigen::VectorXd>& J_vec,
                                  const Eigen::Map<Eigen::MatrixXd>& Y,
                                  double J) {
  // Use BLAS-optimized dense matrix-vector multiplication
  // J_vec^T * Y is equivalent to Y^T * J_vec
  Vec s_0 = Y.transpose() * J_vec;

  // Add constant and normalize by J
  s_0 = (s_0.array()) / J;

  return s_0;
}

//' Compute Library-Size Normalized Expression (Sparse)
//'
//' Normalizes sparse expression matrix Y by dividing each element by the outer product
//' of J_vec and s_0. This is equivalent to Y / (J_vec * s_0^T).
//'
//' @param Y Sparse expression matrix (genes x cells)
//' @param J_vec Vector of ones with length J (number of genes)
//' @param s_0 Vector of cell-specific scaling factors (length N)
//'
//' @return Sparse matrix Yp with normalized expression values
//'
//' @details
//' Efficiently iterates through non-zero elements only, preserving sparsity.
//' Each Y(i,j) is divided by (J_vec(i) * s_0(j)).
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::SparseMatrix<double> compute_Yp(const Eigen::MappedSparseMatrix<double>& Y,
                                       const Eigen::Map<Eigen::VectorXd>& J_vec,
                                       const Eigen::Map<Eigen::VectorXd>& s_0) {
  // Create result sparse matrix with same structure as Y
  SpMat Yp = Y;

  // Iterate through non-zero elements of Y
  for (int k = 0; k < Yp.outerSize(); ++k) {
    for (SpMat::InnerIterator it(Yp, k); it; ++it) {
      double denominator = J_vec(it.row()) * s_0(it.col());
      it.valueRef() = it.value() / denominator;
    }
  }

  Yp.makeCompressed();
  return Yp;
}

//' Compute Gene-Specific Offsets (Dense)
//'
//' Calculates row-wise means of normalized expression matrix for hyperparameter initialization.
//' Computes o_0 = (Yp * N_vec) / N, where each element represents the mean expression
//' of a gene across all cells after library size normalization.
//'
//' @param Yp Normalized dense expression matrix (genes x cells)
//' @param N_vec Vector of ones with length N (number of cells)
//' @param N Number of cells (columns in Yp)
//'
//' @return Vector of length J (number of genes) containing gene-specific offsets
//'
//' @details
//' Used in conjunction with \code{compute_s_0_dense()} and \code{compute_Yp_dense()}
//' to initialize hyperparameters o_0 and s_0 for the GEDI model.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd compute_o_0_dense(const Eigen::Map<Eigen::MatrixXd>& Yp,
                                  const Eigen::Map<Eigen::VectorXd>& N_vec,
                                  double N) {
  // Compute o_0 = Yp * N_vec / N
  Vec o_0 = Yp * N_vec / N;

  return o_0;
}


//' Compute Library-Size Normalized Expression (Dense)
//'
//' Computes normalized expression Yp by subtracting the outer product of J_vec and s_0 from Y.
//' This is equivalent to Yp = Y - J_vec * s_0^T.
//'
//' @param Y Dense expression matrix (genes x cells)
//' @param J_vec Vector of ones with length J (number of genes)
//' @param s_0 Vector of cell-specific scaling factors (length N)
//'
//' @return Dense matrix Yp with normalized expression values
//'
//' @details
//' Alternative formulation of library size normalization for dense matrices.
//' Uses \code{noalias()} for efficient in-place subtraction without temporary objects.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd compute_Yp_dense(const Eigen::Map<Eigen::MatrixXd>& Y,
                                 const Eigen::Map<Eigen::VectorXd>& J_vec,
                                 const Eigen::Map<Eigen::VectorXd>& s_0) {
  Mat Yp = Y;
  // Using noalias() for efficient in-place subtraction
  Yp.noalias() -= J_vec * s_0.transpose();
  return Yp;
}

//' Vector Outer Product
//'
//' Computes the outer product of two vectors: a * b^T, returning a dense matrix.
//' This is used internally for efficient matrix operations in GEDI model fitting.
//'
//' @param a First vector (length m)
//' @param b Second vector (length n)
//'
//' @return Dense matrix of dimensions (m x n) containing the outer product
//'
//' @details
//' Uses BLAS-optimized operations for efficient computation.
//' Equivalent to R's \code{a \%*\% t(b)}.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd VecVecProduct(const Eigen::Map<Eigen::VectorXd>& a,
                                   const Eigen::Map<Eigen::VectorXd>& b) {
  // Direct outer product using BLAS-optimized operations
  return a * b.transpose();
}

//' Dense Matrix-Vector Product
//'
//' Computes the product of a dense matrix A and vector b: A * b.
//' Returns the result as a column vector.
//'
//' @param A Dense matrix (m x n)
//' @param b Vector of length n
//'
//' @return Vector of length m containing the product A * b
//'
//' @details
//' Uses BLAS level-2 optimized matrix-vector multiplication.
//' Equivalent to R's \code{A \%*\% b}.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd MatVecProduct(const Eigen::Map<Eigen::MatrixXd>& A,
                                   const Eigen::Map<Eigen::VectorXd>& b) {
  // Direct matrix-vector multiplication using BLAS
  return A * b;
}

//' Sparse Matrix-Vector Product
//'
//' Computes the product of a sparse matrix A and vector b: A * b.
//' Efficiently handles sparse matrices by only operating on non-zero elements.
//'
//' @param A Sparse matrix (m x n)
//' @param b Vector of length n
//'
//' @return Vector of length m containing the product A * b
//'
//' @details
//' Optimized for sparse matrices, avoiding unnecessary multiplications by zero.
//' Useful for large-scale single-cell data where most entries are zero.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd eigenSparseMatVecProduct(const Eigen::MappedSparseMatrix<double>& A,
                                         const Eigen::Map<Eigen::VectorXd>& b) {
  return A * b;
}


//' Compute Library-Size Normalized Count Matrix (Sparse)
//'
//' Normalizes sparse count matrix M by dividing each element by the outer product
//' of J_vec and s_0. This is equivalent to Mp = M / (J_vec * s_0^T).
//'
//' @param M Sparse raw count matrix (genes x cells)
//' @param J_vec Vector of ones with length J (number of genes)
//' @param s_0 Vector of cell-specific scaling factors (length N)
//'
//' @return Sparse matrix Mp with normalized counts
//'
//' @details
//' Efficiently iterates through non-zero elements only, preserving sparsity structure.
//' Each M(i,j) is divided by (J_vec(i) * s_0(j)).
//' Used for initializing hyperparameters from raw count data.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::SparseMatrix<double> compute_Mp(const Eigen::MappedSparseMatrix<double> M,
                                       const Eigen::Map<Eigen::VectorXd> J_vec,
                                       const Eigen::Map<Eigen::VectorXd> s_0) {
  // Create result sparse matrix with same structure as M
  SpMat Mp = M;

  // Compute the outer product J_vec * s_0^T efficiently
  // We'll divide each element M(i,j) by (J_vec(i) * s_0(j))

  // Iterate through non-zero elements of M
  for (int k = 0; k < Mp.outerSize(); ++k) {
    for (SpMat::InnerIterator it(Mp, k); it; ++it) {
      // it.row() = i, it.col() = j, it.value() = M(i,j)
      double denominator = J_vec(it.row()) * s_0(it.col());
      it.valueRef() = it.value() / denominator;
    }
  }

  Mp.makeCompressed();
  return Mp;
}


//' Compute Cell-Specific Scaling Factors from Counts (Sparse)
//'
//' Calculates column-wise sums of sparse count matrix M for hyperparameter initialization.
//' Computes s_0 = (J_vec^T * M + 0.01) / J, where each element represents the library
//' size (total count) for each cell with a small constant added for numerical stability.
//'
//' @param J_vec Vector of ones with length equal to number of genes (J)
//' @param M Sparse raw count matrix (genes x cells)
//' @param J Number of genes (rows in M)
//'
//' @return Vector of length N (number of cells) containing cell-specific scaling factors
//'
//' @details
//' This is the sparse matrix version used when raw count matrix M is provided.
//' For pre-processed Y matrices, use \code{compute_s_0_dense()} instead.
//' The small constant (0.01) prevents division by zero for cells with very low counts.
//' Uses BLAS-optimized sparse matrix-vector multiplication.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd compute_s_0(const Eigen::Map<Eigen::VectorXd> J_vec,
                            const Eigen::MappedSparseMatrix<double> M,
                            double J) {
  // Use BLAS-optimized sparse-dense multiplication
  // J_vec^T * M is equivalent to M.transpose() * J_vec
  Vec s_0 = M.transpose() * J_vec;

  // Add constant and normalize by J
  s_0 = (s_0.array() + 0.01) / J;

  return s_0;
}


//' Compute Gene-Specific Offsets from Normalized Counts (Sparse)
//'
//' Calculates row-wise sums of library-size normalized sparse matrix Mp for
//' hyperparameter initialization. Computes o_0 = (Mp * N_vec + 1e-5) / N, where
//' each element represents the mean normalized count of a gene across all cells.
//'
//' @param Mp Sparse library-size normalized count matrix (genes x cells)
//' @param N_vec Vector of ones with length N (number of cells)
//' @param N Number of cells (columns in Mp)
//'
//' @return Vector of length J (number of genes) containing gene-specific offsets
//'
//' @details
//' Used in conjunction with \code{compute_s_0()} and \code{compute_Mp()} to initialize
//' hyperparameters o_0 and s_0 from raw count data. The small constant (1e-5) provides
//' numerical stability for genes with very low expression.
//' Uses BLAS-optimized sparse matrix-vector multiplication.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd compute_o_0(const Eigen::MappedSparseMatrix<double> Mp,
                            const Eigen::Map<Eigen::VectorXd> N_vec,
                            double N) {
  // Use BLAS-optimized sparse matrix-vector multiplication
  Vec o_0 = Mp * N_vec;

  // Add constant and normalize by N
  o_0 = (o_0.array() + 1e-5) / N;

  return o_0;
}
