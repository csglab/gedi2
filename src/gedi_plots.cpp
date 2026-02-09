// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <map>
#include <cmath>

using namespace Rcpp;
using namespace Eigen;

//' Compute Single Feature Projection (Internal)
//'
//' Projects a single feature (gene) through the GEDI model to get cell-level values.
//' Computes: (feature_weights * D) %*% B, avoiding full ZDB computation.
//'
//' @param feature_weights Vector of length K (factor loadings for this feature)
//' @param D Scaling vector of length K
//' @param Bi_list List of sample-specific cell projection matrices (K x Ni each)
//' @param verbose Integer verbosity level
//'
//' @return Vector of length N (total cells) with projected values
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd compute_feature_projection(
    const Eigen::Map<Eigen::VectorXd>& feature_weights,
    const Eigen::Map<Eigen::VectorXd>& D,
    const Rcpp::List& Bi_list,
    int verbose = 0
) {
  
  int K = feature_weights.size();
  int numSamples = Bi_list.size();
  
  if (D.size() != K) {
    stop("Dimension mismatch: D must have length K");
  }
  
  if (verbose >= 1) {
    Rcout << "Computing feature projection..." << std::endl;
  }
  
  // Compute weighted_D = feature_weights * D (element-wise)
  VectorXd weighted_D = feature_weights.array() * D.array();
  
  // Count total cells
  int N = 0;
  for (int i = 0; i < numSamples; ++i) {
    MatrixXd Bi = as<MatrixXd>(Bi_list[i]);
    N += Bi.cols();
  }
  
  // Allocate result
  VectorXd result(N);
  int offset = 0;
  
  // For each sample: compute weighted_D^T * Bi
  for (int i = 0; i < numSamples; ++i) {
    MatrixXd Bi = as<MatrixXd>(Bi_list[i]);
    int Ni = Bi.cols();
    
    // result[offset:(offset+Ni)] = weighted_D^T * Bi
    result.segment(offset, Ni) = weighted_D.transpose() * Bi;
    offset += Ni;
  }
  
  if (verbose >= 1) {
    Rcout << "Feature projection computed: " << N << " cells" << std::endl;
  }
  
  return result;
}


//' Compute Multi-Feature Projection (Internal)
//'
//' Projects multiple features through the GEDI model simultaneously.
//' Computes: (feature_weights * D) %*% B for F features.
//'
//' @param feature_weights Matrix K x F (factor loadings for F features)
//' @param D Scaling vector of length K
//' @param Bi_list List of sample-specific cell projection matrices (K x Ni each)
//' @param verbose Integer verbosity level
//'
//' @return Matrix N x F with projected values for each feature
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd compute_multi_feature_projection(
    const Eigen::Map<Eigen::MatrixXd>& feature_weights,
    const Eigen::Map<Eigen::VectorXd>& D,
    const Rcpp::List& Bi_list,
    int verbose = 0
) {
  
  int K = feature_weights.rows();
  int F = feature_weights.cols();
  int numSamples = Bi_list.size();
  
  if (D.size() != K) {
    stop("Dimension mismatch: D must have length K");
  }
  
  if (verbose >= 1) {
    Rcout << "Computing multi-feature projection: " << F << " features..." << std::endl;
  }
  
  // Count total cells
  int N = 0;
  for (int i = 0; i < numSamples; ++i) {
    MatrixXd Bi = as<MatrixXd>(Bi_list[i]);
    N += Bi.cols();
  }
  
  // Allocate result: N x F
  MatrixXd result(N, F);
  int offset = 0;
  
  // For each sample
  for (int i = 0; i < numSamples; ++i) {
    MatrixXd Bi = as<MatrixXd>(Bi_list[i]);
    int Ni = Bi.cols();
    
    // For each feature: compute (feature_weights[:,f] * D)^T * Bi
    for (int f = 0; f < F; ++f) {
      VectorXd weighted_D = feature_weights.col(f).array() * D.array();
      result.block(offset, f, Ni, 1) = (weighted_D.transpose() * Bi).transpose();
    }
    
    offset += Ni;
  }
  
  if (verbose >= 1) {
    Rcout << "Multi-feature projection computed: " << N << " cells x " 
          << F << " features" << std::endl;
  }
  
  return result;
}


//' Aggregate Vector Field into Bins (Internal)
//'
//' Bins vector field data into a grid and computes mean vectors per bin.
//' Used for cleaner arrow plots without overplotting.
//'
//' @param Dim1 Start x-coordinates (length N)
//' @param Dim2 Start y-coordinates (length N)
//' @param To1 End x-coordinates (length N)
//' @param To2 End y-coordinates (length N)
//' @param color Color values (length N), can be numeric or will be converted
//' @param alpha Alpha values (length N)
//' @param n_bins Number of bins per dimension
//' @param min_per_bin Minimum observations required per bin
//'
//' @return Data frame with columns: Dim1, Dim2, To1, To2, deltaDim1, deltaDim2,
//'   Color, Alpha, n (observations per bin)
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
DataFrame aggregate_vectors(
    const Eigen::Map<Eigen::VectorXd>& Dim1,
    const Eigen::Map<Eigen::VectorXd>& Dim2,
    const Eigen::Map<Eigen::VectorXd>& To1,
    const Eigen::Map<Eigen::VectorXd>& To2,
    const Eigen::Map<Eigen::VectorXd>& color,
    const Eigen::Map<Eigen::VectorXd>& alpha,
    int n_bins,
    int min_per_bin
) {
  
  int N = Dim1.size();
  
  if (Dim2.size() != N || To1.size() != N || To2.size() != N ||
      color.size() != N || alpha.size() != N) {
    stop("All input vectors must have the same length");
  }
  
  // Compute bin boundaries
  double x_min = Dim1.minCoeff();
  double x_max = Dim1.maxCoeff();
  double y_min = Dim2.minCoeff();
  double y_max = Dim2.maxCoeff();
  
  double x_step = (x_max - x_min) / n_bins;
  double y_step = (y_max - y_min) / n_bins;
  
  // Map to store bin aggregates
  // Key: bin_x * n_bins + bin_y
  std::map<int, std::vector<int>> bin_indices;
  
  // Assign points to bins
  for (int i = 0; i < N; ++i) {
    int bin_x = std::min(static_cast<int>((Dim1(i) - x_min) / x_step), n_bins - 1);
    int bin_y = std::min(static_cast<int>((Dim2(i) - y_min) / y_step), n_bins - 1);
    int bin_key = bin_x * n_bins + bin_y;
    
    bin_indices[bin_key].push_back(i);
  }
  
  // Aggregate bins
  std::vector<double> agg_Dim1, agg_Dim2, agg_To1, agg_To2;
  std::vector<double> agg_deltaDim1, agg_deltaDim2;
  std::vector<double> agg_Color, agg_Alpha;
  std::vector<int> agg_n;
  
  for (auto& kv : bin_indices) {
    std::vector<int>& indices = kv.second;
    int count = indices.size();
    
    if (count < min_per_bin) continue;
    
    // Compute means
    double sum_Dim1 = 0, sum_Dim2 = 0, sum_To1 = 0, sum_To2 = 0;
    double sum_Color = 0, sum_Alpha = 0;
    
    for (int idx : indices) {
      sum_Dim1 += Dim1(idx);
      sum_Dim2 += Dim2(idx);
      sum_To1 += To1(idx);
      sum_To2 += To2(idx);
      sum_Color += color(idx);
      sum_Alpha += alpha(idx);
    }
    
    double mean_Dim1 = sum_Dim1 / count;
    double mean_Dim2 = sum_Dim2 / count;
    double mean_To1 = sum_To1 / count;
    double mean_To2 = sum_To2 / count;
    
    agg_Dim1.push_back(mean_Dim1);
    agg_Dim2.push_back(mean_Dim2);
    agg_To1.push_back(mean_To1);
    agg_To2.push_back(mean_To2);
    agg_deltaDim1.push_back(mean_To1 - mean_Dim1);
    agg_deltaDim2.push_back(mean_To2 - mean_Dim2);
    agg_Color.push_back(sum_Color / count);
    agg_Alpha.push_back(sum_Alpha / count);
    agg_n.push_back(count);
  }
  
  return DataFrame::create(
    Named("Dim1") = agg_Dim1,
    Named("Dim2") = agg_Dim2,
    Named("To1") = agg_To1,
    Named("To2") = agg_To2,
    Named("deltaDim1") = agg_deltaDim1,
    Named("deltaDim2") = agg_deltaDim2,
    Named("Color") = agg_Color,
    Named("Alpha") = agg_Alpha,
    Named("n") = agg_n
  );
}