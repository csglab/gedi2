// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <vector>
#include <string>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;
using namespace Eigen;

// ============================================================================
// Helper Functions for Yi Initialization from M
// ============================================================================

// Helper to compute log1p for sparse matrix -> dense result
MatrixXd sparse_log1p(const SparseMatrix<double>& M) {
  MatrixXd result = MatrixXd(M.rows(), M.cols());
  result.setZero();
  
  for (int k = 0; k < M.outerSize(); ++k) {
    for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
      result(it.row(), it.col()) = std::log1p(it.value());
    }
  }
  
  return result;
}

// Helper to compute log ratio for paired sparse matrices -> dense result
MatrixXd sparse_log_ratio(const SparseMatrix<double>& M1, const SparseMatrix<double>& M2) {
  MatrixXd result = MatrixXd(M1.rows(), M1.cols());
  
  for (int j = 0; j < M1.cols(); ++j) {
    for (int i = 0; i < M1.rows(); ++i) {
      double m1_val = M1.coeff(i, j);
      double m2_val = M2.coeff(i, j);
      result(i, j) = std::log((m1_val + 1.0) / (m2_val + 1.0));
    }
  }
  
  return result;
}

// ============================================================================
// Timer Class
// ============================================================================

class FunctionTimer {
public:
  FunctionTimer(const std::string& name, int verbose_level)
    : function_name(name), verbose(verbose_level), start_time(std::chrono::high_resolution_clock::now()) {}
  
  ~FunctionTimer() {
    if (verbose >= 4) {
      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
      Rcpp::Rcout << "      [TIMING] " << function_name << ": " << duration.count() << " us" << std::endl;
    }
  }
  
private:
  std::string function_name;
  int verbose;
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
};

// ============================================================================
// GEDI Class
// ============================================================================

class GEDI {
private:
  int J;
  int N;
  int K;
  int P;
  int L;
  int numSamples;
  int num_threads;
  
  std::vector<MatrixXd> params_Bi;
  std::vector<MatrixXd> params_Qi;
  std::vector<VectorXd> params_si;
  std::vector<VectorXd> params_oi;
  VectorXd params_o;
  MatrixXd params_Z;
  MatrixXd params_U;
  VectorXd params_S;
  VectorXd params_D;
  MatrixXd params_A;
  std::vector<MatrixXd> params_Rk;
  MatrixXd params_Ro;
  double params_sigma2;
  
  std::vector<MatrixXd> aux_ZDBi;
  std::vector<MatrixXd> aux_QiDBi;
  std::vector<MatrixXd> aux_Qi_hat;
  std::vector<VectorXd> aux_oi_hat;
  MatrixXd aux_C;
  MatrixXd aux_H;
  MatrixXd aux_diag_K;
  VectorXd aux_J_vec;
  std::vector<VectorXd> aux_Ni_vec;
  VectorXd aux_Ni;
  
  std::vector<MatrixXd> target_Yi;
  std::vector<SparseMatrix<double>> target_Mi;
  std::vector<SparseMatrix<double>> target_M1i;
  std::vector<SparseMatrix<double>> target_M2i;
  std::vector<MatrixXd> target_Xi;
  
  VectorXd hyperparams_S_Qi;
  VectorXd hyperparams_S_oi;
  double hyperparams_S_Z;
  double hyperparams_S_A;
  double hyperparams_S_R;
  double hyperparams_S_Qi_mean;
  double hyperparams_S_oi_mean;
  double hyperparams_S_si;
  double hyperparams_S_o;
  VectorXd hyperparams_o_0;
  std::vector<VectorXd> hyperparams_si_0;
  MatrixXd hyperparams_O;
  
  std::string obs_type;
  std::string mode;
  bool orthoZ;
  bool adjustD;
  bool is_si_fixed;
  int verbose;
  
  MatrixXd workspace_Y_res;
  MatrixXd workspace_B_concat;
  MatrixXd workspace_ZD;
  MatrixXd workspace_CtC_inv;
  MatrixXd workspace_HpHp_inv;
  MatrixXd workspace_Yp;

  // Pre-allocated tracking storage (prevents allocations in loop)
  MatrixXd tracking_prev_Z;
  MatrixXd tracking_prev_A;
  MatrixXd tracking_prev_Ro;
  VectorXd tracking_prev_o;
  std::vector<MatrixXd> tracking_prev_Bi;
  std::vector<MatrixXd> tracking_prev_Qi;
  std::vector<MatrixXd> tracking_prev_Rk;
  std::vector<VectorXd> tracking_prev_si;
  std::vector<VectorXd> tracking_prev_oi;

  std::vector<double> tracking_sigma2;
  std::vector<double> tracking_dZ;
  std::vector<double> tracking_dA;
  std::vector<double> tracking_do;
  std::vector<double> tracking_dRo;
  std::vector<std::vector<double>> tracking_dsi;
  std::vector<std::vector<double>> tracking_doi;
  std::vector<std::vector<double>> tracking_dBi;
  std::vector<std::vector<double>> tracking_dQi;
  std::vector<std::vector<double>> tracking_dRk;
  
  bool is_initialized;
  int total_iterations;

public:
  GEDI(const List& params,
       const List& aux,
       const List& target,
       const List& hyperparams,
       int verbose_ = 1,
       int num_threads_ = 0) :
    num_threads(num_threads_),
    verbose(verbose_),
    is_initialized(false),
    total_iterations(0) {
    
    J = as<int>(aux["J"]);
    N = as<int>(aux["N"]);
    K = as<int>(aux["K"]);
    P = as<int>(aux["P"]);
    L = as<int>(aux["L"]);
    numSamples = as<int>(aux["numSamples"]);
    
    obs_type = as<std::string>(aux["obs.type"]);
    mode = as<std::string>(aux["mode"]);
    orthoZ = as<bool>(aux["orthoZ"]);
    adjustD = as<bool>(aux["adjustD"]);
    is_si_fixed = as<bool>(aux["is_si_fixed"]);
    
    List Bi_list = params["Bi"];
    List Qi_list = params["Qi"];
    List si_list = params["si"];
    List oi_list = params["oi"];
    
    params_Bi.resize(numSamples);
    params_Qi.resize(numSamples);
    params_si.resize(numSamples);
    params_oi.resize(numSamples);
    
    for (int i = 0; i < numSamples; ++i) {
      params_Bi[i] = as<MatrixXd>(Bi_list[i]);
      params_Qi[i] = as<MatrixXd>(Qi_list[i]);
      params_si[i] = as<VectorXd>(si_list[i]);
      params_oi[i] = as<VectorXd>(oi_list[i]);
    }
    
    params_o = as<VectorXd>(params["o"]);
    params_Z = as<MatrixXd>(params["Z"]);
    params_U = as<MatrixXd>(params["U"]);
    params_S = as<VectorXd>(params["S"]);
    params_D = as<VectorXd>(params["D"]);
    params_sigma2 = as<double>(params["sigma2"]);
    
    if (P > 0) {
      params_A = as<MatrixXd>(params["A"]);
    } else {
      params_A = MatrixXd(0, 0);
    }
    
    if (L > 0) {
      List Rk_list = params["Rk"];
      params_Rk.resize(K);
      for (int k = 0; k < K; ++k) {
        params_Rk[k] = as<MatrixXd>(Rk_list[k]);
      }
      params_Ro = as<MatrixXd>(params["Ro"]);
    } else {
      params_Ro = MatrixXd(0, 0);
    }
    
    List ZDBi_list = aux["ZDBi"];
    List QiDBi_list = aux["QiDBi"];
    
    aux_ZDBi.resize(numSamples);
    aux_QiDBi.resize(numSamples);
    
    for (int i = 0; i < numSamples; ++i) {
      aux_ZDBi[i] = as<MatrixXd>(ZDBi_list[i]);
      aux_QiDBi[i] = as<MatrixXd>(QiDBi_list[i]);
    }
    
    if (L > 0) {
      List Qi_hat_list = aux["Qi_hat"];
      List oi_hat_list = aux["oi_hat"];
      aux_Qi_hat.resize(numSamples);
      aux_oi_hat.resize(numSamples);
      for (int i = 0; i < numSamples; ++i) {
        aux_Qi_hat[i] = as<MatrixXd>(Qi_hat_list[i]);
        aux_oi_hat[i] = as<VectorXd>(oi_hat_list[i]);
      }
    }
    
    if (P > 0) {
      aux_C = as<MatrixXd>(aux["C"]);
    } else {
      aux_C = MatrixXd(0, 0);
    }
    
    if (L > 0) {
      aux_H = as<MatrixXd>(aux["H"]);
    } else {
      aux_H = MatrixXd(0, 0);
    }
    
    aux_diag_K = as<MatrixXd>(aux["diag_K"]);
    aux_J_vec = as<VectorXd>(aux["J_vec"]);
    
    List Ni_vec_list = aux["Ni_vec"];
    aux_Ni_vec.resize(numSamples);
    for (int i = 0; i < numSamples; ++i) {
      aux_Ni_vec[i] = as<VectorXd>(Ni_vec_list[i]);
    }
    
    NumericVector Ni_r = aux["Ni"];
    aux_Ni = VectorXd(numSamples);
    for (int i = 0; i < numSamples; ++i) {
      aux_Ni(i) = Ni_r[i];
    }
    
    // ========================================================================
    // OPTION 2: Initialize Yi from M in C++ if Yi is empty
    // ========================================================================
    
    List Yi_list = target["Yi"];
    target_Yi.resize(numSamples);
    
    // Check if Yi is empty (R passed empty lists)
    bool Yi_is_empty = false;
    if (Yi_list.size() == 0) {
      Yi_is_empty = true;
    } else {
      SEXP first_Yi = Yi_list[0];
      if (Rf_isNull(first_Yi)) {
        Yi_is_empty = true;
      } else {
        // Check if it's an empty matrix
        NumericMatrix test_mat = as<NumericMatrix>(first_Yi);
        if (test_mat.nrow() == 0 || test_mat.ncol() == 0) {
          Yi_is_empty = true;
        }
      }
    }
    
    if (Yi_is_empty) {
      // Yi is empty - compute from M based on obs_type
      if (verbose >= 1) {
        Rcout << "Computing Yi from raw data in C++..." << std::endl;
      }
      
      if (obs_type == "M") {
        // Load Mi and compute Yi = log1p(Mi)
        List Mi_list = target["Mi"];
        target_Mi.resize(numSamples);
        for (int i = 0; i < numSamples; ++i) {
          target_Mi[i] = as<SparseMatrix<double>>(Mi_list[i]);
          target_Yi[i] = sparse_log1p(target_Mi[i]);
        }
        
      } else if (obs_type == "M_paired") {
        // Load M1i and M2i, compute Yi = log((M1+1)/(M2+1))
        List M1i_list = target["M1i"];
        List M2i_list = target["M2i"];
        target_M1i.resize(numSamples);
        target_M2i.resize(numSamples);
        for (int i = 0; i < numSamples; ++i) {
          target_M1i[i] = as<SparseMatrix<double>>(M1i_list[i]);
          target_M2i[i] = as<SparseMatrix<double>>(M2i_list[i]);
          target_Yi[i] = sparse_log_ratio(target_M1i[i], target_M2i[i]);
        }
        
      } else if (obs_type == "X") {
        // Load Xi and convert to Yi format
        List Xi_list = target["Xi"];
        target_Xi.resize(numSamples);
        for (int i = 0; i < numSamples; ++i) {
          target_Xi[i] = as<MatrixXd>(Xi_list[i]);
          target_Yi[i] = target_Xi[i]; // Copy
          
          // Convert: 1 -> 1, 0 -> -1, NA -> 0
          for (int j = 0; j < target_Yi[i].rows(); ++j) {
            for (int k = 0; k < target_Yi[i].cols(); ++k) {
              double val = target_Xi[i](j, k);
              if (std::isnan(val)) {
                target_Yi[i](j, k) = 0.0;
              } else if (val == 0.0) {
                target_Yi[i](j, k) = -1.0;
              } else {
                target_Yi[i](j, k) = 1.0;
              }
            }
          }
        }
        
      } else {
        // obs_type == "Y" - but Yi is empty, this is an error
        stop("Yi cannot be empty for observation type 'Y'");
      }
      
      if (verbose >= 1) {
        Rcout << "Yi computation complete." << std::endl;
      }
      
    } else {
      // Yi was provided - use it directly (backward compatibility)
      for (int i = 0; i < numSamples; ++i) {
        target_Yi[i] = as<MatrixXd>(Yi_list[i]);
      }
    }
    
    // Load Mi/M1i/M2i/Xi if not already loaded
    if (obs_type == "M" && target_Mi.size() == 0) {
      List Mi_list = target["Mi"];
      target_Mi.resize(numSamples);
      for (int i = 0; i < numSamples; ++i) {
        target_Mi[i] = as<SparseMatrix<double>>(Mi_list[i]);
      }
    } else if (obs_type == "M_paired" && target_M1i.size() == 0) {
      List M1i_list = target["M1i"];
      List M2i_list = target["M2i"];
      target_M1i.resize(numSamples);
      target_M2i.resize(numSamples);
      for (int i = 0; i < numSamples; ++i) {
        target_M1i[i] = as<SparseMatrix<double>>(M1i_list[i]);
        target_M2i[i] = as<SparseMatrix<double>>(M2i_list[i]);
      }
    } else if (obs_type == "X" && target_Xi.size() == 0) {
      List Xi_list = target["Xi"];
      target_Xi.resize(numSamples);
      for (int i = 0; i < numSamples; ++i) {
        target_Xi[i] = as<MatrixXd>(Xi_list[i]);
      }
    }
    
    hyperparams_S_Qi = as<VectorXd>(hyperparams["S_Qi"]);
    hyperparams_S_oi = as<VectorXd>(hyperparams["S_oi"]);
    hyperparams_S_Z = as<double>(hyperparams["S_Z"]);
    hyperparams_S_o = as<double>(hyperparams["S_o"]);
    hyperparams_S_si = as<double>(hyperparams["S_si"]);
    hyperparams_o_0 = as<VectorXd>(hyperparams["o_0"]);
    hyperparams_O = as<MatrixXd>(hyperparams["O"]);
    
    List si_0_list = hyperparams["si_0"];
    hyperparams_si_0.resize(numSamples);
    for (int i = 0; i < numSamples; ++i) {
      hyperparams_si_0[i] = as<VectorXd>(si_0_list[i]);
    }
    
    if (P > 0) {
      hyperparams_S_A = as<double>(hyperparams["S_A"]);
    }
    if (L > 0) {
      hyperparams_S_R = as<double>(hyperparams["S_R"]);
      hyperparams_S_Qi_mean = as<double>(hyperparams["S_Qi_mean"]);
      hyperparams_S_oi_mean = as<double>(hyperparams["S_oi_mean"]);
    }
    
    int max_cells = 0;
    for (int i = 0; i < numSamples; ++i) {
      max_cells += aux_Ni(i);
    }
    workspace_Y_res = MatrixXd(J, max_cells);
    workspace_B_concat = MatrixXd(K, max_cells);
    workspace_ZD = MatrixXd(J, K);
    workspace_Yp = MatrixXd(J, N);
    
    if (P > 0) {
      precompute_A_inverse();
    }
    
    tracking_dsi.resize(numSamples);
    tracking_doi.resize(numSamples);
    tracking_dBi.resize(numSamples);
    tracking_dQi.resize(numSamples);
    tracking_dRk.resize(K);

    // Pre-allocate tracking storage (CRITICAL FIX for memory growth)
    tracking_prev_Z.resize(J, K);
    tracking_prev_o.resize(J);
    if (P > 0) tracking_prev_A.resize(P, K);
    if (L > 0) tracking_prev_Ro.resize(J, L);

    tracking_prev_Bi.resize(numSamples);
    tracking_prev_Qi.resize(numSamples);
    tracking_prev_si.resize(numSamples);
    tracking_prev_oi.resize(numSamples);

    for (int i = 0; i < numSamples; ++i) {
      tracking_prev_Bi[i].resize(K, aux_Ni(i));
      tracking_prev_Qi[i].resize(J, K);
      tracking_prev_si[i].resize(aux_Ni(i));
      tracking_prev_oi[i].resize(J);
    }

    if (L > 0) {
      tracking_prev_Rk.resize(K);
      for (int k = 0; k < K; ++k) {
        tracking_prev_Rk[k].resize(J, L);
      }
    }
  }
  
  List initialize(bool multimodal = false) {
    
    if (verbose >= 1 && !multimodal) {
      Rcout << "Initializing latent variables..." << std::endl;
    }
    
    {
      FunctionTimer timer("solve_oi_all", verbose);
      solve_oi_all_init();
    }
    
    {
      FunctionTimer timer("compute_residual_Yp", verbose);
      compute_residual_Yp();
    }
    
    {
      FunctionTimer timer("perform_rsvd", verbose);
      perform_rsvd();
    }
    
    {
      FunctionTimer timer("initialize_from_svd", verbose);
      initialize_Z_U_S_from_svd();
      initialize_Qi_from_Z();
      solve_initial_Bi_all();
    }
    
    if (!multimodal) {
      {
        FunctionTimer timer("normalize_B", verbose);
        normalize_B(false);
      }
      {
        FunctionTimer timer("update_ZDBi", verbose);
        update_ZDBi();
      }
      {
        FunctionTimer timer("update_QiDBi", verbose);
        update_QiDBi();
      }
    }
    
    is_initialized = true;
    
    if (verbose >= 1 && !multimodal) {
      Rcout << "Initialization complete." << std::endl;
    }
    
    return wrap_results();
  }
  
  List optimize(int iterations, int track_interval = 5) {
    
    if (!is_initialized) {
      stop("Model must be initialized before optimization. Call initialize() first.");
    }
    
#ifdef _OPENMP
    int original_threads = omp_get_max_threads();
    if (num_threads > 0) {
      omp_set_num_threads(num_threads);
    }
    int actual_threads = omp_get_max_threads();
    if (verbose >= 1) {
      Rcout << "OpenMP enabled: Using " << actual_threads << " threads" << std::endl;
    }
#else
    if (verbose >= 1 && num_threads > 0) {
      Rcout << "OpenMP not available - using single-threaded execution" << std::endl;
    }
#endif
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    tracking_sigma2.resize(iterations);
    tracking_dZ.resize(iterations, NA_REAL);
    tracking_dA.resize(iterations, NA_REAL);
    tracking_do.resize(iterations, NA_REAL);
    tracking_dRo.resize(iterations, NA_REAL);
    
    for (int i = 0; i < numSamples; ++i) {
      tracking_dsi[i].resize(iterations, NA_REAL);
      tracking_doi[i].resize(iterations, NA_REAL);
      tracking_dBi[i].resize(iterations, NA_REAL);
      tracking_dQi[i].resize(iterations, NA_REAL);
    }
    for (int k = 0; k < K; ++k) {
      tracking_dRk[k].resize(iterations, NA_REAL);
    }
    
    if (verbose >= 1) {
      Rcout << "Starting block coordinate descent optimization..." << std::endl;
      Rcout << "Parameters: " << iterations << " iterations, " << numSamples 
            << " samples, " << J << " genes, " << K << " factors" << std::endl;
    }
    
    for (int ite = 0; ite < iterations; ++ite) {

      // Use pre-allocated tracking storage (CRITICAL FIX for memory growth)
      bool do_tracking = (ite % track_interval == 0) || (track_interval == 1);

      if (do_tracking) {
        // Reuse existing allocations - no new memory allocated!
        tracking_prev_Z.noalias() = params_Z;
        tracking_prev_o.noalias() = params_o;
        if (P > 0) tracking_prev_A.noalias() = params_A;
        if (L > 0) tracking_prev_Ro.noalias() = params_Ro;

        for (int i = 0; i < numSamples; ++i) {
          tracking_prev_Bi[i].noalias() = params_Bi[i];
          tracking_prev_Qi[i].noalias() = params_Qi[i];
          tracking_prev_si[i].noalias() = params_si[i];
          tracking_prev_oi[i].noalias() = params_oi[i];
        }

        if (L > 0) {
          for (int k = 0; k < K; ++k) {
            tracking_prev_Rk[k].noalias() = params_Rk[k];
          }
        }
      }
      
      if (verbose >= 1) {
        Rcout << "  Iteration " << (ite + 1) << "/" << iterations << std::endl;
      }
      
      {
        FunctionTimer timer("solve_Bi_all", verbose);
        solve_Bi_all();
      }
      
      {
        FunctionTimer timer("normalize_B", verbose);
        normalize_B(false);
      }
      
      {
        FunctionTimer timer("update_QiDBi (1)", verbose);
        update_QiDBi();
      }
      
      if (orthoZ) {
        FunctionTimer timer("solve_Z_orthogonal", verbose);
        solve_Z_orthogonal();
      } else {
        FunctionTimer timer("solve_Z_regular", verbose);
        solve_Z_regular();
      }
      
      {
        FunctionTimer timer("update_ZDBi", verbose);
        update_ZDBi();
      }
      
      {
        FunctionTimer timer("solve_Qi_all", verbose);
        solve_Qi_all();
      }
      
      {
        FunctionTimer timer("update_QiDBi (2)", verbose);
        update_QiDBi();
      }
      
      {
        FunctionTimer timer("solve_oi_all", verbose);
        solve_oi_all();
      }
      
      if (!is_si_fixed) {
        FunctionTimer timer("solve_si_all", verbose);
        solve_si_all();
      }
      
      {
        FunctionTimer timer("solve_o", verbose);
        solve_o();
      }
      
      {
        FunctionTimer timer("solve_Yi_all", verbose);
        solve_Yi_all();
      }
      
      {
        FunctionTimer timer("solve_sigma2", verbose);
        params_sigma2 = solve_sigma2();
      }
      tracking_sigma2[ite] = params_sigma2;
      
      if (verbose >= 1) {
        Rcout << "    sigma2: " << params_sigma2 << std::endl;
      }
      
      {
        FunctionTimer timer("update_hyperpriors", verbose);
        update_hyperpriors();
      }
      
      if (verbose >= 2) {
        Rcout << "    Mean(o): " << hyperparams_o_0[0] 
              << "; Var(o): " << hyperparams_S_o 
              << "; Var(si): " << hyperparams_S_si << std::endl;
      }
      
      if (P > 0) {
        FunctionTimer timer("solve_A", verbose);
        solve_A();
      }
      
      if (L > 0) {
        {
          FunctionTimer timer("solve_Rk_all", verbose);
          solve_Rk_all();
        }
        {
          FunctionTimer timer("solve_Ro", verbose);
          solve_Ro();
        }
      }
      
      if (do_tracking) {
        FunctionTimer timer("calculate_tracking", verbose);
        calculate_tracking(ite);
      }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
#ifdef _OPENMP
    omp_set_num_threads(original_threads);
#endif
    
    total_iterations += iterations;
    
    if (verbose >= 1) {
      Rcout << "Optimization completed in " << duration.count() << " ms" << std::endl;
      Rcout << "Final sigma2: " << params_sigma2 << std::endl;
    }
    
    return List::create(
      Named("params") = wrap_params(),
      Named("aux") = wrap_aux(),
      Named("hyperparams") = wrap_hyperparams(),
      Named("tracking") = wrap_tracking(),
      Named("convergence") = List::create(
        Named("iterations") = iterations,
        Named("total_iterations") = total_iterations,
        Named("final_sigma2") = params_sigma2,
        Named("total_time_ms") = duration.count(),
        Named("success") = true
      )
    );
  }
  
  List train(int iterations, int track_interval = 5, bool multimodal = false) {
    
    if (verbose >= 1) {
      Rcout << "Starting GEDI training (initialization + optimization)..." << std::endl;
    }
    
    List init_result = initialize(multimodal);
    List opt_result = optimize(iterations, track_interval);
    
    if (verbose >= 1) {
      Rcout << "Training complete." << std::endl;
    }
    
    return opt_result;
  }

private:
  
  void solve_oi_all_init() {
    for (int i = 0; i < numSamples; ++i) {
      double lambda = 1.0 / hyperparams_S_oi(i);
      double Ni_val = aux_Ni(i);

      MatrixXd residual = target_Yi[i] - aux_ZDBi[i] - aux_QiDBi[i];
      residual.rowwise() -= params_si[i].transpose();
      residual.colwise() -= params_o;

      VectorXd row_sums = residual * aux_Ni_vec[i];

      if (L == 0) {
        params_oi[i] = row_sums / (Ni_val + lambda);
      } else {
        params_oi[i] = (row_sums + lambda * aux_oi_hat[i]) / (Ni_val + lambda);
      }
    }
  }

  void compute_residual_Yp() {
    int current_col_offset = 0;
    for (int i = 0; i < numSamples; ++i) {
      const int Ni = aux_Ni(i);
      const VectorXd gene_offset = params_o + params_oi[i];

      MatrixXd sample_residual = target_Yi[i];
      sample_residual.colwise() -= gene_offset;
      sample_residual.rowwise() -= params_si[i].transpose();

      workspace_Yp.block(0, current_col_offset, J, Ni) = sample_residual;
      current_col_offset += Ni;
    }
  }

  void perform_rsvd() {
    List svd_result = rsvd_internal(workspace_Yp, hyperparams_O, K, 2, K, K);
    
    VectorXd svd_d = as<VectorXd>(svd_result["d"]);
    MatrixXd svd_u = as<MatrixXd>(svd_result["u"]);
    
    params_Z = svd_u * svd_d.head(K).asDiagonal() / 2.0;
    params_U = svd_u;
    params_S = svd_d.head(K) / 2.0;
  }

  List rsvd_internal(const MatrixXd& A_const, const MatrixXd& O, int k, 
                     int q = 2, int nu = -1, int nv = -1) {
    
    if (nu < 0) nu = k;
    if (nv < 0) nv = k;
    nu = std::min(k, nu);
    nv = std::min(k, nv);

    MatrixXd A = A_const;
    int m = A.rows();
    int n = A.cols();
    bool flipped = false;

    if (m < n) {
      A = A.transpose();
      std::swap(m, n);
      std::swap(nu, nv);
      flipped = true;
    }

    MatrixXd Y = A * O;

    MatrixXd Q(m, Y.cols());
    if (q > 0) {
      for (int i = 0; i < q; ++i) {
        HouseholderQR<MatrixXd> qr_y(Y);
        Q = qr_y.householderQ() * MatrixXd::Identity(m, Y.cols());
        
        MatrixXd Z = A.transpose() * Q;
        
        HouseholderQR<MatrixXd> qr_z(Z);
        MatrixXd Q_z = qr_z.householderQ() * MatrixXd::Identity(n, Z.cols());
        
        Y = A * Q_z;
      }
    }

    HouseholderQR<MatrixXd> qr_final(Y);
    Q = qr_final.householderQ() * MatrixXd::Identity(m, Y.cols());

    MatrixXd B = Q.transpose() * A;

    JacobiSVD<MatrixXd> svd(B, ComputeThinU | ComputeThinV);
    VectorXd singular_values = svd.singularValues();
    MatrixXd matrix_u = svd.matrixU();
    MatrixXd matrix_v = svd.matrixV();

    MatrixXd final_u;
    MatrixXd final_v;

    if (nu > 0) {
      final_u = Q * matrix_u.leftCols(nu);
    }

    if (nv > 0) {
      final_v = matrix_v.leftCols(nv);
    }

    if (flipped) {
      std::swap(final_u, final_v);
    }

    return List::create(
      Named("d") = singular_values.head(k),
      Named("u") = final_u,
      Named("v") = final_v
    );
  }

  void initialize_Z_U_S_from_svd() {
  }

  void initialize_Qi_from_Z() {
    for (int i = 0; i < numSamples; ++i) {
      params_Qi[i] = params_Z;
    }
  }

  void solve_initial_Bi_all() {
    for (int i = 0; i < numSamples; ++i) {
      MatrixXd Wi(J, K);
      for (int k = 0; k < K; ++k) {
        Wi.col(k) = (params_Z.col(k) + params_Qi[i].col(k)) * params_D(k);
      }

      VectorXd gene_offset = params_o + params_oi[i];
      MatrixXd residual = target_Yi[i];
      residual.colwise() -= gene_offset;
      residual.rowwise() -= params_si[i].transpose();

      MatrixXd WtW = Wi.transpose() * Wi;
      MatrixXd WtY = Wi.transpose() * residual;
      LDLT<MatrixXd> solver(WtW);
      params_Bi[i] = solver.solve(aux_diag_K * WtY);
    }
  }
  
  void solve_Bi_all() {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(numSamples > 1 && num_threads > 1)
#endif
    for (int i = 0; i < numSamples; ++i) {
      // const int Ni = aux_Ni(i);  // Unused local copy

      MatrixXd Wi(J, K);
      for (int k = 0; k < K; ++k) {
        Wi.col(k) = (params_Z.col(k) + params_Qi[i].col(k)) * params_D(k);
      }
      
      VectorXd gene_offset = params_o + params_oi[i];
      MatrixXd residual = target_Yi[i];
      residual.colwise() -= gene_offset;
      residual.rowwise() -= params_si[i].transpose();
      
      MatrixXd WtW = Wi.transpose() * Wi;
      MatrixXd WtY = Wi.transpose() * residual;
      LDLT<MatrixXd> solver(WtW);
      params_Bi[i] = solver.solve(aux_diag_K * WtY);
    }
  }
  
  void normalize_B(bool multimodal = false) {
    if (!multimodal && mode == "Bsphere") {
      for (int i = 0; i < numSamples; ++i) {
        VectorXd col_norms = params_Bi[i].colwise().norm();
        for (int j = 0; j < params_Bi[i].cols(); ++j) {
          if (col_norms(j) > 1e-10) {
            double scale = 1.0 / col_norms(j);
            params_Bi[i].col(j) *= scale;
          }
        }
      }
    }
    
    VectorXd total_row_norms = VectorXd::Zero(K);
    for (int i = 0; i < numSamples; ++i) {
      total_row_norms += params_Bi[i].rowwise().squaredNorm();
    }
    total_row_norms = total_row_norms.cwiseSqrt();
    
    VectorXd row_scaling = VectorXd::Ones(K);
    for (int k = 0; k < K; ++k) {
      if (total_row_norms(k) > 1e-10) {
        row_scaling(k) = 1.0 / total_row_norms(k);
      }
    }
    
    if (!multimodal && mode == "Bl2") {
      for (int i = 0; i < numSamples; ++i) {
        for (int k = 0; k < K; ++k) {
          params_Bi[i].row(k) *= row_scaling(k);
        }
      }
    } else if (adjustD) {
      params_D = row_scaling;
    }
  }
  
  void update_QiDBi() {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(numSamples > 2 && num_threads > 1)
#endif
    for (int i = 0; i < numSamples; ++i) {
      MatrixXd QiD = params_Qi[i] * params_D.asDiagonal();
      aux_QiDBi[i] = QiD * params_Bi[i];
    }
  }
  
  void update_ZDBi() {
    workspace_ZD = params_Z * params_D.asDiagonal();
    
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(numSamples > 2 && num_threads > 1)
#endif
    for (int i = 0; i < numSamples; ++i) {
      aux_ZDBi[i] = workspace_ZD * params_Bi[i];
    }
  }
  
  void solve_Z_orthogonal() {
    int total_cells = 0;
    for (int i = 0; i < numSamples; ++i) {
      total_cells += aux_Ni(i);
    }
    
    std::vector<int> col_offsets(numSamples + 1);
    col_offsets[0] = 0;
    for (int i = 0; i < numSamples; ++i) {
      col_offsets[i + 1] = col_offsets[i] + aux_Ni(i);
    }
    
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(numSamples > 2 && num_threads > 1)
#endif
    for (int i = 0; i < numSamples; ++i) {
      const int Ni = aux_Ni(i);
      const int offset = col_offsets[i];
      VectorXd gene_offset = params_o + params_oi[i];
      
      MatrixXd Yi_res = target_Yi[i] - aux_QiDBi[i];
      Yi_res.colwise() -= gene_offset;
      Yi_res.rowwise() -= params_si[i].transpose();
      workspace_Y_res.block(0, offset, J, Ni) = Yi_res;
      
      workspace_B_concat.block(0, offset, K, Ni) = params_Bi[i];
    }
    
    MatrixXd B = params_D.asDiagonal() * workspace_B_concat.leftCols(total_cells);
    
    double lambda = 1.0 / hyperparams_S_Z;
    
    MatrixXd Y_final, B_final;
    
    if (P > 0) {
      double sqrt_lambda = std::sqrt(lambda);
      MatrixXd CA_term = sqrt_lambda * aux_C * params_A;
      MatrixXd diagK_term = sqrt_lambda * aux_diag_K;
      
      Y_final = MatrixXd(J, total_cells + K);
      Y_final.leftCols(total_cells) = workspace_Y_res.leftCols(total_cells);
      Y_final.rightCols(K) = CA_term;
      
      B_final = MatrixXd(K, total_cells + K);
      B_final.leftCols(total_cells) = B;
      B_final.rightCols(K) = diagK_term;
    } else {
      Y_final = workspace_Y_res.leftCols(total_cells);
      B_final = B;
    }
    
    MatrixXd SB = params_S.asDiagonal() * B_final;
    MatrixXd Y_SB_T = Y_final * SB.transpose();
    
    JacobiSVD<MatrixXd> svd(Y_SB_T, ComputeThinU | ComputeThinV);
    params_U = svd.matrixU() * svd.matrixV().transpose();
    
    MatrixXd YBt = Y_final * B_final.transpose();
    MatrixXd YBtU = YBt.cwiseProduct(params_U);
    
    params_S = YBtU.transpose() * aux_J_vec / (1.0 + lambda);
    params_Z = params_U * params_S.asDiagonal();
  }
  
  void solve_Z_regular() {
    int total_cells = 0;
    for (int i = 0; i < numSamples; ++i) {
      total_cells += aux_Ni(i);
    }
    
    std::vector<int> col_offsets(numSamples + 1);
    col_offsets[0] = 0;
    for (int i = 0; i < numSamples; ++i) {
      col_offsets[i + 1] = col_offsets[i] + aux_Ni(i);
    }
    
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(numSamples > 2 && num_threads > 1)
#endif
    for (int i = 0; i < numSamples; ++i) {
      const int Ni = aux_Ni(i);
      const int offset = col_offsets[i];
      VectorXd gene_offset = params_o + params_oi[i];
      
      MatrixXd Yi_res = target_Yi[i] - aux_QiDBi[i];
      Yi_res.colwise() -= gene_offset;
      Yi_res.rowwise() -= params_si[i].transpose();
      workspace_Y_res.block(0, offset, J, Ni) = Yi_res;
      
      workspace_B_concat.block(0, offset, K, Ni) = params_Bi[i];
    }
    
    double lambda = 1.0 / hyperparams_S_Z;
    
    MatrixXd B = workspace_B_concat.leftCols(total_cells);
    MatrixXd DB = params_D.asDiagonal() * B;
    MatrixXd Y_DB_T = workspace_Y_res.leftCols(total_cells) * DB.transpose();
    
    if (P == 0) {
      MatrixXd gram = DB * DB.transpose() + lambda * aux_diag_K;
      LDLT<MatrixXd> solver(gram);
      params_Z = Y_DB_T * solver.solve(aux_diag_K);
    } else {
      MatrixXd CA_term = lambda * aux_C * params_A;
      MatrixXd gram = DB * DB.transpose() + lambda * aux_diag_K;
      LDLT<MatrixXd> solver(gram);
      params_Z = (Y_DB_T + CA_term) * solver.solve(aux_diag_K);
    }
  }
  
  void solve_Qi_all() {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(numSamples > 1 && num_threads > 1)
#endif
    for (int i = 0; i < numSamples; ++i) {
      double lambda = 1.0 / hyperparams_S_Qi(i);
      
      VectorXd gene_offset = params_o + params_oi[i];
      MatrixXd residual = target_Yi[i] - aux_ZDBi[i];
      residual.colwise() -= gene_offset;
      residual.rowwise() -= params_si[i].transpose();
      
      MatrixXd DBi = params_D.asDiagonal() * params_Bi[i];
      
      MatrixXd gram = DBi * DBi.transpose() + lambda * aux_diag_K;
      LDLT<MatrixXd> solver(gram);
      
      if (L == 0) {
        params_Qi[i] = (residual * DBi.transpose()) * solver.solve(aux_diag_K);
      } else {
        MatrixXd lhs = residual * DBi.transpose() + lambda * aux_Qi_hat[i];
        params_Qi[i] = lhs * solver.solve(aux_diag_K);
      }
    }
  }
  
  void solve_oi_all() {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(numSamples > 2 && num_threads > 1)
#endif
    for (int i = 0; i < numSamples; ++i) {
      double lambda = 1.0 / hyperparams_S_oi(i);
      double Ni_val = aux_Ni(i);
      
      MatrixXd residual = target_Yi[i] - aux_ZDBi[i] - aux_QiDBi[i];
      residual.rowwise() -= params_si[i].transpose();
      residual.colwise() -= params_o;
      
      VectorXd row_sums = residual * aux_Ni_vec[i];
      
      if (L == 0) {
        params_oi[i] = row_sums / (Ni_val + lambda);
      } else {
        params_oi[i] = (row_sums + lambda * aux_oi_hat[i]) / (Ni_val + lambda);
      }
    }
  }
  
  void solve_si_all() {
    double lambda = 1.0 / hyperparams_S_si;
    
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(numSamples > 2 && num_threads > 1)
#endif
    for (int i = 0; i < numSamples; ++i) {
      VectorXd gene_offset = params_o + params_oi[i];
      MatrixXd residual = target_Yi[i] - aux_ZDBi[i] - aux_QiDBi[i];
      residual.colwise() -= gene_offset;
      
      VectorXd col_sums = residual.transpose() * aux_J_vec;
      params_si[i] = (col_sums + lambda * hyperparams_si_0[i]) / (J + lambda);
    }
  }
  
  void solve_o() {
    VectorXd o_sum = VectorXd::Zero(J);
    
    for (int i = 0; i < numSamples; ++i) {
      MatrixXd residual = target_Yi[i] - aux_ZDBi[i] - aux_QiDBi[i];
      residual.rowwise() -= params_si[i].transpose();
      residual.colwise() -= params_oi[i];
      
      o_sum += residual * aux_Ni_vec[i];
    }
    
    double lambda = 1.0 / hyperparams_S_o;
    params_o = (o_sum + hyperparams_o_0 * lambda) / (N + lambda);
  }
  
  void solve_Yi_all() {
    if (obs_type == "Y") {
      return;
    }
    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(numSamples > 1 && num_threads > 1)
#endif
    for (int i = 0; i < numSamples; ++i) {
      if (obs_type == "M") {
        ArrayXXd Yi = target_Yi[i].array();
        ArrayXXd ZDBi = aux_ZDBi[i].array();
        ArrayXXd QiDBi = aux_QiDBi[i].array();
        ArrayXd si = params_si[i].array();
        ArrayXd oi = params_oi[i].array();
        
        ArrayXXd Yi_hat = ((ZDBi + QiDBi).rowwise() + si.transpose()).colwise() + 
          (params_o.array() + oi);
        ArrayXXd logMi = (ArrayXXd(1.0 * target_Mi[i]) + 1e-10).log();
        ArrayXXd solution = Yi;
        ArrayXXd alpha = Yi - Yi_hat - ArrayXXd(params_sigma2 * target_Mi[i]);
        
        const int total_elements = solution.size();
        
        for (int j = 0; j < total_elements; ++j) {
          double exp_Yi, f, fp, fpp, upper, lower;
          
          if (Yi(j) > 0) {
            exp_Yi = std::exp(-Yi(j));
            fpp = params_sigma2;
            fp = params_sigma2 + exp_Yi;
            f = params_sigma2 + exp_Yi * alpha(j);
          } else {
            exp_Yi = std::exp(Yi(j));
            fpp = params_sigma2 * exp_Yi;
            fp = fpp + 1.0;
            f = fpp + alpha(j);
          }
          
          solution(j) -= 2.0 * f * fp / (2.0 * fp * fp - f * fpp);
          
          if (logMi(j) < Yi_hat(j)) {
            upper = Yi_hat(j);
            lower = logMi(j);
          } else {
            upper = logMi(j);
            lower = Yi_hat(j);
          }
          
          if (solution(j) < lower) {
            solution(j) = lower;
          } else if (solution(j) > upper) {
            solution(j) = upper;
          }
        }
        
        target_Yi[i] = solution.matrix();
        
      } else if (obs_type == "M_paired") {
        ArrayXXd Yi = target_Yi[i].array();
        ArrayXXd ZDBi = aux_ZDBi[i].array();
        ArrayXXd QiDBi = aux_QiDBi[i].array();
        ArrayXd si = params_si[i].array();
        ArrayXd oi = params_oi[i].array();
        
        ArrayXXd Yi_hat = ((ZDBi + QiDBi).rowwise() + si.transpose()).colwise() + 
          (params_o.array() + oi);
        ArrayXXd solution = Yi;
        ArrayXXd alpha = Yi_hat - ArrayXXd(params_sigma2 * target_M2i[i]);
        ArrayXXd beta = Yi_hat + ArrayXXd(params_sigma2 * target_M1i[i]);
        
        const int total_elements = solution.size();
        
        for (int j = 0; j < total_elements; ++j) {
          double exp_Yi, f, fp, fpp;
          
          if (Yi(j) > 0) {
            exp_Yi = std::exp(-Yi(j));
            f = Yi(j) - alpha(j) + (Yi(j) - beta(j)) * exp_Yi;
            fp = 1.0 + exp_Yi * (beta(j) - Yi(j) + 1.0);
            fpp = exp_Yi * (Yi(j) - beta(j) - 2.0);
          } else {
            exp_Yi = std::exp(Yi(j));
            f = exp_Yi * (Yi(j) - alpha(j)) + Yi(j) - beta(j);
            fp = exp_Yi + beta(j) - Yi(j) + 1.0;
            fpp = Yi(j) - beta(j) - 2.0;
          }
          
          solution(j) -= 2.0 * f * fp / (2.0 * fp * fp - f * fpp);
          
          if (solution(j) < alpha(j)) {
            solution(j) = alpha(j);
          } else if (solution(j) > beta(j)) {
            solution(j) = beta(j);
          }
        }
        
        target_Yi[i] = solution.matrix();
        
      } else if (obs_type == "X") {
        MatrixXd Yi_pred = (aux_ZDBi[i] + aux_QiDBi[i]).rowwise() + params_si[i].transpose();
        Yi_pred = Yi_pred.colwise() + (params_o + params_oi[i]);
        
        const int J_local = Yi_pred.rows();
        const int Ni = Yi_pred.cols();
        
        for (int j = 0; j < J_local; ++j) {
          for (int n = 0; n < Ni; ++n) {
            double X_val = target_Xi[i](j, n);
            
            if (X_val == 1.0) {
              double Yi_val = Yi_pred(j, n);
              double dnorm_val = -0.5 * std::log(2.0 * M_PI) - 0.5 * Yi_val * Yi_val;
              double pnorm_val = std::log(0.5 * (1.0 + std::erf(Yi_val / std::sqrt(2.0))));
              target_Yi[i](j, n) = Yi_val + std::exp(dnorm_val - pnorm_val);
              
            } else if (X_val == 0.0) {
              double Yi_val = Yi_pred(j, n);
              double dnorm_val = -0.5 * std::log(2.0 * M_PI) - 0.5 * Yi_val * Yi_val;
              double pnorm_val = std::log(0.5 * (1.0 + std::erf(-Yi_val / std::sqrt(2.0))));
              target_Yi[i](j, n) = Yi_val - std::exp(dnorm_val - pnorm_val);
            } else {
              target_Yi[i](j, n) = Yi_pred(j, n);
            }
          }
        }
      }
    }
  }
  
  double solve_sigma2() {
    int denom_N = J * (N + numSamples) +
      numSamples * J * K +
      (K + 1) * J * L +
      K * (J + P);
    
    double S = 0.0;
    
    if (P == 0) {
      S += params_Z.squaredNorm() / hyperparams_S_Z;
    } else {
      MatrixXd Z_residual = params_Z - aux_C * params_A;
      S += Z_residual.squaredNorm() / hyperparams_S_Z;
      S += params_A.squaredNorm() / hyperparams_S_Z / hyperparams_S_A;
    }
    
    for (int i = 0; i < numSamples; ++i) {
      if (L == 0) {
        S += params_oi[i].squaredNorm() / hyperparams_S_oi(i);
      } else {
        VectorXd oi_diff = params_oi[i] - aux_oi_hat[i];
        S += oi_diff.squaredNorm() / hyperparams_S_oi(i);
      }
      
      if (L == 0) {
        S += params_Qi[i].squaredNorm() / hyperparams_S_Qi(i);
      } else {
        MatrixXd Qi_diff = params_Qi[i] - aux_Qi_hat[i];
        S += Qi_diff.squaredNorm() / hyperparams_S_Qi(i);
      }
      
      VectorXd gene_offset = params_o + params_oi[i];
      MatrixXd residual = target_Yi[i] - aux_ZDBi[i] - aux_QiDBi[i];
      residual.colwise() -= gene_offset;
      residual.rowwise() -= params_si[i].transpose();
      
      if (obs_type == "Y") {
        S += residual.squaredNorm();
      } else if (obs_type == "M") {
        S += residual.squaredNorm();
        ArrayXXd Yi_exp = target_Yi[i].array().exp();
        ArrayXXd denominator = Yi_exp + (1.0 / params_sigma2);
        S += (1.0 / denominator).sum();
      } else if (obs_type == "M_paired") {
        S += residual.squaredNorm();
        ArrayXXd expYi = (-target_Yi[i].array().abs()).exp();
        ArrayXXd M_sum = ArrayXXd(target_M1i[i] + target_M2i[i]);
        ArrayXXd one_plus_expYi = 1.0 + expYi;
        ArrayXXd denominator = M_sum * expYi / one_plus_expYi.square() + 
          (1.0 / params_sigma2);
        S += (1.0 / denominator).sum();
      }
    }
    
    if (L > 0) {
      for (int k = 0; k < K; ++k) {
        S += params_Rk[k].squaredNorm() / hyperparams_S_Qi_mean / hyperparams_S_R;
      }
      S += params_Ro.squaredNorm() / hyperparams_S_oi_mean / hyperparams_S_R;
    }
    
    if (obs_type != "X") {
      return S / denom_N;
    } else {
      return 1.0;
    }
  }
  
  void update_hyperpriors() {
    double o_mean = params_o.mean();
    hyperparams_o_0.setConstant(o_mean);
    VectorXd o_diff = params_o - hyperparams_o_0;
    hyperparams_S_o = (o_diff.squaredNorm() + 2.0) / params_sigma2 / (J + 2.0);
    
    double si_sum = 0.0;
    for (int i = 0; i < numSamples; ++i) {
      VectorXd si_diff = params_si[i] - hyperparams_si_0[i];
      si_sum += si_diff.squaredNorm();
    }
    hyperparams_S_si = (si_sum / params_sigma2 + 2.0) / (N + 2.0);
  }
  
  void solve_A() {
    MatrixXd CtZ = aux_C.transpose() * params_Z;
    params_A = workspace_CtC_inv * CtZ;
  }
  
  void solve_Rk_all() {
    MatrixXd Hp(L, numSamples);
    for (int i = 0; i < numSamples; ++i) {
      double scaling = 1.0 / std::sqrt(hyperparams_S_Qi(i));
      Hp.col(i) = aux_H.col(i) * scaling;
    }
    
    double lambda = 1.0 / hyperparams_S_Qi_mean / hyperparams_S_R;
    MatrixXd HpHp_T = Hp * Hp.transpose() + lambda * MatrixXd::Identity(L, L);
    LDLT<MatrixXd> solver(HpHp_T);
    workspace_HpHp_inv = solver.solve(MatrixXd::Identity(L, L));
    
    for (int k = 0; k < K; ++k) {
      MatrixXd Qk(J, numSamples);
      for (int i = 0; i < numSamples; ++i) {
        double scaling = 1.0 / std::sqrt(hyperparams_S_Qi(i));
        Qk.col(i) = params_Qi[i].col(k) * scaling;
      }
      
      MatrixXd QkHp_T = Qk * Hp.transpose();
      params_Rk[k] = QkHp_T * workspace_HpHp_inv;
      
      MatrixXd RkH = params_Rk[k] * aux_H;
      for (int i = 0; i < numSamples; ++i) {
        aux_Qi_hat[i].col(k) = RkH.col(i);
      }
    }
  }
  
  void solve_Ro() {
    MatrixXd O(J, numSamples);
    MatrixXd Hp(L, numSamples);
    
    for (int i = 0; i < numSamples; ++i) {
      double scaling = 1.0 / std::sqrt(hyperparams_S_oi(i));
      O.col(i) = params_oi[i] * scaling;
      Hp.col(i) = aux_H.col(i) * scaling;
    }
    
    double lambda = 1.0 / hyperparams_S_oi_mean / hyperparams_S_R;
    MatrixXd HpHp_T = Hp * Hp.transpose() + lambda * MatrixXd::Identity(L, L);
    LDLT<MatrixXd> solver(HpHp_T);
    MatrixXd HpHp_inv = solver.solve(MatrixXd::Identity(L, L));
    
    MatrixXd OHp_T = O * Hp.transpose();
    params_Ro = OHp_T * HpHp_inv;
    
    MatrixXd RoH = params_Ro * aux_H;
    for (int i = 0; i < numSamples; ++i) {
      aux_oi_hat[i] = RoH.col(i);
    }
  }
  
  void precompute_A_inverse() {
    double lambda = 1.0 / hyperparams_S_A;
    MatrixXd CtC = aux_C.transpose() * aux_C + lambda * MatrixXd::Identity(P, P);
    LDLT<MatrixXd> solver(CtC);
    workspace_CtC_inv = solver.solve(MatrixXd::Identity(P, P));
  }
  
  // Simplified tracking (uses pre-allocated storage)
  void calculate_tracking(int ite) {

    tracking_do[ite] = std::sqrt((params_o - tracking_prev_o).squaredNorm() / params_o.size());
    tracking_dZ[ite] = std::sqrt((params_Z - tracking_prev_Z).squaredNorm() / params_Z.size());

    if (P > 0) {
      tracking_dA[ite] = std::sqrt((params_A - tracking_prev_A).squaredNorm() / params_A.size());
    }

    if (L > 0) {
      tracking_dRo[ite] = std::sqrt((params_Ro - tracking_prev_Ro).squaredNorm() / params_Ro.size());
    }

    for (int i = 0; i < numSamples; ++i) {
      tracking_dsi[i][ite] = std::sqrt((params_si[i] - tracking_prev_si[i]).squaredNorm() /
        params_si[i].size());
      tracking_doi[i][ite] = std::sqrt((params_oi[i] - tracking_prev_oi[i]).squaredNorm() /
        params_oi[i].size());
      tracking_dBi[i][ite] = std::sqrt((params_Bi[i] - tracking_prev_Bi[i]).squaredNorm() /
        params_Bi[i].size());
      tracking_dQi[i][ite] = std::sqrt((params_Qi[i] - tracking_prev_Qi[i]).squaredNorm() /
        params_Qi[i].size());
    }

    if (L > 0) {
      for (int k = 0; k < K; ++k) {
        tracking_dRk[k][ite] = std::sqrt((params_Rk[k] - tracking_prev_Rk[k]).squaredNorm() /
          params_Rk[k].size());
      }
    }
  }
  
  List wrap_results() {
    return List::create(
      Named("params") = wrap_params(),
      Named("aux") = wrap_aux(),
      Named("hyperparams") = wrap_hyperparams()
      // NO target! Saves massive memory
    );
  }
  
  List wrap_params() {
    List Bi_list(numSamples), Qi_list(numSamples);
    List si_list(numSamples), oi_list(numSamples);
    
    for (int i = 0; i < numSamples; ++i) {
      Bi_list[i] = params_Bi[i];
      Qi_list[i] = params_Qi[i];
      si_list[i] = params_si[i];
      oi_list[i] = params_oi[i];
    }
    
    List result = List::create(
      Named("Bi") = Bi_list,
      Named("Qi") = Qi_list,
      Named("si") = si_list,
      Named("oi") = oi_list,
      Named("o") = params_o,
      Named("Z") = params_Z,
      Named("U") = params_U,
      Named("S") = params_S,
      Named("D") = params_D,
      Named("sigma2") = params_sigma2
    );
    
    if (P > 0) {
      result["A"] = params_A;
    }
    
    if (L > 0) {
      List Rk_list(K);
      for (int k = 0; k < K; ++k) {
        Rk_list[k] = params_Rk[k];
      }
      result["Rk"] = Rk_list;
      result["Ro"] = params_Ro;
    }
    
    return result;
  }
  
  List wrap_aux() {
    List result = List::create(
      Named("J") = J,
      Named("N") = N,
      Named("K") = K,
      Named("P") = P,
      Named("L") = L,
      Named("numSamples") = numSamples,
      Named("obs.type") = obs_type,
      Named("mode") = mode,
      Named("orthoZ") = orthoZ,
      Named("adjustD") = adjustD,
      Named("is_si_fixed") = is_si_fixed
    );
    
    return result;
  }
  
  List wrap_hyperparams() {
    List si_0_list(numSamples);
    for (int i = 0; i < numSamples; ++i) {
      si_0_list[i] = hyperparams_si_0[i];
    }
    
    List result = List::create(
      Named("S_Qi") = hyperparams_S_Qi,
      Named("S_oi") = hyperparams_S_oi,
      Named("S_Z") = hyperparams_S_Z,
      Named("S_o") = hyperparams_S_o,
      Named("S_si") = hyperparams_S_si,
      Named("o_0") = hyperparams_o_0,
      Named("si_0") = si_0_list
    );
    
    if (P > 0) {
      result["S_A"] = hyperparams_S_A;
    }
    
    if (L > 0) {
      result["S_R"] = hyperparams_S_R;
      result["S_Qi_mean"] = hyperparams_S_Qi_mean;
      result["S_oi_mean"] = hyperparams_S_oi_mean;
    }
    
    return result;
  }
  
  List wrap_tracking() {
    List dsi_list(numSamples), doi_list(numSamples);
    List dBi_list(numSamples), dQi_list(numSamples);
    
    for (int i = 0; i < numSamples; ++i) {
      dsi_list[i] = wrap(tracking_dsi[i]);
      doi_list[i] = wrap(tracking_doi[i]);
      dBi_list[i] = wrap(tracking_dBi[i]);
      dQi_list[i] = wrap(tracking_dQi[i]);
    }
    
    List result = List::create(
      Named("sigma2") = tracking_sigma2,
      Named("dZ") = tracking_dZ,
      Named("do") = tracking_do,
      Named("dsi") = dsi_list,
      Named("doi") = doi_list,
      Named("dBi") = dBi_list,
      Named("dQi") = dQi_list
    );
    
    if (P > 0) {
      result["dA"] = tracking_dA;
    }
    
    if (L > 0) {
      result["dRo"] = tracking_dRo;
      List dRk_list(K);
      for (int k = 0; k < K; ++k) {
        dRk_list[k] = wrap(tracking_dRk[k]);
      }
      result["dRk"] = dRk_list;
    }
    
    return result;
  }
};

// ============================================================================
// Rcpp Exports
// ============================================================================

//' Create New GEDI Model Object (Internal)
//'
//' Constructs the C++ GEDI model object and returns an external pointer.
//' This function is called internally by the R6 GEDI class during model setup
//' and should not be called directly by users.
//'
//' @param params List containing initialized parameter matrices:
//'   \itemize{
//'     \item Bi: List of cell projection matrices (K x Ni) for each sample
//'     \item Qi: List of sample-specific metagene matrices (J x K)
//'     \item si, oi: Cell and sample-specific offset vectors
//'     \item o: Global gene offset vector (length J)
//'     \item Z: Shared metagene matrix (J x K)
//'     \item U, S, D: SVD components and scaling factors
//'     \item sigma2: Initial variance estimate
//'     \item A: Pathway-latent factor connection matrix (if C provided)
//'     \item Rk, Ro: Sample covariate effect matrices (if H provided)
//'   }
//' @param aux List containing auxiliary variables and dimensions:
//'   \itemize{
//'     \item J, N, K, P, L: Dimensions (genes, cells, factors, pathways, covariates)
//'     \item numSamples: Number of samples
//'     \item obs.type: Observation type ("M", "M_paired", "Y", or "X")
//'     \item mode: Normalization mode ("Bl2" or "Bsphere")
//'     \item orthoZ, adjustD, is_si_fixed: Boolean flags
//'     \item C, H: Prior matrices (gene pathways, sample covariates)
//'     \item diag_K, J_vec, Ni_vec, Ni: Dimension vectors
//'     \item ZDBi, QiDBi: Precomputed product matrices
//'   }
//' @param target List containing target data matrices:
//'   \itemize{
//'     \item Yi: Log-transformed expression by sample (computed by C++ if empty)
//'     \item Mi: Raw count matrix by sample (for obs.type = "M")
//'     \item M1i, M2i: Paired count matrices (for obs.type = "M_paired")
//'     \item Xi: Binary indicator matrix (for obs.type = "X")
//'   }
//' @param hyperparams List containing hyperparameters for regularization:
//'   \itemize{
//'     \item S_Qi, S_oi, S_Z, S_A, S_R: Shrinkage parameters
//'     \item S_Qi_mean, S_oi_mean, S_si, S_o: Mean shrinkage parameters
//'     \item o_0, si_0: Prior mean values
//'     \item O: Random matrix for rSVD initialization
//'   }
//' @param verbose Integer verbosity level (0=silent, 1=info, 2=debug, 4=timing)
//' @param num_threads Number of OpenMP threads (0=auto, uses all available)
//'
//' @return External pointer (SEXP) to C++ GEDI object
//'
//' @details
//' This function creates a stateful C++ GEDI object that stores all model parameters
//' and data. The C++ object performs all heavy computation, while the R6 wrapper
//' remains lightweight (~1 KB). Memory-efficient design: if Yi is not provided,
//' C++ computes it from raw counts M, eliminating duplicate storage in R.
//'
//' The returned external pointer is managed by R's garbage collector and will
//' automatically free the C++ object when no longer referenced.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP GEDI_new(List params, List aux, List target, List hyperparams,
              int verbose = 1, int num_threads = 0) {
  GEDI* model = new GEDI(params, aux, target, hyperparams, verbose, num_threads);
  XPtr<GEDI> ptr(model, true);
  return ptr;
}

//' Initialize Latent Variables (Internal)
//'
//' Initializes latent variables Z, Qi, and Bi using randomized SVD on residual
//' expression data. This is the first step in GEDI model fitting and must be
//' called before optimization.
//'
//' @param model_ptr External pointer to C++ GEDI object (from GEDI_new)
//' @param multimodal Logical, if TRUE skips normalization steps for multi-modal
//'   integration (used internally when combining multiple data modalities)
//'
//' @return List containing initialized model state:
//'   \itemize{
//'     \item params: Updated parameter matrices (Z, Qi, Bi, etc.)
//'     \item aux: Auxiliary variables and dimensions
//'     \item hyperparams: Hyperparameters
//'   }
//'
//' @details
//' Initialization procedure:
//' \enumerate{
//'   \item Solve for initial oi (sample-specific offsets)
//'   \item Compute residual Yp after removing o, oi, si effects
//'   \item Perform randomized SVD on Yp to obtain initial Z
//'   \item Initialize Qi = Z for all samples
//'   \item Solve for initial Bi given Z and Qi
//'   \item Normalize B matrices and update auxiliary products (ZDBi, QiDBi)
//' }
//'
//' Uses randomized SVD (rSVD) for computational efficiency with large matrices.
//' The multimodal flag is used internally for advanced integration workflows.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List GEDI_initialize(SEXP model_ptr, bool multimodal = false) {
  XPtr<GEDI> ptr(model_ptr);
  return ptr->initialize(multimodal);
}

//' Optimize GEDI Model via Block Coordinate Descent (Internal)
//'
//' Performs block coordinate descent optimization to fit the GEDI model.
//' Iteratively updates all model parameters (Bi, Z, Qi, oi, si, o, sigma2)
//' until convergence or maximum iterations reached.
//'
//' @param model_ptr External pointer to initialized C++ GEDI object
//' @param iterations Integer, number of optimization iterations to perform
//' @param track_interval Integer, interval for tracking convergence metrics
//'   (e.g., track_interval=5 tracks every 5th iteration; track_interval=1 tracks all)
//'
//' @return List containing optimized model results:
//'   \itemize{
//'     \item params: Optimized parameter matrices (Z, Bi, Qi, sigma2, etc.)
//'     \item aux: Auxiliary variables and dimensions
//'     \item hyperparams: Updated hyperparameters (o_0, S_o, S_si)
//'     \item tracking: Convergence metrics (dZ, dBi, dQi, sigma2 trajectory)
//'     \item convergence: Summary statistics (iterations, total_time_ms, final_sigma2)
//'   }
//'
//' @details
//' Block coordinate descent iteration:
//' \enumerate{
//'   \item Solve Bi for all samples (parallelized with OpenMP)
//'   \item Normalize B matrices
//'   \item Solve Z (orthogonal or regular, depending on orthoZ flag)
//'   \item Solve Qi for all samples (parallelized)
//'   \item Solve oi, si, o (offset parameters)
//'   \item Update Yi from raw counts M (for count data)
//'   \item Solve sigma2 (variance parameter)
//'   \item Update hyperpriors (o_0, S_o, S_si)
//'   \item Solve A, Rk, Ro (if priors C or H are provided)
//'   \item Track convergence metrics at specified intervals
//' }
//'
//' OpenMP parallelization is automatically enabled if available. Progress and
//' timing information printed based on verbosity level.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List GEDI_optimize(SEXP model_ptr, int iterations, int track_interval = 5) {
  XPtr<GEDI> ptr(model_ptr);
  return ptr->optimize(iterations, track_interval);
}

//' Train GEDI Model (Initialize + Optimize) (Internal)
//'
//' Convenience function that performs both initialization and optimization in a
//' single call. Equivalent to calling GEDI_initialize() followed by GEDI_optimize().
//'
//' @param model_ptr External pointer to C++ GEDI object (from GEDI_new)
//' @param iterations Integer, number of optimization iterations to perform
//' @param track_interval Integer, interval for tracking convergence metrics
//' @param multimodal Logical, if TRUE skips normalization during initialization
//'   (for multi-modal integration workflows)
//'
//' @return List containing trained model results (same format as GEDI_optimize):
//'   \itemize{
//'     \item params: Optimized parameter matrices
//'     \item aux: Auxiliary variables
//'     \item hyperparams: Updated hyperparameters
//'     \item tracking: Convergence metrics
//'     \item convergence: Summary statistics
//'   }
//'
//' @details
//' This is the recommended way to fit a GEDI model in a single step. Internally:
//' \enumerate{
//'   \item Calls initialize() to set up latent variables via rSVD
//'   \item Calls optimize() to perform block coordinate descent
//'   \item Returns the final optimized model state
//' }
//'
//' For more control over the fitting process, use GEDI_initialize() and
//' GEDI_optimize() separately. The multimodal parameter is for advanced use cases.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List GEDI_train(SEXP model_ptr, int iterations, int track_interval = 5, bool multimodal = false) {
  XPtr<GEDI> ptr(model_ptr);
  return ptr->train(iterations, track_interval, multimodal);
}
