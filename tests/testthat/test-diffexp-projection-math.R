# tests/testthat/test-diffexp-projection-math.R
# -----------------------------------------------------------------------------
# Numerical-correctness ("math invariant") tests for the differential
# projections wired in for issue #26:
#   * model$diffADB(contrast)               (differential pathway activity, P x N)
#   * plot_features(..., projection = "diffexp"/"diffadb")
#
# These go beyond "is it a ggplot": they validate that the values produced by
# the C++ kernel getDiffADB_cpp and the efficient plot_features paths match an
# independent pure-R reconstruction of the documented formulae.
#
# Core identities being checked (derivation, with Cmat = inputC %*% C.rotation):
#   ZDB     = Z %*% DB                         (shared manifold)
#   ADB     = C.rotation %*% A %*% DB          (pathway activity)
#   A       = solveA(Z)                        (solve_A is linear in Z)
#   diffExp = diffQ_Z %*% DB                   (differential expression, J x N)
#   diffADB = C.rotation %*% solveA(diffQ_Z) %*% DB
#           = ADB(Z + diffQ_Z) - ADB(Z)        (TRUE differential of ADB)
# where solveA(X) = (Cmat^T Cmat + (1/S_A) I)^{-1} Cmat^T X applies the same
# shrinkage solve_A() uses (gedi_core.cpp). diffQ_Z[, k] = Rk[[k]] %*% (H.rot %*% contrast).
# -----------------------------------------------------------------------------

make_diff_math_fixture <- function() {
  set.seed(202606)
  n_genes    <- 50
  n_cells    <- 80
  n_samples  <- 4
  n_pathways <- 6
  K          <- 4

  M <- Matrix::Matrix(
    matrix(stats::rpois(n_genes * n_cells, lambda = 5), n_genes, n_cells),
    sparse = TRUE
  )
  rownames(M) <- paste0("Gene_", seq_len(n_genes))
  colnames(M) <- paste0("Cell_", seq_len(n_cells))

  samples <- factor(rep(paste0("S", seq_len(n_samples)),
                        each = n_cells / n_samples))

  # 3 covariates, linearly dependent (cov3 = cov1 + cov2) -> rank-deficient H,
  # so L < num_covariates and we also exercise the H.rotation projection.
  H <- matrix(
    c(1, 1, 0, 0,
      0, 0, 1, 1,
      1, 1, 1, 1),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("batch", "condition", "combined"),
                    paste0("S", seq_len(n_samples)))
  )

  # C pathway membership (J x n_pathways), 0/1, each gene in >= 1 pathway.
  set.seed(11)
  C <- matrix(0L, nrow = n_genes, ncol = n_pathways,
              dimnames = list(rownames(M), paste0("PW", seq_len(n_pathways))))
  for (i in seq_len(n_genes)) {
    C[i, sample(seq_len(n_pathways), 1)] <- 1L
  }
  # guarantee every pathway has at least one gene
  for (p in seq_len(n_pathways)) if (sum(C[, p]) == 0) C[p, p] <- 1L

  model <- CreateGEDIObject(
    Samples = samples, M = M, K = K, H = H, C = C, verbose = 0
  )
  model$train(iterations = 20, track_interval = 5)

  list(model = model, H = H, C = C, num_cov = nrow(H),
       n_pathways = n_pathways, n_genes = n_genes, n_cells = n_cells, K = K)
}

# Pure-R helpers rebuilt from stored model state ------------------------------

.rebuild <- function(m) {
  Z   <- m$params$Z                      # J x K
  D   <- m$params$D                      # K
  Bi  <- m$params$Bi                     # list of K x Ni
  Rk  <- m$params$Rk                     # list of J x L
  A   <- m$params$A                      # P x K
  Cin <- as.matrix(m$priors$inputC)      # J x num_pathways
  Cr  <- as.matrix(m$priors$C.rotation)  # num_pathways x P
  Hr  <- as.matrix(m$priors$H.rotation)  # L x num_cov
  S_A <- m$hyperparams$S_A
  K   <- m$aux$K
  P   <- m$aux$P

  B   <- do.call(cbind, Bi)              # K x N
  DB  <- diag(D, nrow = K) %*% B         # K x N
  Cmat <- Cin %*% Cr                     # J x P  (= aux_C)
  CtC <- t(Cmat) %*% Cmat + (1 / S_A) * diag(P)   # P x P
  solveA <- function(X) solve(CtC, t(Cmat) %*% X) # P x ncol(X)

  list(Z = Z, D = D, DB = DB, Rk = Rk, A = A, Cr = Cr, Hr = Hr,
       S_A = S_A, K = K, P = P, Cmat = Cmat, solveA = solveA)
}

.diffQ_Z <- function(r, contrast) {
  H_c <- as.vector(r$Hr %*% contrast)                       # length L
  vapply(r$Rk, function(Rk_k) as.vector(Rk_k %*% H_c),
         numeric(nrow(r$Rk[[1]])))                          # J x K
}

# -----------------------------------------------------------------------------
test_that("solve_A identity: stored A == (Cmat^T Cmat + (1/S_A) I)^{-1} Cmat^T Z", {
  fx <- make_diff_math_fixture(); m <- fx$model
  r  <- .rebuild(m)
  A_from_Z <- r$solveA(r$Z)                                 # P x K
  expect_equal(unname(r$A), unname(A_from_Z), tolerance = 1e-8)
})

# -----------------------------------------------------------------------------
test_that("ADB == C.rotation %*% A %*% DB (stored pathway activity)", {
  fx <- make_diff_math_fixture(); m <- fx$model
  r  <- .rebuild(m)
  ADB_R  <- r$Cr %*% r$A %*% r$DB                           # num_pathways x N
  ADB_m  <- m$projections$ADB
  expect_equal(dim(ADB_m), c(fx$n_pathways, fx$n_cells))
  expect_equal(unname(as.matrix(ADB_m)), unname(ADB_R), tolerance = 1e-8)
})

# -----------------------------------------------------------------------------
test_that("diffADB == C.rotation %*% solveA(diffQ_Z) %*% DB (cross-impl check)", {
  fx <- make_diff_math_fixture(); m <- fx$model
  r  <- .rebuild(m)
  contrast <- c(0.7, -0.4, 0.2)

  dQZ      <- .diffQ_Z(r, contrast)                         # J x K
  diffADB_R <- r$Cr %*% r$solveA(dQZ) %*% r$DB              # num_pathways x N
  diffADB_m <- m$diffADB(contrast)

  expect_equal(unname(as.matrix(diffADB_m)), unname(diffADB_R), tolerance = 1e-7)
})

# -----------------------------------------------------------------------------
test_that("diffADB is the TRUE differential of ADB: ADB(Z+dQ) - ADB(Z) == diffADB", {
  fx <- make_diff_math_fixture(); m <- fx$model
  r  <- .rebuild(m)
  contrast <- c(0.7, -0.4, 0.2)

  dQZ <- .diffQ_Z(r, contrast)                              # J x K
  ADB_base <- r$Cr %*% r$A %*% r$DB                         # = model$projections$ADB
  ADB_pert <- r$Cr %*% r$solveA(r$Z + dQZ) %*% r$DB         # perturb gene loadings
  diff_via_perturb <- ADB_pert - ADB_base

  expect_equal(unname(as.matrix(m$diffADB(contrast))),
               unname(diff_via_perturb), tolerance = 1e-7)
})

# -----------------------------------------------------------------------------
test_that("diffADB carries the solve_A shrinkage: diffADB == gamma * M * diffExp", {
  # Regression guard: an implementation that forgot the solve_A shrinkage would
  # produce diffADB = M %*% diffExp (too large by 1/gamma = 1 + 1/S_A).
  fx <- make_diff_math_fixture(); m <- fx$model
  r  <- .rebuild(m)
  contrast <- c(0.5, 0.3, -0.2)

  # Cmat is orthonormal (SVD basis) => CtC = (1 + 1/S_A) I => solveA = gamma * Cmat^T
  gamma <- 1 / (1 + 1 / r$S_A)                              # = S_A / (S_A + 1)
  Mop   <- r$Cr %*% t(r$Cmat)                               # num_pathways x J
  diffExp_noO <- m$diffExp(contrast, include_O = FALSE)     # J x N

  diffADB_scaled <- gamma * (Mop %*% as.matrix(diffExp_noO))
  expect_equal(unname(as.matrix(m$diffADB(contrast))),
               unname(diffADB_scaled), tolerance = 1e-7)

  # And explicitly confirm the shrinkage is NOT unity by default (A_shrinkage=1 => gamma=0.5)
  expect_lt(gamma, 0.9)
})

# -----------------------------------------------------------------------------
test_that("diffADB is linear in the contrast", {
  fx <- make_diff_math_fixture(); m <- fx$model
  c1 <- c(1, 0, 0); c2 <- c(0, 1, 0)

  d1  <- m$diffADB(c1)
  d2  <- m$diffADB(c2)
  d12 <- m$diffADB(c1 + c2)
  d2x <- m$diffADB(2 * c1)

  expect_equal(unname(as.matrix(d12)), unname(as.matrix(d1 + d2)), tolerance = 1e-8)
  expect_equal(unname(as.matrix(d2x)), unname(as.matrix(2 * d1)),  tolerance = 1e-8)
})

# -----------------------------------------------------------------------------
test_that("diffExp == diffQ_Z %*% DB and relates to ZDB by swapping Z -> diffQ_Z", {
  fx <- make_diff_math_fixture(); m <- fx$model
  r  <- .rebuild(m)
  contrast <- c(0.5, -0.3, 0.7)

  dQZ <- .diffQ_Z(r, contrast)
  diffExp_R <- dQZ %*% r$DB                                 # J x N
  expect_equal(unname(as.matrix(m$diffExp(contrast, include_O = FALSE))),
               unname(diffExp_R), tolerance = 1e-8)

  # ZDB = Z %*% DB (same DB map applied to Z instead of diffQ_Z)
  ZDB_R <- r$Z %*% r$DB
  expect_equal(unname(as.matrix(m$projections$ZDB)), unname(ZDB_R), tolerance = 1e-8)
})

# -----------------------------------------------------------------------------
test_that("plot_features(projection='diffexp') plotted values == full diffExp rows", {
  fx <- make_diff_math_fixture(); m <- fx$model
  contrast <- c(0.5, -0.3, 0.7)
  genes <- m$metadata$geneIDs[c(2, 7, 23)]
  emb <- cbind(seq_len(fx$n_cells), rev(seq_len(fx$n_cells)))  # N x 2, cell order

  full <- m$diffExp(contrast, include_O = FALSE)              # J x N
  p <- plot_features(m, features = genes, embedding = emb,
                     projection = "diffexp", contrast = contrast,
                     randomize = FALSE)
  expect_s3_class(p, "ggplot")
  pd <- p$data
  for (g in genes) {
    vals <- pd$Value[as.character(pd$Feature) == g]           # N values, cell order
    expect_equal(unname(vals), unname(as.vector(full[g, ])), tolerance = 1e-7)
  }

  # include_O path matches the offset-augmented full matrix
  full_O <- m$diffExp(contrast, include_O = TRUE)
  pO <- plot_features(m, features = genes, embedding = emb,
                      projection = "diffexp", contrast = contrast,
                      include_O = TRUE, randomize = FALSE)
  pdO <- pO$data
  for (g in genes) {
    vals <- pdO$Value[as.character(pdO$Feature) == g]
    expect_equal(unname(vals), unname(as.vector(full_O[g, ])), tolerance = 1e-7)
  }
})

# -----------------------------------------------------------------------------
test_that("plot_features(projection='diffadb') plotted values == diffADB rows", {
  fx <- make_diff_math_fixture(); m <- fx$model
  contrast <- c(0.5, -0.3, 0.7)
  pws <- c("PW1", "PW4", "PW6")
  emb <- cbind(seq_len(fx$n_cells), rev(seq_len(fx$n_cells)))

  full <- m$diffADB(contrast)                                 # num_pathways x N
  p <- plot_features(m, features = pws, embedding = emb,
                     projection = "diffadb", contrast = contrast,
                     randomize = FALSE)
  expect_s3_class(p, "ggplot")
  pd <- p$data
  for (pw in pws) {
    vals <- pd$Value[as.character(pd$Feature) == pw]
    expect_equal(unname(vals), unname(as.vector(full[pw, ])), tolerance = 1e-7)
  }
})
