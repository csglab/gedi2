#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n",
    "=======================================================================================================\n",
    "  Welcome to GEDI v2 — Developed by the Computational and Statistical Genomics Laboratory\n",
    "  at McGill University, Montreal. Distributed under the MIT License.\n",
    "  For documentation and examples, see the package vignette.\n",
    "-------------------------------------------------------------------------------------------------------\n",
    "  Bienvenue dans GEDI v2 — Developpe par le Laboratoire de genomique computationnelle et statistique\n",
    "  de l'Universite McGill, a Montreal. Distribue sous licence MIT.\n",
    "  Pour la documentation et des exemples, consultez la vignette du package.\n",
    "=======================================================================================================\n"
  )
}

.onLoad <- function(libname, pkgname) {
  # Check for required packages
  required_pkgs <- c("Rcpp", "RcppEigen", "Matrix", "R6")
  missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

  if (length(missing_pkgs)) {
    packageStartupMessage(
      "Installing missing dependencies: ",
      paste(missing_pkgs, collapse = ", ")
    )
    utils::install.packages(missing_pkgs, dependencies = TRUE)
  }

  # Load the DLL
  library.dynam(pkgname, pkgname, libname)
}

.onUnload <- function(libpath) {
  library.dynam.unload("gedi", libpath)
}
