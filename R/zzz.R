#' @importFrom utils packageVersion

.onAttach <- function(libname, pkgname) {
  # Get package version from DESCRIPTION
  version <- utils::packageVersion(pkgname)
  
  packageStartupMessage(
    "\n",
    "=======================================================================================================\n",
    "  Welcome to GEDI v", version, " - Developed by the Computational and Statistical Genomics Laboratory\n",
    "  at McGill University, Montreal. Distributed under the MIT License.\n",
    "  For documentation and examples, see the package vignette.\n",
    "-------------------------------------------------------------------------------------------------------\n",
    "  Bienvenue dans GEDI v", version, " - Developpe par le Laboratoire de genomique computationnelle\n",
    "  et statistique de l'Universite McGill, a Montreal. Distribue sous licence MIT.\n",
    "  Pour la documentation et des exemples, consultez la vignette du package.\n",
    "=======================================================================================================\n"
  )

  # Check for optional dependencies
  missing_pkgs <- character(0)
  
  # Check H5 support
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, "hdf5r (for H5AD/H5 file reading)")
  }
  
  # Check UMAP support
  if (!requireNamespace("uwot", quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, "uwot (for UMAP embeddings)")
  }
  
  # Check digest for strict validation
  if (!requireNamespace("digest", quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, "digest (for strict M validation)")
  }
  
  if (length(missing_pkgs) > 0) {
    packageStartupMessage(
      "\nOptional dependencies not installed:\n",
      paste0("  - ", missing_pkgs, collapse = "\n"), "\n",
      "\nInstall them manually if needed, e.g.:\n",
      "  install.packages(c('hdf5r', 'uwot', 'digest'))\n"
    )
  }
}

# .onLoad <- function(libname, pkgname) {
#   # Load the DLL
#   library.dynam(pkgname, pkgname, libname)
# }

.onUnload <- function(libpath) {
  library.dynam.unload("gedi2", libpath)
}


#' List Optional GEDI Dependencies
#'
#' Reports which optional packages are needed and provides install commands.
#' Does **not** install anything automatically.
#'
#' @param which Character vector specifying which dependency groups to query.
#'   Options: "h5" (hdf5r), "umap" (uwot), "validation" (digest), or "all" (default).
#' @param verbose Logical, whether to print messages (default: TRUE)
#'
#' @return A named logical vector indicating which packages are installed (invisibly).
#'
#' @examples
#' # Show which optional packages are missing
#' install_optional_dependencies()
#'
#' @export
install_optional_dependencies <- function(which = "all", verbose = TRUE) {

  # Define package groups
  pkg_groups <- list(
    h5 = c("hdf5r"),
    umap = c("uwot"),
    validation = c("digest")
  )

  # Determine which packages to check
  if ("all" %in% which) {
    packages <- unlist(pkg_groups, use.names = FALSE)
  } else {
    invalid <- setdiff(which, names(pkg_groups))
    if (length(invalid) > 0) {
      stop("Invalid dependency group(s): ", paste(invalid, collapse = ", "), "\n",
           "Valid options: ", paste(names(pkg_groups), collapse = ", "), ", all",
           call. = FALSE)
    }
    packages <- unlist(pkg_groups[which], use.names = FALSE)
  }

  installed <- vapply(packages, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  }, logical(1))
  to_install <- packages[!installed]

  if (length(to_install) > 0) {
    if (verbose) {
      message("The following optional packages are not installed:")
      message("  ", paste(to_install, collapse = ", "))
      message("\nTo install them, run:")
      message("  install.packages(c(", paste0("'", to_install, "'", collapse = ", "), "))")
    }
  } else {
    if (verbose) message("All requested optional packages are already installed.")
  }

  invisible(installed)
}


#' Check GEDI Optional Dependencies
#'
#' Check which optional dependencies are installed and display their status.
#'
#' @return Named logical vector indicating which optional packages are installed
#'
#' @examples
#' check_optional_dependencies()
#'
#' @export
check_optional_dependencies <- function() {
  
  deps <- c(
    "hdf5r" = "H5AD/H5 file reading",
    "uwot" = "UMAP embeddings",
    "digest" = "Strict M validation (cryptographic hash)"
  )
  
  status <- vapply(names(deps), function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  }, logical(1))
  
  message("=== GEDI Optional Dependencies ===\n")
  for (i in seq_along(deps)) {
    pkg <- names(deps)[i]
    desc <- deps[i]
    installed <- status[i]
    
    status_msg <- if (installed) "\u2713 Installed" else "\u2717 Not installed"
    message(sprintf("%-12s %-40s %s", pkg, paste0("(", desc, ")"), status_msg))
  }
  
  if (!all(status)) {
    message("\nTo install missing dependencies: install_optional_dependencies()")
  }
  
  invisible(status)
}