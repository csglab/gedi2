#' @importFrom utils packageVersion installed.packages

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
      "\nInstall with: gedi::install_optional_dependencies()\n"
    )
  }
}

# .onLoad <- function(libname, pkgname) {
#   # Load the DLL
#   library.dynam(pkgname, pkgname, libname)
# }

.onUnload <- function(libpath) {
  library.dynam.unload("gedi", libpath)
}


#' Install Optional GEDI Dependencies
#'
#' Convenience function to install optional dependencies for extended GEDI functionality.
#' This includes packages for H5 file reading, UMAP embeddings, and strict validation.
#'
#' @param which Character vector specifying which dependencies to install.
#'   Options: "h5" (hdf5r), "umap" (uwot), "validation" (digest), or "all" (default).
#' @param repos Character vector. The CRAN repository to use for installation.
#'   Default is the user's configured repository.
#'
#' @return NULL (invisibly)
#'
#' @examples
#' \dontrun{
#' # Install all optional dependencies
#' install_optional_dependencies()
#'
#' # Install only H5 support
#' install_optional_dependencies(which = "h5")
#'
#' # Install UMAP and validation support
#' install_optional_dependencies(which = c("umap", "validation"))
#' }
#'
#' @export
install_optional_dependencies <- function(which = "all", repos = getOption("repos")) {
  
  # Define package groups
  pkg_groups <- list(
    h5 = c("hdf5r"),
    umap = c("uwot"),
    validation = c("digest")
  )
  
  # Determine which packages to install
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
  
  message("Checking optional GEDI dependencies...")
  
  # Check which packages are not installed
  installed <- packages %in% rownames(utils::installed.packages())
  to_install <- packages[!installed]
  
  if (length(to_install) > 0) {
    message("Installing: ", paste(to_install, collapse = ", "))
    utils::install.packages(to_install, repos = repos)
    message("\n=== Installation complete! ===")
    
    # Verify installation
    newly_installed <- to_install %in% rownames(utils::installed.packages())
    if (all(newly_installed)) {
      message("All packages successfully installed.")
    } else {
      failed <- to_install[!newly_installed]
      warning("Failed to install: ", paste(failed, collapse = ", "), call. = FALSE)
    }
  } else {
    message("All requested packages are already installed.")
  }
  
  invisible(NULL)
}


#' Check GEDI Optional Dependencies
#'
#' Check which optional dependencies are installed and display their status.
#'
#' @return Named logical vector indicating which optional packages are installed
#'
#' @examples
#' \dontrun{
#' # Check dependency status
#' check_optional_dependencies()
#' }
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