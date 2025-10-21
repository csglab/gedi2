#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n",
    "=======================================================================================================\n",
    "  Welcome to GEDI - Developed by the Computational and Statistical Genomics Laboratory\n",
    "  at McGill University, Montreal. Distributed under the MIT License.\n",
    "  For documentation and examples, see the package vignette.\n",
    "-------------------------------------------------------------------------------------------------------\n",
    "  Bienvenue dans GEDI - Developpe par le Laboratoire de genomique computationnelle et statistique\n",
    "  de l'Universite McGill, a Montreal. Distribue sous licence MIT.\n",
    "  Pour la documentation et des exemples, consultez la vignette du package.\n",
    "=======================================================================================================\n"
  )

  # Check for optional H5 reading dependencies
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    packageStartupMessage(
      "\nNote: To use H5AD/H5 file reading functions (read_h5ad, read_h5), ",
      "please install the following package:\n",
      "  install.packages('hdf5r')\n",
      "Or use: gedi::install_h5_dependencies()\n"
    )
  }
}

.onLoad <- function(libname, pkgname) {
  # Load the DLL
  library.dynam(pkgname, pkgname, libname)
}

.onUnload <- function(libpath) {
  library.dynam.unload("gedi", libpath)
}


#' Install H5 file reading dependencies
#'
#' A convenience function to install the required dependencies for reading
#' H5AD and H5 files. This includes the hdf5r package for HDF5 file support.
#'
#' @param repos Character vector. The CRAN repository to use for installation.
#'   Default is the user's configured repository.
#'
#' @return NULL (invisibly)
#'
#' @examples
#' \dontrun{
#' # Install H5 reading dependencies
#' install_h5_dependencies()
#' }
#'
#' @export
install_h5_dependencies <- function(repos = getOption("repos")) {
  message("Installing H5 file reading dependencies...")

  # Packages needed for H5/H5AD support
  packages <- c("hdf5r")

  # Check which packages are not installed
  installed <- packages %in% rownames(installed.packages())
  to_install <- packages[!installed]

  if (length(to_install) > 0) {
    message("Installing: ", paste(to_install, collapse = ", "))
    install.packages(to_install, repos = repos)
    message("Installation complete!")
  } else {
    message("All required packages are already installed.")
  }

  invisible(NULL)
}
