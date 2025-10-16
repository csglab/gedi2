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
}

.onLoad <- function(libname, pkgname) {
  # Load the DLL
  library.dynam(pkgname, pkgname, libname)
}

.onUnload <- function(libpath) {
  library.dynam.unload("gedi", libpath)
}
