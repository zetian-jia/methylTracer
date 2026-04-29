.onAttach <- function(libname, pkgname) {
  
  packageStartupMessage("Welcome to methylTracer")
  packageStartupMessage(
    "[1] ISSUES: https://github.com/zetian-jia/methylTracer")
  
  invisible()
}

.onLoad <- function(libname, pkgname) {
  # invisible()
}
