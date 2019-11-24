.onLoad <- function(libname, pkgname) {
  disclaimer="MAGeCKFlute is an open source software package.
Users are required to formally cite the original MAGeCKFlute paper in publications."
  packageStartupMessage(paste(strwrap(disclaimer, 80), collapse="\n"))
}

.onAttach <- function(libname, pkgname) {
  env <- as.environment(paste0("package:", pkgname))
  env[[".conflicts.OK"]] <- TRUE
}
