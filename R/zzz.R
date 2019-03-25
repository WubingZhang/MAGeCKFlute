.onLoad <- function(libname, pkgname) {
  disclaimer="MAGeCKFlute is an open source software package. Particullary,
users are required to formally cite the original MAGeCKFlute paper in publications
or products. For details, do citation(\"MAGeCKFlute\") within R.\n\n"
  suppressMessages(require("pathview"))
  packageStartupMessage(paste(strwrap(disclaimer, 80), collapse="\n"))
}

.onAttach <- function(libname, pkgname) {
  env <- as.environment(paste0("package:", pkgname))
  env[[".conflicts.OK"]] <- TRUE
  suppressMessages(require("pathview"))
}
