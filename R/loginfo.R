#' Log message
#'
#' Write log message on console
#'
#' @docType methods
#' @name loginfo
#' @rdname loginfo
#' @aliases log
#'
#' @param msg a character, specifying log message.
#'
#' @author Binbin Wang
#'
#' @note
#' The source can be found by typing \code{MAGeCKFlute:::loginfo}
#' or \code{getMethod("loginfo")}, or
#' browsed on github at \url{https://github.com/WubingZhang/MAGeCKFlute/tree/master/R/loginfo.R}
#' Users should find it easy to customize this function.
#'
#' @examples
#' loginfo("Test the loginfo function ...")
#'
#' @export
#'

loginfo <- function(msg) {message(paste(Sys.time()), ' # ', msg)}
