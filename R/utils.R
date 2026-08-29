#' Suppress stdout, messages, and warnings
#'
#' Evaluates an expression while suppressing both stdout and message output.
#' This is necessary because some dependencies (like RAS) emit verbose output
#' directly to stderr or through R's message system.
#'
#' @param expr Expression to evaluate.
#' @return The result of evaluating `expr`.
#' @export
quiet <- function(expr) {
  nc <- file(nullfile(), open = "wt")
  sink(nc, type = "output"); sink(nc, type = "message")
  on.exit({ sink(type = "message"); sink(type = "output"); close(nc) })
  suppressMessages(suppressWarnings(force(expr)))
}