#' Suppress stdout and messages
#'
#' Evaluates an expression while suppressing both stdout and message output.
#'
#' @param expr Expression to evaluate.
#' @return The result of evaluating `expr`.
#' @export
quiet <- function(expr) {
  nc <- file(nullfile(), open = "wt")
  sink(nc, type = "output"); sink(nc, type = "message")
  on.exit({ sink(type = "message"); sink(type = "output"); close(nc) })
  force(expr)
}