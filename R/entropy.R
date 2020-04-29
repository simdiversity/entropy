#' Compute entropy
#'
#' @param D A D matrix
#' @param f Weigh array
#' @example
#' M <- matrix(seq(9), nrow = 3)
#' f <- rowSums(M) / sum(M)
#' D <- dist(f*M)
#' entropy(D,f)
#' @export
entropy <- function(D, f) {
  if (class(D) == "dist") {
    D <- as.matrix(D)
  } else {
    D <- as.matrix(stats::as.dist(D))
  }
  s <- 1.0 - D/sum(D) # n times n similarity matrix
  b <- s %*% f # banalities
  R <- -sum(f %*% log(b)) # 0.5077364 reduced entropy
  exp(R) # 1.661526 number of effective varieties
}
