
entropy <- function(dissimilarity, weights) {
  s = 1 - dissimilarity # n times n similarity matrix
  b = s %*% weights # banalities
  R = -sum(weights %*% log(b)) # 0.5077364 reduced entropy
  exp(R) # 1.661526 number of effective varieties
}
