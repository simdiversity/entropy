#' Compute effective entropy
#'
#' @param D A dissimilarity matrix
#' @param f Weigh array
#' @param Nloop Number of loops
#' @param Nfine finesse of powers
#' @param pa initial power
#' @param pb final power
#' @example
#' M <- matrix(seq(9), nrow = 3)
#' f <- rowSums(M) / sum(M)
#' D <- dist(f*M)
#' effective_entropy(D,f)
#' @export
effective_entropy <- function(D, f,
                              Nloop = 4000, Nfine=300,
                              pa = -4, pb = 3
                              ) {

  if (class(D) == "dist") {
    D <- as.matrix(D)
  } else {
    D <- as.matrix(stats::as.dist(D))
  }
  n <- nrow(D)

  power_selection <- pa + (seq(1:Nfine) - 1) * (pb - pa) / (Nfine - 1)
  Nsteps_beta <- length(power_selection)
  beta_rel <- c()

  counter_seq <- seq(length(power_selection))

  Delta <- as.numeric(0.5 * t(f) %*% D %*% f)
  Dif <- D %*% f - Delta

  iterations <- lapply(counter_seq, function(counter) {
    power <- power_selection[[counter]]
    beta_rel <- 10^power
    beta <- beta_rel * as.numeric(1 / Delta) # fixes the inverse temperature
    S <- as.matrix(exp(-beta * D)) # creates a D matrix
    b <- as.vector(S %*% f) #  banality

    ######### usual approach by memberships z_{ij}, with usual EM-iteration
    index <- which.min(D %*% f) # index de l'observation la plus proche
    indic_min <- as.numeric(1:n == index) # 0 partout, sauf 1 sur l'observation "centrale"

    rho <- as.vector(indic_min)
    E <- -sum(f * log(S %*% rho)) # effective entropy
    R <- -sum(f * log(S %*% f)) # reduced entropy
    HR <- 0 # group entropy
    Ty <- 1

    # Z <- diag(n) # initialisation du clustering soft, efficient for beta large
    #(and Niter, the number of iterations, can be small, convergence occurs rapidly: pure stayers)

    # initialisation alternative, BIEN meilleure pour les hautes temperatures
    Z <- matrix(0, n, n)
    Z[, index] <- 1
    eps10 <- 1e-20
    Z <- eps10 * matrix(1, n, n) + (1 - eps10) * Z

    # number of iterates (soft clustering)
    ones <- rep(1,n)
    res = soft_clustering(f, Z, S, Nloop)
    rho = res$rho
    Z = res$Z

    list(
      "power" = power,
      "rho" = rho,
      "E" = -sum(f *  log( S %*% rho)), # effective entropy
      "R" = -sum(f * log( S %*% f)), # reduced entropy
      "HR" = -sum(rho * log(rho + 10^(-13))), # group entropy
      "Ty" = sum(isFALSE(all.equal(rho,0))),
      "banalities" = b,
      "beta_rel" = beta_rel,
      "beta" = beta,
      "S" = S
    )
  })
  iterations
  result = list()
  for (name in c(
    "power", "E", "R", "HR", "Ty", "banalities", "beta", "beta_rel"
  )) {
    result[[name]] <- as.vector(do.call(rbind,
                                        lapply(iterations, function(el){el[[name]]})
    ))
  }
  result$S <-  lapply(iterations,  function(el){el$S})
  result$rho <- do.call(cbind, lapply(iterations, function(el){el$rho}))
  result
}
