#' This functions perform the projection of the weight vectors on the space defined
#' by the constraints in the function netSGCCA. It is used in the netSGCCA function.
#' @author Vincent Guillemot and Arnaud Gloaguen
#' @export

normL1 <- function(vec) {
  return(sum(abs(vec)))
}

normL2 <- function(vec) {
  return(sqrt(sum(vec**2)))
}

proxL1 <- function(vec, lambda) {
  return(sign(vec)*pmax(0, abs(vec) - lambda))
}


projL1L2 <- function(vec, rds) {
  norm2_x <- normL2(vec)
  if (norm2_x < .Machine$double.eps) {
    warning("Projecting a vector with a very small L2-norm!")
    return(list(x = 0 * vec, lambda = NA, k = NaN))
  }
  if (sum(abs(vec / norm2_x)) <= rds)
    return(list(
      x = vec / norm2_x,
      lambda = max(abs(vec)),
      k = NaN
    ))
  uneq <- vec != 0
  L <- sum(!uneq)
  p <- abs(vec[uneq])
  # Check if multiple maximum
  MAX <- max(p)
  bMAX <- p == MAX
  nMAX <- sum(bMAX)
  if (rds < sqrt(nMAX)) {
    warning("Minimum radius is: ", sqrt(nMAX))
    x_soft        <- rep(0, length(vec))
    x_soft[abs(vec) == MAX]  <- 1 / sqrt(nMAX)
    return(list(x = x_soft, lambda = MAX, k = NaN))
  } else if (rds == sqrt(nMAX)) {
    warning("radius is equal to sqrt(nMAX)")
    x_soft        <- rep(0, length(vec))
    x_soft[abs(vec) == MAX]  <- 1 / sqrt(nMAX)
    # x_soft[bMAX]  <- 1 / sqrt(nMAX)
    return(list(x = x_soft, lambda = MAX, k = NaN))
  }
  # 1. Take the absolute value of $\x$
  #   and sort its elements in decreasing order
  # to get $\widetilde{\x}$\;
  xtilde <- sort(abs(vec), decreasing = TRUE)
  psi_xtilde <- psi(xtilde, xtilde)
  # 2. Find i such that
  # $\psi(\widetilde{x}_{i+1})\leq c<\psi(\widetilde x_{i})$\;
  i <- max(which(psi_xtilde <= rds))
  # 3. Let $\displaystyle \delta = \frac{\normTwo{S(\widetilde{\x}	,
  #                                                 \widetilde x_i)}}{i}\left( c\sqrt{\frac{i-\psi(\widetilde x_i)^2}{i-c^2}}
  #                                                                            - \psi(\widetilde x_i)\right)$\;
  t1 <- normL2(proxL1(xtilde, xtilde[i])) / i
  t2 <- (i - psi_xtilde[i]^2) / (i - rds^2)
  t3 <- psi_xtilde[i]
  delta <- t1 * (rds * sqrt(t2) - t3)
  # 4. Compute $S(\x, \lambda)$ with $\lambda = \widetilde x_i - \delta$
  lambda <- xtilde[i] - delta
  x_soft <- sign(vec) * pmax(0, abs(vec) - lambda)
  return(list(
    x = x_soft / normL2(x_soft) ,
    lambda = lambda,
    k = NaN
  ))
}

phi <- Vectorize(function(x, lambda) {
  return(normL1(proxL1(vec = x, lambda = lambda)))
}, vectorize.args = "lambda")

psi <- Vectorize(function(x, lambda) {
  x_soft <- proxL1(vec = x, lambda = lambda)
  return(normL1(x_soft) / normL2(x_soft))
}, vectorize.args = "lambda")

