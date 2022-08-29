#' Laplacian_matrix
#' @description Laplacian_matrix is a function which calculates the Laplacian matrix of a graph and which apply
#' a filter (sigma) on the eigenvalues of the Laplacian matrix in order to reduce high frequencies induced by the GraphNet penalty.
#' @return This function returns the Laplacian matrix and low frequencies Laplacian matrix associated with the graph of the chosen block.

#' @param graph is the graph associated with the block, the initial graph.
#'
#' @param sigma is the eigenvalues filter, between 0 and 2 (Eigenvalues of the Laplacian matrix normalised are between 0 and 2).
#' @export

#' @examples
#' library(igraph)
#' W <- matrix(rnorm(200, mean = 9, sd = 3), nrow=10)
#' g1 <- sample_grg(ncol(W),0.4)
#'
#' Laplacian_matrix(g1, 0.6)


Laplacian_matrix <- function(graph, sigma) {
  L <- igraph::graph.laplacian(graph, normalized =TRUE)

  p <- igraph::vcount(graph)
  evL <- eigen(L)
  D <- evL$values

  #Filter of eigenvalues
  Dlow <- diag(D)[(p-sum(D<sigma)+1):p,(p-sum(D<sigma)+1):p]
  Plow <- evL$vectors[, (p-sum(D<sigma)+1):p, drop = FALSE]
  lowfreq_matrix <- Plow%*%Dlow%*%t(Plow)
  dimnames(lowfreq_matrix) <- dimnames(L)

  return(list(Laplacian=L, filter_Laplacian=lowfreq_matrix))
}
