#' netSGCCA
#' @description netSGCCA function performs netSGCCA method which is a multiblock method for the study of the relationship
#' between blocks and allows the integration of an information network reflecting the interactions among the variables
#' within a given data block. It returns a list of weight vectors and the component associated with each block.
#' @return The function netSGCCA returns a list of weight vectors and the component associated with each block.
#'
#' @param data is a list with all the blocks.
#' @param design is the design matrix.

#' @param s is a list of the quantity of sparsity for each block.  It is the radius of the l1 norm ball for weight vector of the each block.
#' For blocks associated with graph, s must be strictly superior to 1. If a block is not associated with a graph, put the square root of his column number.

#' @param gamma is a regulatory parameter between the GraphNet penalties. It's a vector with the same size as the data list.

#' @param lambda is a regulatory parameter between the interaction between blocks and the GraphNet penalty.  0<= lambda <1

#' @param Laplacian_list is a list of the Laplacian matrices in the same order as the data list.
#' If one block doesn't have a graph, so it doesn't have a Laplacian matrix, put a square matrix of zeros of the same size as the block in question.

#' @param epsilon is a stop condition.
#' @param iter_max is the maximal number of iteration.
#' @export

#' @examples
#' library(igraph)
#' ## Data :
#' W <- matrix(rnorm(200, mean = 9, sd = 3), nrow=10)
#' g1 <- sample_grg(ncol(W),0.4)
#'
#' Q <-  matrix(rnorm(100, mean = 5, sd = 3), nrow=10)
#' g2 <-  sample_grg(ncol(Q),0.3)
#'
#' P <- matrix(rnorm(150, mean = 2, sd = 1), nrow=10)
#' g3 <-  sample_grg(ncol(P),0.6)
#'
#' R= matrix(rnorm(30, 0), nrow=10)
#'
#' ## Laplacian Matrix
#' L= Laplacian_matrix(g1, 0.3)$Laplacian
#' M = matrix(0, ncol=ncol(R), nrow=ncol(R))
#' N =Laplacian_matrix(g2, 0.3)$Laplacian
#' O =Laplacian_matrix(g3, 0.3)$Laplacian
#'
#' ## Parameters
#' data= list(W, R, Q, P)
#' Laplacian_list=list(L, M, N, O)
#'
#' gamma=rep(0.5, length(data))
#' design=matrix(c(0,1,1,1, 1,0,1,1, 1,1,0, 1, 1,1,1,0), ncol=length(data))
#' lambda=0.5
#' s= c(0.5*sqrt(ncol(W)),sqrt(ncol(R)), 2.5, sqrt(ncol(P)) )
#' graphs <- list(g1, NA, g2, g3)
#' iter_max=100
#' epsilon=1e-10
#'
#' netSGCCA(data, design, s, gamma, lambda, Laplacian_list, epsilon, iter_max)



netSGCCA <- function(data, design, s, gamma, lambda, Laplacian_list, epsilon, iter_max) {

  #Weight vectors initialisation
  mat <- lapply(data, FUN=svd)  #svd on all blocks
  a_new <- mapply(function(i) mat[[i]]$v[,1], 1:length(mat))

  a_old <- a_new
  grad <- list()

  for(iter in 1:iter_max)  {

    for (k in 1:length(a_new)){
      Zk <- 0

      for (i in 1:length(a_new)) {
          Zk <- Zk+ design[i,k]*t(data[[k]])%*%data[[i]]%*%a_new[[i]]}

      grad[[k]] <- (1-lambda)*Zk + lambda*gamma[[k]]*2*Laplacian_list[[k]]%*%a_new[[k]]
      a_new[[k]] <- projL1L2(grad[[k]], s[k])$x
    }

    if (lapply(mapply("-", a_new, a_old, SIMPLIFY = FALSE), FUN= normL2) < epsilon) {
      break}

    else{ a_old <- a_new }
  }

  comp <- mapply(function(data, a_new) data%*%a_new, data, a_new, SIMPLIFY=F)

  return(list(weight_vectors=a_new, components=comp))
}

