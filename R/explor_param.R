#' Search of optimal parameters
#' @return a dataframe with a search grid and the results for the three criteria.
#'
#' @param data is a list with all the blocks.
#'
#' @param graphs is a list with graphs. If a block doesn't associated with a graph, put NA.
#'
#' @param lambda is a vector with different values of lambda between 0 and 1.
#'
#' @param sigma is a vector with different values of sigma, the eigenvalues filter,  0 <= sigma <= 2.
#'
#' @param s is a list of vectors of different values for the quantity of sparsity for each block. s > 1
#' For blocks associated with graph, s must be strictly greater than 1. For a block without a graph, s=sqrt(ncol(block)).
#' @export
#' @importFrom dplyr %>%
#' @examples
#' library(igraph)
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
#' # Example 1 (with a block without a graph)
#' data= list(W, R, Q, P)
#' graphs <- list(g1, NA, g2, g3)
#' lambda <- c(0, 0.5, 0.7, 0.9, 0.99)
#' gamma=rep(0.5, length(data))
#' design=matrix(c(0,1,1,1, 1,0,1,1, 1,1,0,1, 1,1,1,0), ncol=length(data))
#' sigma <- c(0.3, 0.5, 0.9)
#' s <- list(c(2, 0.5*sqrt(ncol(W))),sqrt(ncol(R)), c(4, 3.6, 0.5*sqrt(ncol(Q))), c(1.4,0.5*sqrt(ncol(P))))
#' iter_max=100
#' epsilon=1e-10
#'
#' explor_param(data, graphs, gamma, lambda, sigma, s)
#'
#'
#' # Example 2 (all blocks have an associated graph)
#' data= list(W, Q, P)
#' design=matrix(c(0,1,1,1,0,1,1,1,0), ncol=3)
#' iter_max=100
#' epsilon=1e-10
#' gamma=rep(0.3, length(data))
#' lambda <- c(0, 0.2, 0.5, 0.9)
#' sigma <- c(0.3, 0.5, 0.9, 1, 1.3)
#' s <- list(c(6, 0.5*sqrt(ncol(W)), sqrt(ncol(W))),c(1.5, 3, 0.5*sqrt(ncol(Q))), c(1.4,0.5*sqrt(ncol(P))))
#' graphs <- list(g1, g2, g3)
#'
#' explor_param(data, graphs, gamma, lambda, sigma, s)



explor_param <- function(data, graphs, gamma, lambda, sigma, s){

  #Build a dataframe in order to test each combinaison
  df <- expand.grid(lambda, sigma) %>% stats::setNames(c("lambda", "sigma"))
  kk=purrr::cross(s)
  ff=lapply(kk, rbind)
  df= merge(df, matrix(unlist(ff), ncol=length(data), byrow=T) )
  for (i in 1:length(s)){
    names(df)[2+i] <- paste0("s", i)
  }
  df[,(ncol(df)+1):(ncol(df)+(1+length(data)*2)*2)] <- NA  #For HF and LF

  no_graph <- c()


  for (i in 1:nrow(df)){

    #Construction of Laplacian matrices
    laplacianHF <- list()
    laplacianLF <- list()

    s <- as.numeric(df[i, 3:(length(data)+2)])

    for (j in 1:length(graphs)){

      # Laplacian Matrix for blocks associated with a graph
      if(is.na(graphs[[j]])!=TRUE){

        laplacianHF[[j]] <- Laplacian_matrix(graphs[[j]], df[i,2] )$Laplacian
        laplacianLF[[j]] <- Laplacian_matrix(graphs[[j]], df[i,2] )$filter_Laplacian
        }else  # Laplacian Matrix for block without a graph
          {
      laplacianHF[[j]] <- matrix(0, ncol=ncol(data[[j]]), nrow=ncol(data[[j]]))
      laplacianLF[[j]] <- matrix(0, ncol=ncol(data[[j]]), nrow=ncol(data[[j]]))
      no_graph <- c(no_graph, j)}
    }

    ### Sparsity ###
    netsgccaHF <- netSGCCA(data, design, s, gamma, df[i,1], laplacianHF, epsilon, iter_max )$weight_vectors
    netsgccaLF <- netSGCCA(data, design, s, gamma, df[i,1], laplacianLF, epsilon, iter_max )$weight_vectors

    df[i, (length(data)+3):(2*length(data)+2)] <- lapply(netsgccaHF, function(x) sum(x == 0))[1:length(data)]
    df[i, (length(data)*3+4):(length(data)*4+3)] <- lapply(netsgccaLF, function(x) sum(x == 0))[1:length(data)]


    for (k in 1:length(graphs)){
      names(df)[(length(data)+2)+k] <- paste0("sparsity_HF_block", k)}

    for (k in 1:length(graphs)){
      names(df)[length(data)*3+3+k] <- paste0("sparsity_LF_block", k)}

    ### Adequacy of data ###
    netsgccacompHF <- netSGCCA(data, design, s, gamma, df[i,1], laplacianHF, epsilon, iter_max )$components
    netsgccacompLF <- netSGCCA(data, design, s, gamma, df[i,1], laplacianLF, epsilon, iter_max )$components

    if(TRUE %in% is.na(graphs)){
      clini <- which(is.na(graphs)==TRUE)

      if (length(data[-clini])>1){#PLS regression if there is a block without a graph and several blocks associated with a graph
      regHF <- PLSR_SVD(sapply(lapply(netsgccacompHF[-clini], FUN=as.numeric), FUN=cbind), data[[clini]], 1)

      regLF <- PLSR_SVD(sapply(lapply(netsgccacompLF[-clini], FUN=as.numeric), FUN=cbind), data[[clini]], 1)

      df[i, 2*length(data)+3] <- regHF$R2y
      df[i, (length(data)*4+4)] <- regLF$R2y

      }else  #Case with just one block associated with a graph and a block without a graph
        {
        df[i, 2*length(data)+3] <- sum(cor(as.numeric(netsgccacompHF[[-clini]]), data[[clini]])**2)
        df[i, (length(data)*4+4)] <- sum(cor(as.numeric(netsgccacompLF[[-clini]]), data[[clini]])**2)
      }

    }else   #Case with no block without a graph
    {
      mean_varHF <- c()
      mean_varLF <- c()

      df[i, 2*length(data)+3] <- mean(sapply(mapply(function(a,b) cor(a,b)**2, data, lapply(netsgccacompHF, as.numeric)), mean))

      df[i, (length(data)*4+4)] <- mean(sapply(mapply(function(a,b) cor(a,b)**2, data, lapply(netsgccacompLF, as.numeric)), mean))

    }

    names(df)[2*length(data)+3] <- "Adequacy_data_HF"
    names(df)[(length(data)*4+4)] <- "Adequacy_data_LF"

    ### Graph Regularity ###
    if (length(no_graph)==0){ #No block without graph
      graphs_true=graphs
      netsgcca2HF <- netsgccaHF
      netsgcca2LF <- netsgccaLF

    }else  #if there is a block without graph
      {
      graphs_true <- graphs[-no_graph]
      netsgcca2HF <- netsgccaHF[-no_graph]
      netsgcca2LF <- netsgccaLF[-no_graph]}

    communities <- lapply(graphs_true, FUN=igraph::cluster_louvain)  #To find communities for each graph
    mem <- lapply(communities, FUN=igraph::membership)

    for (n in 1:length(graphs_true)){
      df[i, 2*length(data)+3+n] <- mean(tapply(netsgcca2HF[[n]], mem[[n]], var), na.rm=T)
      df[i,  (length(data)*4+4)+n] <- mean(tapply(netsgcca2LF[[n]], mem[[n]], var), na.rm=T)

      names(df)[2*length(data)+3+n] <- paste0("meanVarCommuHFBlock", n)
      names(df)[(length(data)*4+4)+n] <- paste0("meanVarCommuLFBlock", n)}

  }

  df=df[,apply(df, 2, function(x) !all(is.na(x)))]
  readr::write_csv(df,"search_parameters_Simulation.csv")

  return(df)

}

