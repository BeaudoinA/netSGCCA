#' This function is used to perform a PLS regression for the adequacy of data criterion.
#' Code from the package data4PCCAR on GitHub  https://github.com/HerveAbdi/data4PCCAR/blob/master/R/PLS_jack_svds_HA.R
#' @author Hervé Abdi
#' @export


PLSR_SVD <- function(X, Y,
                     nfactor,
                     inference = TRUE,
                     displayJack = TRUE){
  # First: An internal function to mimic MATLAB's repmat
  #___________________________________________________________________
  repmat <-  function (a, n, m){kronecker(matrix(1, n, m), a)}
  #___________________________________________________________________
  X <- as.data.frame(X) # fixes a strange error
  Y <- as.data.frame(Y)
  X <- as.matrix(X)
  # Make sure that we are dealing with matrices
  Y <- as.matrix(Y)
  obs.names  <- rownames(X)
  Xvar.names <- colnames(X)
  Yvar.names <- colnames(Y)
  X_ori <- X
  Y_ori <- Y

  maxfac <- qr(X)$rank-1;
  if(!exists("nfactor")){	nfactor < -maxfac}
  if (nfactor > maxfac){  nfactor <- maxfac}
  M_X <- apply(X_ori,2,mean)
  M_Y <- apply(Y_ori,2,mean)
  S_X <- apply(X_ori,2,sd)
  S_Y <- apply(Y_ori,2,sd)
  X <- scale(X_ori)
  Y <- scale(Y_ori)
  n <- nrow(Y)
  nn <- nrow(X)
  np <- ncol(X)
  nq <- ncol(Y)
  if (nn != n){
    stop("Incompatible # of rows for X and Y")
    return(NULL)
  }
  # Precision for Convergence.
  # not used in this version
  #	epsilon<-.Machine$double.eps
  # num of components kept
  # Initialization ----
  # The Y set
  U <- matrix(0,n,nfactor)
  C <- matrix(0,nq,nfactor)
  # The X set
  T <- matrix(0,n,nfactor)
  P <- matrix(0,np,nfactor)
  W <- matrix(0,np,nfactor)
  b <- matrix(0,1,nfactor)
  R2_X <- matrix(0,1,nfactor)
  R2_Y <- matrix(0,1,nfactor)
  RESS <- matrix(0,1,nfactor)
  PRESS <- matrix(0,1,nfactor)
  # Yhat4Press is a cube
  # of the jackknifed reconstitution of Y
  # Needed to compute PRESS
  Yjack4Press  <- array(0,dim = c(n,nq,nfactor))
  Yhat4Res_rec <- array(0,dim = c(n,nq,nfactor))
  Xres <- X
  Yres <- Y
  SS_X <- sum(X^2)
  SS_Y <- sum(Y^2)
  for (i in 1:nfactor){
    myplsr <- svd(t(Xres) %*% Yres)
    #		delta1 <- myplsr$d
    w <- myplsr$u[,1]
    c <- myplsr$v[,1]
    #tt <- Xres%*%w
    #t <- tt/sqrt(sum(tt^2))
    t <- normaliz(Xres %*% w)
    u <- Yres %*% c
    p <- t(Xres) %*% t
    b_l <- t(u) %*% t
    # Store in matrices
    b[i] <- b_l
    P[,i] <- p
    W[,i] <- w
    #T[,i]<- t
    T[,i] <- as.matrix(t)
    U[,i] <- u
    C[,i] <- c
    # deflation of X and Y
    Xres <- Xres -t %*% t(p)
    Yres <- Yres -
      repmat(b[i], n, nq)*(t %*% t(c))
    #R2_X[i] <- (t(t) %*% t) * (t(p)%*%p)/SS_X  #T has an unit norm
    R2_X[i] <- (t(p) %*% p) / SS_X
    R2_Y[i] <- (b[i]^2) / SS_Y
    #
    if (nq == 1){ # Problem with how R handles vectors vs matrices
      Yhat4Res_rec[,,i] <- (T[,1:i]*
                              repmat(t(b[1:i]),n,1)) %*%
        as.matrix(C[,1:i]) *
        repmat(t(S_Y),n,1)+
        repmat(t(M_Y),n,1)
    } else {
      Yhat4Res_rec[,,i] <-
        (T[,1:i] *
           repmat(t(b[1:i]),n,1)) %*%
        t(C[,1:i]) *
        repmat(t(S_Y),n,1) +
        repmat(t(M_Y),n,1)
    }
    RESS[i] <- sum((Y_ori - Yhat4Res_rec[,,i])^2)
  }
  X_rec<- T %*% t(P)
  Y_rec<- (T *
             repmat(b,n,1)) %*% t(C)
  #Bring back X and Y to their original units
  Xori_rec <- X_rec*
    repmat(t(S_X),n,1) +
    repmat(t(M_X),n,1)
  Yori_rec <- Y_rec*
    repmat(t(S_Y),n,1)+
    repmat(t(M_Y),n,1)
  #The Wstar weights give T = X*Wstar
  Wstar <- W %*% (solve(t(P) %*% W))
  # B_pls=Wstar*diag(b)*C'
  B_pls <- (Wstar *
              repmat(b,np,1)) %*% t(C)
  Bpls_star <- t(
    repmat(t(S_X^-1),nq,1)) * B_pls *
    repmat(t(S_Y), np, 1)
  Bpls_star <- rbind(-M_X %*% Bpls_star, Bpls_star)
  Bpls_star[1,] <- Bpls_star[1,] + M_Y
  #Y_pred <- cbind(1,X_ori)%*%Bpls_star
  if (isTRUE(inference)){
    # Now go to the jackknifed version (cross-validation)
    Yjack <- matrix(0,n,nq)
    print(paste0('Fixed Model Done. Start Jackknife for ', n, ' iterations.'))
    for (i in 1:n){# if inference ----
      if (displayJack){
        print(c('Jackniffing row #:',toString(i))) }
      X4j <- X_ori
      Y4j <- Y_ori
      X4j <- X4j[-i,]
      Y4j <- Y4j[-i,]
      leyhat <- PLS4jack(X4j,Y4j,X_ori[i,],nfactor)
      Yjack[i,] <- leyhat[nfactor,]

      for (l in 1:nfactor){
        Yjack4Press[i,,l] <- leyhat[l,]
      }
      r2_random <- matrix(0,1,nfactor)
      rv_random <- matrix(0,1,nfactor)
      for (l in 1:nfactor){
        PRESS[l] <- sum((Y_ori - Yjack4Press[,,l])^2)
        CORR <- corrcoef4mat(Y_ori, Yjack4Press[,,l])
        r2_random[l]<-CORR$r2_random
        rv_random[l]<-CORR$rv_random
      }
      Q2 <- 1 - PRESS[nfactor-1]/RESS[nfactor]
      # From Tenenhaus (1988), p83, 138 for Q2(l)
      Q2 <- c(PRESS[1] / (nq*(nn-1)),Q2)
    }
    # print('After Jackknife and Q')
  }
  if (isTRUE(inference)){	print('Jackknife done.') }
  #print('nfactor')
  #print(nfactor)

  Fact.names <- paste0("Factor ",1:nfactor)
  #print('dimT')
  #print(dim(T))
  rownames(T) <- obs.names
  colnames(T) <- Fact.names
  rownames(P) <- Xvar.names
  #print('dimP')
  #print(dim(P))
  colnames(P) <- Fact.names
  rownames(W) <- Xvar.names
  #print('dimW')
  #print(dim(W))
  colnames(W) <- Fact.names
  rownames(Wstar) <- Xvar.names
  #print('dim-Wstar')
  #print(dim(Wstar))
  colnames(Wstar) <- Fact.names
  rownames(U) <- obs.names
  colnames(U) <- Fact.names
  rownames(b) <- "b"
  colnames(b) <- Fact.names
  rownames(C) <- Yvar.names
  colnames(C) <- Fact.names
  rownames(B_pls) <- Xvar.names
  colnames(B_pls) <- Yvar.names
  rownames(Bpls_star) <- c('Intercept',Xvar.names)
  colnames(Bpls_star) <- Yvar.names
  rownames(Xori_rec)  <- obs.names
  colnames(Xori_rec)  <- Xvar.names
  rownames(Yori_rec)  <- obs.names
  colnames(Yori_rec)  <- Yvar.names

  rownames(R2_X) <- "R2x"
  colnames(R2_X) <- Fact.names
  rownames(R2_Y) <- "R2y"
  colnames(R2_Y) <- Fact.names
  rownames(RESS) <- "RESSy"
  colnames(RESS) <- Fact.names
  if (inference){ # if inference ----
    # if inference:
    #print('dim-PRESS')
    #print(dim(PRESS))
    rownames(Yjack)     <- obs.names
    colnames(Yjack)     <- Yvar.names
    rownames(PRESS) <- "PRESSy"
    #print(PRESS)
    colnames(PRESS) <- Fact.names
    #print(PRESS)
    Q2 <- t(as.matrix(Q2))
    #print('dim-Q2')
    #print(dim(Q2))
    rownames(Q2) <- "Q2"
    colnames(Q2) <- Fact.names[1:length(Q2)]
    rownames(r2_random) <- "Random-r2"
    colnames(r2_random) <- Fact.names
    #print('dim-rv_random')
    #print(dim(rv_random))
    rownames(rv_random) <- "Random-rv"
    colnames(rv_random) <- Fact.names
    dimnames(Yjack4Press)[[1]] <- obs.names
    dimnames(Yjack4Press)[[2]] <- Yvar.names
    dimnames(Yjack4Press)[[3]] <- paste0(1:nfactor,' Factor Solution')
    dimnames(Yhat4Res_rec)[[1]] <- obs.names
    dimnames(Yhat4Res_rec)[[2]] <- Yvar.names
    dimnames(Yhat4Res_rec)[[3]] <- paste0(1:nfactor,' Factor Solution')
    #
  } # end if inference ----

  # Return list 4 PLSR_SVD ----
  if (inference){
    return.list=structure(list(T = T, P = P, W = W,
                               Wstar = Wstar, U = U, B = b, C = C,
                               Bpls = B_pls, Bpls_star = Bpls_star,
                               Xhat = Xori_rec,
                               Yhat = Yori_rec,
                               R2x = R2_X, R2y = R2_Y,
                               RESSy = RESS,
                               Yhat4Ress = Yhat4Res_rec,
                               # Random part ----
                               Yjack = Yjack,
                               PRESSy = PRESS,
                               Q2 = Q2,
                               r2y_random = r2_random,
                               rv_random = rv_random,
                               Yhat4Press = Yjack4Press
    ),
    class = 'PLSR_SVD'
    )

  } else {
    return.list=structure(list(T = T, P = P,
                               W = W, Wstar = Wstar, U = U,
                               B = b, C = C,
                               Bpls = B_pls, Bpls_star = Bpls_star,
                               Xhat = Xori_rec,
                               Yhat = Yori_rec, #,
                               #
                               #Yjack=Yjack,
                               R2x = R2_X, R2y = R2_Y,
                               RESSy=RESS,
                               #PRESSy=PRESS,
                               #Q2=Q2,r2y_random=r2_random, rv_random=rv_random,
                               #Yhat4Press=Yjack4Press,
                               Yhat4Ress = Yhat4Res_rec
    ),
    class = 'PLSR_SVD'
    )
  }
  # Note the class "PLSR_SVD"
  # it is "linked" to the
  # function print.PLSR_SVD
  # that modifies the print function so that it will
  # print the description of the results
  # when the print function is used
  return(return.list)
} # End of function PLSR_SVD


normaliz <- function (F)
{
  repmat <- function(a, n, m) {
    kronecker(matrix(1, n, m), a)
  }
  v <- matrix(0, 1, ncol(F))
  for (j in 1:ncol(F)) {
    v[j] <- sqrt(sum(F[, j]^2))
  }
  f <- F/repmat(v, nrow(F), 1)
  return(f)
}

corrcoef4mat<-function (Y1, Y2)
{
  y1 <- c(Y1)
  y2 <- c(Y2)
  rv <- (y1 %*% y2)^2/((y1 %*% y1) * (y2 %*% y2))
  y1 <- y1 - mean(y1)
  y2 <- y2 - mean(y2)
  r2 <- (y1 %*% y2)^2/((y1 %*% y1) * (y2 %*% y2))
  return.list = list(rv_random = rv, r2_random = r2)
  return(return.list)
}


PLS4jack <- function (X, Y, xsup, nfactor) {
  repmat <- function(a, n, m) {
    kronecker(matrix(1, n, m), a)
  }
  M_X <- apply(X, 2, mean)
  S_X <- apply(X, 2, sd)
  if (NCOL(Y) == 1) {
    M_Y <- mean(Y)
    S_Y <- sd(Y)
  }
  else {
    M_Y <- apply(Y, 2, mean)
    S_Y <- apply(Y, 2, sd)
  }
  X <- scale(X)
  Y <- scale(Y)
  n <- nrow(Y)
  nn <- nrow(X)
  np <- ncol(X)
  nq <- ncol(Y)
  if (nn != n) {
    print("Incompatible # of rows for X and Y")
  }
  Yhatsup <- matrix(0, nfactor, nq)
  epsilon <- .Machine$double.eps
  U <- matrix(0, n, nfactor)
  C <- matrix(0, nq, nfactor)
  T <- matrix(0, n, nfactor)
  P <- matrix(0, np, nfactor)
  W <- matrix(0, np, nfactor)
  b <- matrix(0, 1, nfactor)
  Xres <- X
  Yres <- Y
  for (j in 1:nfactor) {
    plsr <- svd(t(Xres) %*% Yres)
    w <- plsr$u[, 1]
    c <- plsr$v[, 1]
    t <- normaliz(Xres %*% w)
    u <- Yres %*% c
    p <- t(Xres) %*% t
    b_l <- t(u) %*% t
    b[j] <- b_l
    P[, j] <- p
    W[, j] <- w
    T[, j] <- t
    U[, j] <- u
    C[, j] <- c
    Xres <- Xres - t %*% t(p)
    Yres <- Yres - repmat(b[j], n, nq) * (t %*% t(c))
    Wstar <- W %*% (MASS::ginv(t(P) %*% W))
    B_pls <- (Wstar * repmat(b, np, 1)) %*% t(C)
    Bpls_star <- t(repmat(t(S_X^-1), nq, 1)) * B_pls * repmat(t(S_Y),
                                                              np, 1)
    Bpls_star <- rbind(-M_X %*% Bpls_star, Bpls_star)
    Bpls_star[1, ] <- Bpls_star[1, ] + M_Y
    Yhatsup[j, ] <- c(1, xsup) %*% Bpls_star
  }
  return(Yhatsup)
}
