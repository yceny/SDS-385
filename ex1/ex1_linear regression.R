library(Matrix)



nn = 10
pp = 5
X = matrix(runif(nn*pp),nrow = nn, ncol = pp)
W = diag(runif(nn),nrow = nn, ncol = nn)
y = rnorm(nn)


## solve linear equations using matrix inverse
matrixInv <- function(X,W,y){
  A = t(X)%*%W%*%X
  b = t(X)%*%W%*%y
  beta_hat = solve(A,b)
  return (beta_hat)
}


## solve linear equations using Cholesky decomposition

choleskyDecomp <- function(X,W,y){
  A = t(X)%*%W%*%X
  b = t(X)%*%W%*%y
  upperA = chol(A)
  Y = solve(t(upperA),b)
  beta_hat = solve(upperA,Y)
  return(beta_hat)
}


## solve linear equations using LU decomposition
LUDecomp <- function(X,W,y){
  A = t(X)%*%W%*%X
  b = t(X)%*%W%*%y
  LUDA = lu(A)
  eLUDA = expand(LUDA)
  lowerA = eLUDA$L
  upperA = eLUDA$U
  Y = solve(lowerA,b)
  beta_hat = solve(upperA,Y)
  return(beta_hat)
}


## solve linear equations using QR decomposition
QRDecomp <- function(X,W,y){
  A = sqrt(W) %*% X
  Aqr = qr(A)
  Q = qr.Q(Aqr)
  R = qr.R(Aqr)
  b = t(Q) %*% sqrt(W) %*% y
  beta_hat = solve(R,b)
  return(beta_hat)
}


## solve linear equations using SVD
SVDecomp <- function(X,W,y){
  A = t(X)%*%W%*%X
  b = t(X)%*%W%*%y
  Asvd = svd(A)
  U = Asvd$u
  V = Asvd$v
  S = diag(Asvd$d, nrow = pp)
  b_new = t(U) %*% b
  Y = solve(S,b_new)
  beta_hat = solve(t(V),Y)
  return(beta_hat)
}

## ------------------------------------------------------------------------------------------------
## Start Simulation
speed1 = matrix(,nrow = 15, ncol = 20)
speed2 = matrix(,nrow = 15, ncol = 20)
speed3 = matrix(,nrow = 15, ncol = 20)
speed4 = matrix(,nrow = 15, ncol = 20)
speed5 = matrix(,nrow = 15, ncol = 20)

for (pp in seq(100,800,50)){
  for (nn in seq(100,2000,100)){
    # generate data
    set.seed(258)
    X = matrix(rnorm(nn*pp), nrow = nn)
    W = diag(runif(nn),nrow = nn, ncol = nn)
    y = rnorm(nn)
    
    x = try(matrixInv(X,W,y),silent = TRUE)
    if ('try-error' %in% class(x)){
      i = (pp-100)/50 + 1
      j = nn / 100
      speed1[i,j] = 'NA'
    } else {
      tMatrixInverse = microbenchmark(matrixInv(X,W,y))
      tMatrixInverse = summary(tMatrixInverse)$mean
      i = (pp-100)/50 + 1
      j = nn / 100
      speed1[i,j] = tMatrixInverse
    }
  
    x = try(choleskyDecomp(X,W,y),silent = TRUE)
    if ('try-error' %in% class(x)){
      i = (pp-100)/50 + 1
      j = nn / 100
      speed2[i,j] = 'NA'
    } else {
      tCholeskyDecomp = microbenchmark(choleskyDecomp(X,W,y))
      tCholeskyDecomp = summary(tCholeskyDecomp)$mean
      i = (pp-100)/50 + 1
      j = nn / 100
      speed2[i,j] = tCholeskyDecomp
    }
    
    x = try(LUDecomp(X,W,y),silent = TRUE)
    if ('try-error' %in% class(x)){
      i = (pp-100)/50 + 1
      j = nn / 100
      speed3[i,j] = 'NA'
    } else {
      tLUDecomp = microbenchmark(LUDecomp(X,W,y))
      tLUDecomp = summary(tLUDecomp)$mean
      i = (pp-100)/50 + 1
      j = nn / 100
      speed3[i,j] = tLUDecomp
    }
    
    x = try(QRDecomp(X,W,y),silent = TRUE)
    if ('try-error' %in% class(x)){
      i = (pp-100)/50 + 1
      j = nn / 100
      speed4[i,j] = 'NA'
    } else {
      tQRDecomp = microbenchmark(QRDecomp(X,W,y))
      tQRDecomp = summary(tQRDecomp)$mean
      i = (pp-100)/50 + 1
      j = nn / 100
      speed4[i,j] = tQRDecomp
    }
 
    x = try(SVDecomp(X,W,y),silent = TRUE)
    if ('try-error' %in% class(x)){
      i = (pp-100)/50 + 1
      j = nn / 100
      speed5[i,j] = 'NA'
    } else {
      tSVDecomp = microbenchmark(SVDecomp(X,W,y))
      tSVDecomp = summary(tSVDecomp)$mean
      i = (pp-100)/50 + 1
      j = nn / 100
      speed5[i,j] = tSVDecomp
    }
  }
}


## ------------------------------------------------------------------------------------------------------------
## Start Simulation for sparse matrix
speed1 = matrix(,nrow = 15, ncol = 20)
speed2 = matrix(,nrow = 15, ncol = 20)
speed3 = matrix(,nrow = 15, ncol = 20)
speed4 = matrix(,nrow = 15, ncol = 20)
speed5 = matrix(,nrow = 15, ncol = 20)

for (pp in seq(100,800,50)){
  for (nn in seq(100,2000,100)){
    # generate data
    set.seed(258)
    X = matrix(rnorm(nn*pp), nrow = nn)
    mask = matrix(rbinom(nn*pp,1,0.05),nrow = nn)
    X = mask*X
    W = diag(runif(nn),nrow = nn, ncol = nn)
    y = rnorm(nn)
    
    x = try(matrixInv(X,W,y),silent = TRUE)
    if ('try-error' %in% class(x)){
      i = (pp-100)/50 + 1
      j = nn / 100
      speed1[i,j] = 'NA'
    } else {
      tMatrixInverse = microbenchmark(matrixInv(X,W,y))
      tMatrixInverse = summary(tMatrixInverse)$mean
      i = (pp-100)/50 + 1
      j = nn / 100
      speed1[i,j] = tMatrixInverse
    }
    
    x = try(choleskyDecomp(X,W,y),silent = TRUE)
    if ('try-error' %in% class(x)){
      i = (pp-100)/50 + 1
      j = nn / 100
      speed2[i,j] = 'NA'
    } else {
      tCholeskyDecomp = microbenchmark(choleskyDecomp(X,W,y))
      tCholeskyDecomp = summary(tCholeskyDecomp)$mean
      i = (pp-100)/50 + 1
      j = nn / 100
      speed2[i,j] = tCholeskyDecomp
    }
    
    x = try(LUDecomp(X,W,y),silent = TRUE)
    if ('try-error' %in% class(x)){
      i = (pp-100)/50 + 1
      j = nn / 100
      speed3[i,j] = 'NA'
    } else {
      tLUDecomp = microbenchmark(LUDecomp(X,W,y))
      tLUDecomp = summary(tLUDecomp)$mean
      i = (pp-100)/50 + 1
      j = nn / 100
      speed3[i,j] = tLUDecomp
    }
    
    x = try(QRDecomp(X,W,y),silent = TRUE)
    if ('try-error' %in% class(x)){
      i = (pp-100)/50 + 1
      j = nn / 100
      speed4[i,j] = 'NA'
    } else {
      tQRDecomp = microbenchmark(QRDecomp(X,W,y))
      tQRDecomp = summary(tQRDecomp)$mean
      i = (pp-100)/50 + 1
      j = nn / 100
      speed4[i,j] = tQRDecomp
    }
    
    x = try(SVDecomp(X,W,y),silent = TRUE)
    if ('try-error' %in% class(x)){
      i = (pp-100)/50 + 1
      j = nn / 100
      speed5[i,j] = 'NA'
    } else {
      tSVDecomp = microbenchmark(SVDecomp(X,W,y))
      tSVDecomp = summary(tSVDecomp)$mean
      i = (pp-100)/50 + 1
      j = nn / 100
      speed5[i,j] = tSVDecomp
    }
  }
}






