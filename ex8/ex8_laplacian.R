library(Matrix)
library(ggplot2)

# Construct the first-difference matrix, D, over a grid graph, stolen from James
makeD2_sparse = function (dim1, dim2)  {
    require(Matrix)
    D1 = bandSparse(dim1 * dim2, m = dim1 * dim2, k = c(0, 1), 
                    diagonals = list(rep(-1, dim1 * dim2), rep(1, dim1 * 
                                                                   dim2 - 1)))
    D1 = D1[(seq(1, dim1 * dim2)%%dim1) != 0, ]
    D2 = bandSparse(dim1 * dim2 - dim1, m = dim1 * dim2, k = c(0, 
                                                               dim1), diagonals = list(rep(-1, dim1 * dim2), rep(1, 
                                                                                                                 dim1 * dim2 - 1)))
    return(rBind(D1, D2))
}


# a direct solver using sparse LU factorization of C
lu <- function(C, b){
    # input C should be a sparse matrix
    # sparseLU (see ?lu):
    # decomposition: A = P' L U Q
    # problem: A x = b
    # --------
    # P' L U Q x = b
    # P' L y = b
    # P P' L y = L y= P b ***
    # U Q x = y
    # U z = y ***
    # Q' Q x = x = Q' z ***
    # --------
    # So the results we need are
    # L y = P b
    # U z = y
    # x = Q' z
    lu.C = expand(Matrix::lu(C))
    y = Matrix::solve(lu.C$L, lu.C$P %*% b, system = 'L')
    z = Matrix::solve(lu.C$U, y, system = 'L')
    x = t(lu.C$Q) %*% z
    return(x)
}

# Gauss Seidel method
GS <- function(C,b, num.iteration = 1000, tol = 1e-8){
    # input C should be a sparse matrix
    p = dim(C)[2]
    L_star = tril(C) # lower triangle of sparse matrix C, including diagonal values
    U = triu(C, 1) # strict upper triangle of sparse matrix C, excluding diagonal values
    x_old = rep(0, p)
    for (i in 1:num.iteration){
        x = Matrix::solve(L_star, b - U %*% x_old, system = 'L')
        if (sum((x - x_old)^2) < tol){
            break
        } else {
            x_old = x
        }
    }
    return(x)
}

# Jacobi method
JCB <- function(C,b, num.iteration = 1000, tol = 1e-8){
    # input C should be a sparse matrix
    p = dim(C)[2]
    invD = 1 / diag(C)
    R = C - Diagonal(x = diag(C))
    
    x_old = rep(0, p)
    for (i in 1:num.iteration){
        x = invD * (b - R %*% x_old)
        if (sum((x - x_old)^2) < tol){
            break
        } else {
            x_old = x
        }
    }
    return(x)
}

# compare different methods
fmri = as.matrix(fmri)
y = Matrix::Matrix(as.vector(fmri))
D = makeD2_sparse(nrow(fmri), ncol(fmri))
L = t(D) %*% D

# change lambda to see what happens
nn = dim(L)[1]
lambda = 2
C = Diagonal(nn) + lambda * L
x = lu(C, c(fmri))
gs.x = GS(C, y)
jcb.x = JCB(C, y)

xMat = matrix(x, nrow = nrow(fmri))
gsxMat = matrix(gs.x, nrow = nrow(fmri))
jcbxMat = matrix(jcb.x, nrow = nrow(fmri))

# Plot
par(mfrow = c(2, 2))
image(fmri, main = "Raw")
image(xMat, main = "LU")
image(gsxMat, main = "Gauss Seidel")
image(jcbxMat, main = "Jacobi")

