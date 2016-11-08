library(Matrix)
library(ggplot2)
library(microbenchmark)

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


#####################################################
# objective function: 0.5||y-x||^2 + gamma * ||Dx||
#####################################################

obj <- function(y, beta, gamma, D){
    l = 0.5 * crossprod(y - beta)
    phi = gamma * sum(abs(D %*% beta))
    r = l + phi
    return(as.numeric(r))
}

# proximal operator, which is basically soft thresholding method
prox <- function (u, penalty){
    a = abs(u) - penalty
    b = cbind(rep(0,length(a)), a)
    r = sign(u) * apply(b,1,max)
    return(r)
}

# ADMM algorithm with constant rho 
admm <- function(y, D, rho = 0.01, gamma = 0.01, num.iteration = 1000, ABSTOL = 1e-4, RELTOL = 1e-2){
    p = dim(D)[2]
    n = dim(D)[1]
    
    # initialize values
    beta = rep(0,p) # coefficients of interest
    z = rep(0,n)
    # lambada = rep(0,n) # lagrangian multiplier
    u = rep(0, n) # scaled dual variables
    
    # initialize matrix to hold data
    beta.path = matrix(0, nrow = num.iteration,ncol = p)
    # z.path = matrix(0, nrow = num.iteration, ncol = n)
    u.path = matrix(0, nrow = num.iteration, ncol = n)
    
    # create empty vector to stor objective value
    objective = c()
    
    # when rho is constant, we cache matrix inverse used in the first step
    inv_cache = solve(diag(p) + rho * crossprod(D))
    
    for (i in 1:num.iteration){
        # step 1
        # beta = inv_cache %*% (y + rho * (t(D) %*% z) - t(D) %*% lambda)
        beta = inv_cache %*% (y - rho * t(D) %*% (z - u) )
        
        
        # step 2
        # z = prox(D %*% beta + lambda / rho, gamma / rho)
        z = prox(D %*% beta + u, gamma / rho)
        
        # step 3
        # lambda = lambda + rho * (D %*% beta - z)
        u = u + D %*% beta - z
        
        beta.path[i,] = as.vector(beta)
        # z.path[i,] = as.vector(z)
        u.path[i,] = as.vector(u)
        objective = c(objective, obj(y, beta, gamma, D))
        
        if (i >= 2){
            # dual.res = - rho * (z.path[i,] - z.path[(i-1),])
            dual.res = - rho * t(D) %*% (u.path[i,] - u.path[(i - 1), ])
            prime.res = D %*% beta - z
            dual.tol = sqrt(n) * ABSTOL + RELTOL * sqrt(sum((rho * t(D) %*% u) ^ 2))
            prime.tol = sqrt(p) * ABSTOL + RELTOL * max(sqrt(sum((D %*% beta) ^ 2)), sqrt(sum(z ^ 2)))
            if (sqrt(sum(dual.res ^ 2)) < dual.tol & sqrt(sum(prime.res ^ 2)) < prime.tol){
                break
            }
        }
        
    }
    
    return(list(beta = beta, iter = i, objective = objective))
    
}


# more efficient ADMM
fast_admm <- function(y, D, rho = 0.01, gamma = 0.1, num.iteration = 1000, ABSTOL = 1e-4, RELTOL = 1e-2){
    p = dim(D)[2]
    n = dim(D)[1]
    
    # initialize values
    beta = rep(0,p) # coefficients of interest
    z = rep(0,n)
    r = rep(0,p) # slack variable 1
    s = rep(0,n) # slack variable 2
    t = rep(0,n) # where to scale
    u = rep(0,p) # where to scale
    
    # initialize matrix to hold data
    beta.path = matrix(0, nrow = num.iteration,ncol = p)
    # z.path = matrix(0, nrow = num.iteration, ncol = n)
    t.path = matrix(0, nrow = num.iteration, ncol = n)
    
    # create empty vector to stor objective value
    objective = c()
    
    # when rho is constant, we cache matrix inverse used in the first step
    inv_cache = solve((1 + rho) * diag(p))
    inv_cache1 = solve(diag(p) + crossprod(D))
    
    for (i in 1:num.iteration){
        # step 1
        beta = inv_cache %*% (y + rho * r)
        
        # step 2
        z = prox(s - t, gamma / rho) # shall we divide t by rho
        
        # step 3
        w = beta + u
        v = z + t
        r = inv_cache1 %*% (w + t(D) %*% v)
        s = D %*% r
        
        # step 4
        u = u + beta - r
        t = t + z -s
        
        beta.path[i,] = as.vector(beta)
        # z.path[i,] = as.vector(z)
        t.path[i,] = as.vector(t)
        objective = c(objective, obj(y, beta, gamma, D))
        
        if (i >= 2){
            # dual.res = - rho * (z.path[i,] - z.path[(i-1),])
            dual.res = - rho * t(D) %*% (t.path[i,] - t.path[(i - 1), ])
            prime.res = D %*% beta - z
            dual.tol = sqrt(n) * ABSTOL + RELTOL * sqrt(sum((rho * t(D) %*% t) ^ 2))
            prime.tol = sqrt(p) * ABSTOL + RELTOL * max(sqrt(sum((D %*% beta) ^ 2)), sqrt(sum(z ^ 2)))
            if (sqrt(sum(dual.res ^ 2)) < dual.tol & sqrt(sum(prime.res ^ 2)) < prime.tol){
                break
            }
        }
        
    }
    
    return(list(beta = beta, betas = beta.path, iter = i, objective = objective))
    
}


fmri = as.matrix(fmri)
y = Matrix::Matrix(as.vector(fmri))
D = makeD2_sparse(nrow(fmri), ncol(fmri))
fit_admm = admm(y, D)
fit_fast.admm = fast_admm(y,D)

admmMatrix = matrix(fit_admm$beta, nrow = nrow(fmri))
fast_admmMatrix = matrix(fit_fast.admm$beta, nrow = nrow(fmri))

# plot results
par(mfrow = c(1, 3))
image(fmri, main = "Raw")
image(admmMatrix, main = "ordinary ADMM")
image(fast_admmMatrix, main = "efficient ADMM")

# plot convergence
par(mfrow = c(2,1))
plot(fit_admm$objective)
plot(fit_fast.admm$objective)

microbenchmark::microbenchmark(admm(y, D), fast_admm(y,D), times = 10L)
