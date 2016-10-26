########################################
# glmnet method
########################################
library(glmnet)

x = scale(dx, center = T) # after scaling and centering the data, we do not need to include constant 
y = scale(dy, center = T)
fit_glmnet = glmnet(x,y, family = 'gaussian')
coef_glmnet = fit_glmnet$beta[,which.min(fit_glmnet$lambda)]
fit_glmnet.cv = cv.glmnet(x,y)
coef_glmnet.cv = fit_glmnet.cv$glmnet.fit$beta[,which(fit_glmnet.cv$lambda == fit_glmnet.cv$lambda.min)]

###############################
# ADMM
###############################

# objective value for the 0.5/nrow(X) * ||y-X*beta||^2 + gamma * ||beta||
obj <- function(y,X,beta,gamma){
  A = y - X %*% beta
  # l = 0.5 * crossprod(A) 
  l = (0.5 / nrow(X)) * crossprod(A)
  phi = gamma * sum(abs(beta))
  obj = l + phi
  return(as.numeric(obj))
}

# proximal operator, which is basically soft thresholding method
prox <- function (u, penalty){
  a = abs(u) - penalty
  b = cbind(rep(0,length(a)), a)
  r = sign(u) * apply(b,1,max)
  return(r)
}

# ADMM algorithm with constant rho 
admm <- function(y, X, rho = 0.01, gamma = 0.1, num.iteration = 1000, ABSTOL = 1e-4, RELTOL = 1e-2){
  p = ncol(X)
  
  # initialize values
  beta = rep(0,p) # coefficients of interest
  z = rep(0,p)
  lambda = rep(0,p) # lagrangian multiplier
  
  # initialize matrix to hold data
  beta.path = matrix(0, nrow = num.iteration,ncol = p)
  z.path = matrix(0, nrow = num.iteration, ncol = p)
  
  # create empty vector to stor objective value
  objective = c()
  
  # when rho is constant, we cache matrix inverse used in the first step
  inv_cache = solve((1 / nrow(X)) * crossprod(X) + rho * diag(p))
  
  for (i in 1:num.iteration){
    # step 1
    beta = inv_cache %*% ((1 / nrow(X)) * t(X) %*% y + rho * z - lambda)
    
    # step 2
    z = prox(beta + lambda / rho, gamma / rho)
    
    # step 3
    lambda = lambda + rho * (beta - z)
    
    beta.path[i,] = beta
    z.path[i,] = z
    objective = c(objective, obj(y, X, beta, gamma))
    
    if (i >= 2){
      dual.res = - rho * (z.path[i,] - z.path[(i-1),])
      prime.res = beta - z
      dual.tol = sqrt(p) * ABSTOL + RELTOL * sqrt(sum(lambda ^ 2))
      prime.tol = sqrt(p) * ABSTOL + RELTOL * max(sqrt(sum(beta ^ 2)), -sqrt(sum(z ^ 2)))
      if (sqrt(sum(dual.res ^ 2)) < dual.tol & sqrt(sum(prime.res ^ 2)) < prime.tol){
        break
      }
    }
    
  }
  
  return(list(beta = beta, betas = beta.path, iter = i, objective = objective))
  
}

fit_admm.glm = admm(y,x, gamma = fit_glmnet.cv$lambda.min)
coef_admm = fit_admm.glm$beta

# ADMM algorithm with changing rho 
Ad_admm <- function(y, X, maxRho = 5, gamma = 0.1, num.iteration = 1000, ABSTOL = 1e-4, RELTOL = 1e-2, mu = 10, tau.incr = 2, tau.decr = 2){
  p = ncol(X)
  
  # initialize values
  rho = 0.0001
  beta = rep(0,p) # coefficients of interest
  z = rep(0,p)
  lambda = rep(0,p) # lagrangian multiplier
  
  # initialize matrix to hold all betas
  beta.path = matrix(,nrow = num.iteration,ncol = p)
  z.path = matrix(0, nrow = num.iteration, ncol = p)
  
  # create empty vector to stor objective value
  objective = c()
  
  
  for (i in 1:num.iteration){
    # step 1
    inv_cache = solve((1 / nrow(X)) * crossprod(X) + rho * diag(p))
    beta = inv_cache %*% ((1 / nrow(X)) * t(X) %*% y + rho * z - lambda)
    
    # step 2
    z = prox(beta + lambda / rho, gamma / rho)
    
    # step 3
    lambda = lambda + rho * (beta - z)
    
    beta.path[i,] = beta
    z.path[i,] = z
    objective = c(objective, obj(y, X, beta, gamma))
    
    
    if (i >= 2){
      dual.res = - rho * (z.path[i,] - z.path[(i-1),])
      prime.res = beta - z
      dual.tol = sqrt(p) * ABSTOL + RELTOL * sqrt(sum(lambda ^ 2))
      prime.tol = sqrt(p) * ABSTOL + RELTOL * max(sqrt(sum(beta ^ 2)), -sqrt(sum(z ^ 2)))
      if (sqrt(sum(dual.res ^ 2)) < dual.tol & sqrt(sum(prime.res ^ 2)) < prime.tol){
        break
      }
      
      # update rho
      if (sqrt(sum(prime.res ^ 2)) > mu * sqrt(sum(dual.res ^ 2))) {
        rho = tau.incr * rho
      } else if (sqrt(sum(dual.res ^ 2)) > sqrt(sum(prime.res ^ 2))) {
        rho = rho / tau.decr
      } else {
        rho = rho
      }
    }
    
  }
  
  return(list(beta = beta, betas = beta.path, iter = i, objective = objective))
  
}

fit_ad_admm.glm = Ad_admm(y,x, gamma = fit_glmnet.cv$lambda.min)
coef_ad_admm = fit_ad_admm.glm$beta


####################
# Plot
####################
plot(fit_admm.glm$objective, type = 'l')
lines(fit_ad_admm.glm$objective, col = 'red')
legend("top", c('ADMM with constant rho', 'ADMM with changing rho'), bty = 'n', fill = c('black','red'))






