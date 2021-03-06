# generate data
set.seed(128)
n = 1000
x1 = rnorm(n)
x2 = rnorm(n)
x3 = rnorm(n)
x4 = rnorm(n)
x5 = rnorm(n)
z = 1 + 1 * x1 + 2 * x2 + 3 * x3 + 4 * x4 + 5 * x5
pr = 1/(1+exp(-z))
y = rbinom(n,1,pr)
X = cbind(1,x1,x2,x3,x4,x5)
p = 6
y = as.matrix(y, ncol = 1)

# using glm
fit_glm = glm(y ~ x1 + x2 + x3 + x4 +x5, family = binomial(link = 'logit'))
summary(fit_glm)

# calculate log likelihood function
L <- function(X,y,beta) {
  m = matrix(rep(1,n),ncol = 1)
  L = - t(y) %*% log(1/(1+exp(-X %*% beta))) - t(m- y) %*% log(1 - 1/(1+exp(-X %*% beta))) ## log might incur log0 problem, add a small variance to log(w + e) or log(1 - w + e)
  return (L)
}

# calculate Gradient 
grad <- function(X,y,beta){
  grad = - t(X) %*% (y - 1/(1+exp(-X %*% beta)))
  return(grad)
}

# calculate Hessian
hess <- function(X,y,beta){
  m = matrix(rep(1,n),ncol = 1)
  W = (1/(1+exp(-X %*% beta))) %*% t(m - 1/(1+exp(-X %*% beta)))
  diag_W = Diagonal(n,diag(W))
  hess = t(X) %*% diag_W %*% X
  return(hess)
}

##################################################################
# Gradient Descent 
##################################################################

GD = function(X,y, alpha = 0.01, num.iterations = 5000, threshold = 1e-5){
  
  
  # initialize the parameters
  beta = matrix(rep(0,p),ncol = 1)
  
  # update parameters iteratively
  beta.path1 = matrix(,nrow = num.iterations,ncol = p)
  logLike = rep(0, num.iterations)
  
  for (i in 1:num.iterations){
    beta = beta - alpha*grad(X,y,beta) # update
    logLike[i] = L(X,y,beta)
    
    if (all(is.na(beta))){
      break
    } else {
      beta.path1[i,] = t(beta)
    }
    
    #if (i >= 2){
    #  if (all(abs(t(beta) - beta.path1[i-1,]) < threshold)){
    #    break
    #  }
    #}
    
    if (i >= 2){
      if (abs(logLike[i] - logLike[i-1]) < threshold){
        break
      }
    }
  }
  
  return(list(coef = beta, iter=i, loglik=logLike))
}

###################################################################
## Netwon's method
###################################################################
Newton = function(X, y, num.iterations = 500, threshold = 1e-5){
  
  
  # initialize the parameters
  beta = matrix(rep(0,p),ncol = 1)
 
  
  # update parameters iteratively
  beta.path2 = matrix(,nrow = num.iterations,ncol = p)
  logLike = rep(0, num.iterations)

  for (i in 1:num.iterations){
    A = hess(X,y,beta)
    invA = solve(A) 
    beta = beta - invA %*% grad(X,y,beta) 
    logLike[i] = L(X,y,beta)
    
    if (all(is.na(beta))){
      break
    } else {
      beta.path2 = rbind(beta.path2,beta)
    }
    
    #if (i >= 2){
    #  if (all(abs(t(beta) - beta.path1[i-1,]) < threshold)){
    #    break
    #  }
    #}
    
    if (i >= 2){
      if (abs(logLike[i] - logLike[i-1]) < threshold){
        break
      }
    }
  }
  return(list(coef = beta, iter=i, loglik=logLike))
}





