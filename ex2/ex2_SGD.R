## Generate some silly data
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
y = as.matrix(y,ncol = 1)

## Stochastic Gradient Descent algorithm
# calculate log likelihood function
L <- function(X,y,beta) {
  m = matrix(rep(1,n),ncol = 1)
  L = - t(y) %*% log(1/(1+exp(-X %*% beta)) + 1e-4) - t(m- y) %*% log(1 - 1/(1+exp(-X %*% beta)) + 1e-4)
  return (L)
}

# calculate Gradient 
grad <- function(X,y,beta){
  grad = - t(X) %*% (y - 1/(1+exp(-X %*% beta)))
  return(grad)
}

# Stochastic Gradient descent algorithm
SGD = function(x,y, alpha = 0.01, num.iterations = 5000, threshold = 1e-5){
  
  
  # initialize the parameters
  beta = matrix(rep(0,p),ncol = 1)
  
  
  # update parameters iteratively
  beta.path1_sgd = matrix(,nrow = num.iterations,ncol = p)
  logLike = rep(0, num.iterations)
  
  for (i in 1:num.iterations){
    # randomly sample one data point  
    sampleindex = sample(1:n,1)
    xvector = matrix(X[sampleindex,],nrow = 1)
    yscalar = y[sampleindex]
    beta = beta - alpha*grad(xvector,yscalar,beta) # update
    logLike[i] = L(X,y,beta)
    
    if (all(is.na(beta))){
      break
    } else {
      beta.path1_sgd[i,] = t(beta)
    }
    
    
    if (i >= 2){
      if (all(abs(t(beta) - beta.path1_sgd[i-1,]) < threshold)){
        break
      }
    }
    
  }
  return(list(coef = beta, iter=i, loglik=logLike[1:i]))
}


## Stochastic Gradient Descent algorithm using Robbins-Monro rule

# Stochastic Gradient descent algorithm
SGDRM = function(x,y, C = 0.5 , RMparm = -0.8, t0=1, num.iterations = 100000, threshold = 1e-5){
  
  # initialize the parameters
  beta = matrix(rep(0,p),ncol = 1)

  # update parameters iteratively
  beta.path1_sgdrm = matrix(,nrow = num.iterations,ncol = p)
  logLike = rep(0, num.iterations)
  
  for (i in 1:num.iterations){
    # randomly sample one data point  
    sampleindex = sample(c(1:n),1)
    xvector = matrix(X[sampleindex,], nrow = 1)
    yscalar = y[sampleindex]
    # use Robbins-Monro rule to update step size
    alpha = C * (i + t0)^RMparm
    # update beta
    beta = beta - alpha*grad(xvector,yscalar,beta)
    logLike[i] = L(X,y,beta)
    
    if (all(is.na(beta))){
      break
    } else {
      beta.path1_sgdrm[i,] = t(beta)
    }
    
    
    if (i >= 2){
      if (all(abs(t(beta) - beta.path1_sgdrm[i-1,]) < threshold)){
        break
      }
    }
  }
  return(list(coef = beta, iter=i, loglik=logLike[1:i]))
}


## Stochastic Gradient Descent algorithm using Polyak-Ruppert averaging
# Stochastic Gradient descent algorithm
SGDPR = function(x,y,alpha = 0.01, num.iterations = 20000, threshold = 1e-4){
  
  # initialize the parameters
  beta = matrix(rep(0,p),ncol = 1)
  
  # update parameters iteratively
  beta.path1_sgdpr = matrix(,nrow = num.iterations,ncol = p)
  logLike = rep(0, num.iterations)
  beta.PRavg = t(beta)
  beta.sum = t(beta)
  
  for (i in 1:num.iterations){
    # randomly sample one data point
    beta.old = beta
    sampleindex = sample(c(1:n),1)
    xvector = matrix(X[sampleindex,], nrow = 1)
    yscalar = y[sampleindex]
    beta = beta.old - alpha*grad1(xvector,yscalar,beta = beta.old) # update
    logLike[i] = L(X,y,beta)
    
    if (all(is.na(beta))){
      break
    } else {
      # report Polyak-Ruppert average value
      beta.sum = beta.sum + t(beta.old)
      beta.avg = beta.sum / i
      beta.path1_sgdpr[i,] = t(beta)
      beta.PRavg = rbind(beta.PRavg,beta.avg)
    }

    
    
    if (i >= 2){
      if (all(abs(t(beta) - beta.path1_sgdpr[i-1,]) < threshold)){
        break
      }
    }
  }
  return(list(coef = beta, iter=i, loglik=logLike[1:i]))
}

