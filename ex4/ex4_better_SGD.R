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

## Stochastic Gradient Descent algorithm using mini batch
# input: sample.period: how often do a mini batch optimization
#        batch.size: when conducting a mini batch optimization, how much data do we sample, in terms of the proportion of the whole data
        
SGD_mini = function(x,y, sample.period = 50, batch.size = 0.07, num.iterations = 50000, threshold = 1e-5){
  
  N = nrow(x)
  
  # initialize the parameters
  beta = matrix(rep(0,p),ncol = 1)
  alpha = 0.01
  
  # update parameters iteratively
  beta.path1_sgdmini = matrix(,nrow = num.iterations,ncol = p)
  logLike = rep(0, num.iterations)
  avglogLike = rep(0, num.iterations)
  contributions = rep(0, num.iterations)
  
  for (i in 1:num.iterations){
    
    # optimize alpha over mini batch once every sample.period
    if (i %% sample.period == 0){
      size = batch.size * n
      gradient = matrix(,nrow = 6, ncol = size)
      for (j in 1:size){
        index = sample(1:n,1)
        xsample = matrix(x[index,],nrow = 1)
        ysample = y[index]
        gradient[,j] = grad(xsample,ysample,beta)
      }
      grad_avg = rowMeans(gradient)
      logLike_mini = rep(0, 1000)
      for (k in 1:1000){
        alpha = 0.0001 + (k-1) * 0.0001
        beta_mini = beta - alpha * grad_avg 
        logLike_mini[k] = L(X,y,beta_mini)
      }
      
      k = which.min(logLike_mini)
      alpha = 0.0001 + (k-1) * 0.0001
    }
    
    # basice SGD randomly sampling one data point  
    sampleindex = sample(1:n,1, replace = T)
    xvector = matrix(x[sampleindex,],nrow = 1)
    yscalar = y[sampleindex]
    beta = beta - alpha*grad(xvector,yscalar,beta)
    logLike[i] = L(X,y,beta)
    
    # Log-likelihood samples
    contribution = N * L(xvector, yscalar, beta)
    contributions[i] = contribution
    avglogLike[i] = contribution
    # Simple moving average (scale previous entry to the sum of all previous
    # entries, add new entry, divide by new size)
    if (i > 1) {
      avglogLike[i] = ((i - 1) * avglogLike[i - 1] + avglogLike[i]) / i
    }
    
    
    if (all(is.na(beta))){
      break
    } else {
      beta.path1_sgdmini[i,] = t(beta)
    }
    
    
    # convergence test using both acutal log likelihood and moving average log likelihood from samples
    if (i >= 2){
      if (abs(logLike[i] - logLike[i-1]) < threshold & abs(avglogLike[i] - avglogLike[i-1]) < threshold){
        break
      }
    }
    
  }
  return(list(coef = beta, iter=i, loglik=logLike, avglik = avglogLike))
}


## Stochastic Gradient Descent algorithm using AdaGrad
AdaGrad = function(x,y, alpha = 1, fudge_factor = 1e-6, num.iterations = 30000, threshold = 1e-5){
  
  N = nrow(x)
  
  # initialize the parameters
  beta = matrix(rep(0,p),ncol = 1)
  historical_grad = rep(0,p)
  
  
  # update parameters iteratively
  beta.path= matrix(,nrow = num.iterations,ncol = p)
  logLike = rep(0, num.iterations)
  avglogLike = rep(0, num.iterations)
  contributions = rep(0, num.iterations)
  
  for (i in 1:num.iterations){
    # randomly sample one data point  
    sampleindex = sample(1:n,1)
    xvector = matrix(x[sampleindex,],nrow = 1)
    yscalar = y[sampleindex]
    g = grad(xvector,yscalar,beta)
    historical_grad = historical_grad + g^2
    adjusted_alpha = alpha / (sqrt(historical_grad) + fudge_factor)
    beta = beta - adjusted_alpha * g # update
    logLike[i] = L(X,y,beta)
    
    # Log-likelihood samples
    contribution = N * L(xvector, yscalar, beta)
    contributions[i] = contribution
    avglogLike[i] = contribution
    # Simple moving average (scale previous entry to the sum of all previous
    # entries, add new entry, divide by new size)
    if (i > 1) {
      avglogLike[i] = ((i - 1) * avglogLike[i - 1] + avglogLike[i]) / i
    }
    
    if (all(is.na(beta))){
      break
    } else {
      beta.path[i,] = t(beta)
    }
    
    
    # convergence test using both acutal log likelihood and moving average log likelihood from samples
    if (i >= 2){
        if (abs(logLike[i] - logLike[i-1]) < threshold & abs(avglogLike[i] - avglogLike[i-1]) < threshold){
            break
        }
    }
  }
  return(list(coef = beta, iter=i, loglik=logLike, avglik = avglogLike))
}
