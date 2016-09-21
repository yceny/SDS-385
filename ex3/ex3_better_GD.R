## generate some silly data
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
    L = - t(y) %*% log(1/(1+exp(-X %*% beta)) + 1e-4) - t(m- y) %*% log(1 - 1/(1+exp(-X %*% beta)) + 1e-4)
    return (L)
}

# calculate Gradient 
grad <- function(X,y,beta){
    grad = - t(X) %*% (y - 1/(1+exp(-X %*% beta)))
    return(grad)
}

# backtrack line search condition
backtrack = function (X,y,beta,oldLogL,c,alpha,gL,direction){
  new.beta = beta + alpha * direction
  newCond = L(X,y,new.beta)
  change.term = c * alpha * crossprod(gL,direction)
  oldCond = oldLogL + change.term
  return (newCond <= oldCond)
}

############################################################
## Gradient descent algorithm using backtrack line search
############################################################

BGD = function(X,y,alpha0 = 1, rho = 0.7, c =1e-3, num.iterations = 10000, threshold = 1e-5){
    
    # initialize the parameters
    beta = matrix(rep(0,p),ncol = 1)
    
    # update parameters iteratively
    beta.path1_backtracking = matrix(,nrow = num.iterations,ncol = p)
    logLike = rep(0, num.iterations)
    
    for (i in 1:num.iterations){
        # backtracking line search to determine alpha
        alpha = alpha0
        gL = grad(X,y,beta)
        direction = -gL
        oldLogL = L(X,y,beta)
        while (!backtrack(X,y,beta,oldLogL,c,alpha,gL,direction)){
          alpha = alpha * rho
        }
        

        # update new beta
        beta = beta - alpha*grad(X,y,beta) 
        logLike[i] = L(X,y,beta)
        
        if (all(is.na(beta))){
            break
        } else {
            beta.path1_backtracking[i,] = t(beta)
        }
        
        #if (i >= 2){
        #    if (all(abs(t(beta) - beta.path1_backtracking[i-1,]) < threshold)){
        #        break
        #    }
        #}
        
        if (i >= 2){
          if (abs(logLike[i] - logLike[i-1]) < threshold){
            break
          }
        }
    }
    return(list(coef = beta, iter=i, loglik=logLike))
}

############################################################################
## Quasi Newton using BFGS
############################################################################

qNewton = function(X,y,alpha = 0.01, num.iterations = 5000, threshold = 1e-5){
  # initialize the parameters
  beta = matrix(rep(0,p,ncol = 1))
  B = diag(p)
  
  # update parameters iteratively
  beta.path1_qn = matrix(,nrow = num.iterations,ncol = p)
  logLike = rep(0, num.iterations)
  
  beta.old = beta
  
  for (i in 1:num.iterations){
    beta.new = beta.old - alpha * solve(B, grad(X,y,beta.old))
    logLike[i] = L(X,y,beta.new)
    
    # update Hessian using BFGS
    if (all(is.na(beta.new))){
      break
    } else {
      beta.path1_qn[i,] = t(beta.new)
      yy = grad(X,y,beta.new) - grad(X,y, beta.old)
      s = beta.new - beta.old
      beta.old = beta.new
      f1 = as.numeric(t(s) %*% B %*% s)
      A1 = B %*% s %*% t(s) %*% B / f1
      f2 = as.numeric(t(yy) %*% s)
      A2 = crossprod(t(yy)) / f2
      B = B - A1 + A2
    }
    
    if (i >= 2){
      if (abs(logLike[i] - logLike[i-1]) < threshold){
        break
      }
    }
    
  }
  return(list(coef = beta.old, iter=i, loglik=logLike))
}


###################################################################################################
## qausi-Newton's method using BFGS update and backtrack line search
###################################################################################################
qNewton_bb = function(X,y,alpha0 = 1, rho = 0.7, c =1e-3, num.iterations = 10000, threshold = 1e-5){
    # initialize the parameters
    beta = matrix(rep(0,p,ncol = 1))
    B = diag(p)
    
    # update parameters iteratively
    beta.path1_BFGS = matrix(,nrow = num.iterations,ncol = p)
    logLike = rep(0, num.iterations)
    
    beta.old = beta
    for (i in 1:num.iterations){
        
        # backtracking line search to determine alpha
        alpha = alpha0
        gL = grad(X,y,beta.old)
        direction = -gL
        oldLogL = L(X,y,beta.old)
        while (!backtrack(X,y,beta.old,oldLogL,c,alpha,gL,direction)){
          alpha = alpha * rho
        }
        
        # update estimates based on step length (alpha) 
        beta.new = beta.old - alpha * solve(B, grad(X,y,beta.old))
        logLike[i] = L(X,y,beta.new)
        
        # update Hessian using BFGS
        if (all(is.na(beta.new))){
          break
        } else {
          beta.path1_BFGS[i,] = t(beta.new)
          yy = grad(X,y,beta.new) - grad(X,y, beta.old)
          s = beta.new - beta.old
          beta.old = beta.new
          f1 = as.numeric(t(s) %*% B %*% s)
          A1 = B %*% s %*% t(s) %*% B / f1
          f2 = as.numeric(t(yy) %*% s)
          A2 = crossprod(t(yy)) / f2
          B = B - A1 + A2
        }
        
        #if (i >= 2){
        #    if (all(abs(beta.new - beta.path1_BFGS[i-1,]) < threshold)){
        #        break
        #    }
        #}
        
        if (i >= 2){
          if (abs(logLike[i] - logLike[i-1]) < threshold){
            break
          }
        }
        
    }
    return(list(coef = beta.old, iter=i, loglik=logLike))
}

