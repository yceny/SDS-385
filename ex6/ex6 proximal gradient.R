########################################
# glmnet method
########################################
x = scale(dx, center = T)
y = scale(dy, center = T)
fit_glmnet = glmnet(x,y, family = 'gaussian')
coef_glmnet = fit_glmnet$beta[,which.min(fit_glmnet$lambda)]
fit_glmnet.cv = cv.glmnet(x,y)
coef_glmnet.cv = fit_glmnet.cv$glmnet.fit$beta[,which(fit_glmnet.cv$lambda == fit_glmnet.cv$lambda.min)]


########################################
# proximal gradient algorithm for lasso
#######################################
obj <- function(y,X,beta,lambda){
        A = y - X %*% beta
        l = (0.5 / nrow(X)) * crossprod(A) 
        phi = lambda * sum(abs(beta))
        obj = l + phi
        return(as.numeric(obj))
}

gradl <- function (y, X, beta){
    grad = (1 / nrow(X)) * (crossprod(X) %*% beta - t(X) %*% y) 
    return(grad)
}

prox <- function (x, gamma, lambda){
    r = sign(x) * pmax(rep(0, length(x)), abs(x) - gamma * lambda)
    return(r)
}

PG <- function(y,X,gamma = 0.01, lambda = 0.1, num.iteration = 1000, tol = 1e-5){
    p = ncol(X)
    
    # initialize beta
    beta_old = rep(0, p)
    
    # initialize matrix to hold all betas
    beta.path = matrix(NA,nrow = num.iteration,ncol = p)
    
    # create empty vector to stor objective value
    objective = c()
    
    for (t in 1:num.iteration){
        u = beta_old - gamma * gradl(y,X,beta_old) # update u
        beta_new = prox(u,gamma, lambda) # update beta
        beta.path[t,] = beta_new # store beta
        objective = c(objective, obj(y,X,beta_new, lambda))  # compute objective value
        
        if (abs(obj(y,X,beta_new, lambda) - obj(y,X,beta_old,lambda)) < tol){
            break
        }
        
        beta_old = beta_new
    }
    return(list(beta = beta_new, objective = objective, iter = t, betas = beta.path))
}

fit_pg = PG(y,x)
coef_pg = fit_pg$beta
    
###########################################
# Accerlerated proximal gradient algorithm
###########################################
APG <- function(y,X,gamma = 0.01, lambda = 0.1, num.iteration = 1000, tol = 1e-5){
    p = ncol(X)
    
    # initialize values
    beta_old = rep(0, p)
    z_old = rep(0,p)
    
    # initialize matrix to hold all betas
    beta.path = matrix(NA, nrow = num.iteration,ncol = p)
    beta.path[1,] = beta_old
    z.path = matrix(NA, nrow = num.iteration,ncol = p)
    z.path[1,] = z_old
    s = array(NA, dim = num.iteration)
    s[1] = 1
    
    # create empty vector to stor objective value
    objective = c()
    
    for (j in 2:num.iteration){
      s[j] = 0.5 * (1 + sqrt(1 + 4 * s[j - 1] ^2))
    }
    
    for (iter in 2: num.iteration){
        gradient <- gradl(y, X, z.path[iter-1, ])
        u <- z.path[iter-1, ] - gamma * gradient
        
        # Update beta
        beta.path[iter, ] <- prox(u, gamma, lambda)
        s[iter] <- (1 + sqrt(1 + 4 * s[iter-1]^2))/2
        z.path[iter, ] <- beta.path[iter, ] + ((s[iter-1] - 1)/s[iter]) * (beta.path[iter, ] - beta.path[iter-1, ])
        
        # Compute log-likelihood
        objective = c(objective, obj(y,X, beta.path[iter, ], lambda))
        
        if (abs(obj(y,X,beta.path[iter, ], lambda) - obj(y,X,beta.path[iter-1, ],lambda)) < tol){
           break
        }
        

    }
    return(list(beta = beta.path[iter,], objective = objective, iter = iter, betas = beta.path))
}

fit_apg = APG(y,x)
coef_apg = fit_apg$beta

