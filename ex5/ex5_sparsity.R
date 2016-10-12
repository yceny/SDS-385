# objective function
f <- function(y, lambda, theta){
    f = 0.5 * (y - theta)^2 + lambda * abs(theta)
    return(f)
}

# soft thresholding
st <- function(y,lambda){
        v = sign(y) * max(0, abs(y) - lambda)
    return(v)
}

# hard thresholding
ht <- function(y, lambda){
        if (y >= lambda){
            return(y)
        } else{
            return(0)
        }
}

# generate sparse vecta theta
set.seed(2016)
toy <- function(n, sparsity.rate, lambda){
# function input:
    # n = sample size
    # sparsity.rate = percentage of theta that equals to 0
    # lambda = vector of lambda values to test
# function output:
    # theta = true values of generated theta
    # theta_est = estimated theta values
    theta = matrix(,nrow = n, ncol = 1)
    sparse.index = sample(c(1:n), n * sparse.rate)
    theta[sparse.index, ] = 0
    theta[-sparse.index, ] = runif((1 - sparse.rate) * n)
    sigma = theta # set variance equal to mean
    
    theta_hat = matrix(, nrow = n, ncol = length(lambda))
    MSE = rep(0, length(lambda))
    for (k in 1:length(lambda)){
        theta_est = c()
        for (i in 1:n){
            y = rnorm(1, mean = theta[i], sd = sqrt(sigma[i]))
            lambda.new = lambda[k] * sigma[i]
            theta_est = c(theta_est, st(y, lambda.new))
        }
        theta_hat[ , k] = theta_est 
        MSE[k] = sum((theta_est - theta)^2) / n
    }
    
    return(list(theta = theta, lambda = lambda, theta_hat = theta_hat, mse = MSE, sparsity = sparsity.rate))
}

# begin tests

sparsity = c(0.2, 0.4, 0.6, 0.8)
lambda = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
L = length(sparsity)
for (l in 1:L){
    fit = toy(n = 1000, sparsity.rate = sparsity[l], lambda = lambda)
    par(mfrow = c(3,4))
    for (j in 1 : length(fit$lambda)){
        plot(fit$theta, fit$theta_hat[,j], main = paste('sparsity = ', sparsity[l],'lambda =', fit$lambda[j]), xlab = 'true theta', ylab = 'estimated theta')
    }
    plot(fit$lambda, fit$mse, type = 'l', main = paste('Optimal lambda for sparsity = ', sparsity[l]))
}


# using glmnet to fit diabetes data
x = scale(dx)
y = scale(dy)
fit_glm = glmnet(x,y, family = 'gaussian') # setting parameter, alpha, to specify regularized term, default is alpha = 1, which is lasso; x must be a matrix and y must be a vector
plot(fit_glm, xvar = 'lambda')

# in sample mse
Lambda = fit_glm$lambda
betas = fit_glm$beta
N = length(Lambda)
mse.insample = rep(0, N)
for (i in 1:N){
    mse.insample[i] = sum((y - x %*% betas[,i])^2) / nrow(x)
}
plot(Lambda, mse.insample)
line(log(Lambda), mse.insample)

# cross validation to choose optimal lambda
# benchmark model using cv.glmnet in glm package
fit_cvglm = cv.glmnet(x,y)
fit_cvglm$lambda.min

# build cross validation function (stolen from Cross Validated "http://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation")
data = cbind(y,x)
test_lambda = fit_glm$lambda
cross_validation <- function(data, fold = 10, test_lambda){
    folds <- cut(seq(1,nrow(data)),breaks=fold,labels=FALSE)
    
    pred_error = matrix(,nrow = fold, ncol = length(test_lambda))
    
    #Perform 10 fold cross validation
    for(i in 1:fold){
        #Segement your data by fold using the which() function 
        testIndexes <- which(folds==i,arr.ind=TRUE)
        testData <- data[testIndexes, ]
        trainData <- data[-testIndexes, ]
        
        train_fit = glmnet(x = trainData[,-1], y = trainData[,1], family = 'gaussian',lambda = test_lambda)
        test_pred = predict(train_fit, newx = testData[,-1], s = test_lambda)
        prederr = apply(test_pred, 2, function(yhat) sum((yhat - testData[,1])^2) / nrow(testData))
        pred_error[i,] = prederr                
    }
    
    return(list(lambda = test_lambda, MOOSE = colMeans(pred_error)))
}


# C_p statistics to choose optimal lambda
Lambda = fit_glm$lambda
betas = fit_glm$beta
N = length(Lambda)
n = nrow(x)
Cp = c()
for (i in 1:N){
    df = fit_glm$df[i]
    mse = sum((y - x %*% betas[,i])^2) / n
    sigma2 = var(y - x %*% betas[,i])
    Cp = c(Cp, mse + 2 * df * sigma2 / n)
}

plot(Lambda,Cp, type = 'l')
Lambda[which.min(Cp)]

# plot 
plot(log(Lambda),mse.insample, type = 'l', xlab = 'log(lambda)', ylab = 'Error')
lines(log(Lambda), cv$MOOSE, col = 2)
lines(log(Lambda), Cp, col = 3)
legend('topleft', legend = c('ISMSE','MOOSE','Cp'), col = 1:3, lty = 1)
