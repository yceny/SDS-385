library(Matrix)
library(inline)
library(Rcpp)
library(RcppEigen)
sourceCpp('sgd_lasso_logit.cpp')

n = length(y_label)
p = ncol(X)

# column-oriented storage of each observation, matrix operation in Eigen is defaulted as column major
tX = t(X)
init_beta = rep(0.0, p)

# test on small dataset
fit_test = sgd_lasso_logit(tX[1:10,1:1000], y_label[1:1000], beta0 = init_beta[1:10], stepsize = 1, discount = 0.01)

# run on final dataset
fit = sgd_lasso_logit(tX, y_label, beta0 = init_beta, stepsize = 1, discount = 0.01)