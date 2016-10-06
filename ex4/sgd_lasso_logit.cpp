#include <RcppEigen.h>
#include <algorithm>  // std::max

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::MatrixXi;
using Eigen::Upper;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::SparseVector;
typedef Eigen::MappedSparseMatrix<double>  MapMatd;
typedef Map<MatrixXi>  MapMati;
typedef Map<VectorXd>  MapVecd;
typedef Map<VectorXi>  MapVeci;


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::export]]
SEXP sgd_lasso_logit(MapMatd X, VectorXd Y, VectorXd beta0, double stepsize, double fudge_factor = 1e-6, double lambda = 1.0, double threshold = 1e-5, double discount = 0.01){
    // X is the design matrix stored in column-major format
    // i.e. with features for case i stores in column i
    // Y is the vector of counts
    // Thus Y[i] ~ Binomial( 1, w[i])  )
    // w[i] = 1/{1+exp(- x[i] dot beta)}
    // where Beta is the regression vector we want to estimate
    // lambda is the regularization parameter. The default is 1.0, which corresponds to lasso
    // discount is used to calculate the average of negative log likelihood
    int num_obs = X.cols();
    int num_feature = X.rows();
    SparseVector<double> x(num_feature);
    double phi, y_hat, weight, error, skip, h, new_stepsize, mu;
    double avg_nll = 0.0;
    NumericVector nll(num_obs, 0.0);
    double g0squared = 0.0;
    double grad = 0.0;
    int j, k;
    
    // keep track last time the feature i was updated (used for lazy update)
    NumericVector last_update(num_feature, 0.0);
    
    //initialize parameters
    double w = (Y.sum() + 1.0) / (num_obs + 2.0);
    double alpha = log( w /(1.0 - w));
    VectorXd beta(num_feature);
    VectorXd Gsquared(num_feature);
    for(int j=0; j<num_feature; j++) {
        Gsquared(j) = 1e-3;
        beta(j) = beta0(j);
    }
    
    
    //global counter
    k = 0;
    
    //loop over each observation
    for (int i=0; i < num_obs; i++){
        //extract non zero features
        x = X.innerVector(i);
        phi = alpha + x.dot(beta);
        y_hat = exp(phi);
        
        //Update average negative log likelihood
        avg_nll = (1.0 - discount) * avg_nll + discount * (log(1.0 + exp(phi)) - Y[i] * phi);
        nll[k] = avg_nll;
        
        //Update intercept
        error = Y[i] - y_hat;
        g0squared += error * error;
        alpha += stepsize * error / (sqrt(g0squared) + fudge_factor);
        
        //Start updating data, iterate over non zero value in the feature
        for (SparseVector<double>::InnerIterator it(x); it; ++it){
            // the index of non zero value in the feature
            j = it.index();
            
            // weighting for lasso penalty
            weight = 1.0/(1.0 + fabs(beta(j)));
            
            // Step 1: lazy update : aggregate all the penalty-only updates since the last time we updated this feature.
            skip = k - last_update(j);
            h = sqrt(Gsquared(j) + fudge_factor);
            new_stepsize = skip * stepsize / h;
            beta(j) = sgn(beta(j))*fmax(0.0, fabs(beta(j)) - new_stepsize * weight * lambda);
            
            //update last_update vector
            last_update[j] = i;
            
            //Step 2: gradient descent update (using AdaGrad)
            grad = - error * it.value();
            Gsquared(j) += grad * grad;
            h = sqrt(Gsquared(j)) + fudge_factor;
            new_stepsize = stepsize / h;
            mu = beta(j) - new_stepsize * grad;
            beta(j) = sgn(mu) * fmax(0.0, fabs(mu) - new_stepsize * weight * lambda);
        }
        k++;
    }
    
    // At the very end, apply the accumulated penalty for the variables we haven't touched recently
    for(int j=0; j< num_feature; j++) {
            skip = k - last_update(j);
            h = sqrt(Gsquared(j)) + fudge_factor;
            new_stepsize = skip * stepsize / h;
            beta(j) = sgn(beta(j))*fmax(0.0, fabs(beta(j)) - new_stepsize * weight * lambda);
    }
    
    return List::create(Named("alpha") = alpha,
                        Named("beta") = beta,
                        Named("nll_tracker") = nll);

}
