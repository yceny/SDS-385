#################################
# rank 1 matrix approximation
###################################################################################
# objective: max uXv s.t. ||u||^2_2 <= 1, ||v||^2_2 <= 1, ||u|| <= c1, ||v|| <= c2
###################################################################################

# proximal operator, which is basically soft thresholding method
prox <- function (u, penalty){
    aa = abs(u) - penalty
    bb = cbind(rep(0,length(aa)), aa)
    r = sign(u) * apply(bb,1,max)
    return(r)
}

binary_search <- function(u,c){

    if(l2.norm(u)==0 || sum(abs(u/sqrt(sum(u^2)))) <= c) {
        return(0)
    }
    lower = 0
    # upper = max(abs(u))-1e-5
    upper = max(abs(u))
    for (i in 1:150){
        su = prox(u, (lower + upper)/2)
        if(sum(abs(su/sqrt(sum(su^2)))) < c){
            upper = (lower + upper)/2
        } else {
            lower = (lower + upper)/2
        } 
        if((upper - lower) < 1e-6) {
            return((lower + upper)/2)
        }
    }
    warning("Didn't quite converge")
    return((lower + upper)/2)
}

r1_apprx <- function(X, c1, c2, num.iteration = 1000, tol = 1e-6){
    p = dim(X)[2]
    n = dim(X)[1]
    
    # Step 1: initialize v to have l2 norm 1, using the first right singular vector of X
    # v = svd(X)$v[,1]
    # could use a random vector as the initialization
    v <- rep(1, p)
    v <- v / l2.norm(v)

    # tol1 = convergenceCriteria * n
    # tol2 = convergenceCriteria * p
    
    # Step 2: update u and v until convergence
    # initialize delta1 and delta2, how ?
    delta1 = 0.5
    delta2 = 0.5
    u.path = matrix(0, nrow = num.iteration, ncol = n)
    v.path = matrix(0, nrow = num.iteration, ncol = p)
    
    # start loop
    for (i in 1:num.iteration){
        # update u, given v
        a = X %*% v
        delta1 = binary_search(a, c1)
        proximal.u = prox(a, delta1)
        u = proximal.u / sqrt(sum(proximal.u^2))
        u.path[i,] = u
        
        
        # update v, given u
        b = t(X) %*% u
        delta2 = binary_search(b, c2)
        proximal.v = prox(b, delta2)
        v = proximal.v / sqrt(sum(proximal.v^2))
        v.path[i, ] = v
        
        if (i >= 2){
            if (sum(abs(u.path[i, ] - u.path[i,])) <= tol && sum(abs(v.path[i, ] - v.path[i,])) <= tol){
                break
            } 
        }
    }
    
    # Step 3: update d
    d = t(u) %*% X %*% v
    
    return(list(u = u,v = v,d = d, fullu = u.path, fullv = v.path, iter = i))
    
}


# Create a random matrix with rank 1
A = runif(100)
B = runif(100)
X = A %o% B
# c should be between 0 and 1
# imit of u is c1, and the largest possible is sqrt(n)
c <- 0.5
c1 <- c*sqrt(n)
c2 <- c*sqrt(p)

##########################
# rank k approximations
##########################
rk_apprx <- function(X, k, c1, c2, max_iteration = 10000, tol = 1e-6){
	if(k == 1){
		return(r1_apprx(X, c1, c2))
	}
	result = list()
	d = vector()
	Xleft = X
	for(i in 1:k){
		r1 = r1_apprx(Xleft, c1, c2)
		d = c(d, r1$d)
		Xleft = Xleft - r1$u %*% t(r1$v)* as.numeric(r1$d)
		result$u = cbind(result$u, r1$u)
		result$v = cbind(result$v, r1$v)
	}
	result$d = diag(d)
	# result$Xhat  = result$u %*% result$d %*% t(result$v)
	# result$Fdistance = sum((result$Xhat - X) ^ 2)
	return(result)
}

# Create a random matrix with rank k
k = 3
X = matrix(0, nrow = 100, ncol = 100)
for (i in 1:k) {
      A = runif(100)
      B = runif(100)
      X = X + A %o% B
}

######################################
# application to social marketing data
######################################
tweet_data = read.csv('...')
cbind(colmeans=sort(colMeans(sqrt(tweet_data[,-1])),decreasing=T))

socialMarketing <- function(k, c1, c2){
	
	tweet_data = as.matrix(tweet_data[,-1])
	#scale X: 1.scale each row to sum =1 2. scale each row with max = 1  3. scale each row to mean = 0 and same standard deviation 4. tfidf 
	scaled_tweet_data = sqrt(tweet_data)
	
	
	result = rk_apprx(scaled_tweet_data, k, c1, c2)
	
	return(result)	
}

# interpretation of the penalty parameter c1 and c2 ?
# u is a row vector, associated with person
# v is a column vector, associated with topic
sm_result5 = socialMarketing(1, c1 = 5, c2 = 5)
sm_result3 = socialMarketing(1, c1 = 3, c2 = 3)
sm_result1.5 = socialMarketing(1, c1 = 1.5, c2 = 1.5)
