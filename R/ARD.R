#' Automatic Relevance Determination
#'
#' \code{ARD} returns the estimated coefficients, estimated hyper-parameters, number of iterations and likelihood scores by Automated Relevance Determination (ARD).
#'
#' ARD based on three optimization algorithms -- EM algorithm, Fixed Point update rule and Reweighted-l1 algorithm.
#' If use method 'reweighted-l1', please provide sigma.square otherwise an external function was called to estimate sigma square
#' if you have prefered initial value for gamma, set gamma.init, otherwise gamma.init was initialized to be rep(1,ncol(X))
#' If use method: EM or Fixed Point update,sigma.square.init was used,the default initial value is variance of Y,if you have prefered initial value for eta, set eta.init, otherwise eta.init was initialized to be rep(1,ncol(X))
#' @param X predictor matrix of dimension \eqn{n*q}.
#' @param Y continuous outcome vector of dimension \eqn{p}.
#' @param sigma.square variance of noise, used for method 'Reweighted l1'. Default value is estimated by \code{estimateVariance(X,Y)}
#' @param sigma.square.init variance of noise, used for method 'FixedPoint' and 'EM'. Default value is \code{var(Y)}
#' initial_val initial value for \eqn{\alpha}. Default value is a vector of \eqn{0}s.
#' @param gamma.init initial value for hyper parameter, only for method 'Reweighted l1'
#' @param eta.init initial value for hyper parameter, only for method 'EM' and 'FixedPoint'
#' @param threshold_eta cut off point for very large eta
#' @param maxstep max number of iterations
#' @param margin iteration stoping creteria
#' @param verbosity print current iteraction number
#' @param compute.likelihood whether compute likelihood in each step
#' @param method chose from method "FixedPoint","EM" and "L1"
#' @return returns the estimated coefficients, estimated hyper-parameters, number of iterations and likelihood scores by Automated Relevance Determination (ARD).
#' @examples
#' set.seed(99)
#' n = 100
#' p = 200
#' x <- matrix(rnorm(2*n*p,0,1),nrow=2*n,ncol=p)
#' betas=rnorm(n = p, s = 1/exp(z_design%*%alpha))
#' y <- x%*%betas + rnorm(2*n,0,1)
#' x_train = x[1:n,]
#' y_train = y[1:n]
#' x_test = x[(n+1):(2*n),]
#' y_test = y[(n+1):(2*n)]

#' ARD.fit = ARD(x_train,y_train)
#' 1-get_mse(x_test%*%ARD_fit$coef,y_test)/var(y_test)

ARD<-function(X,Y,sigma.square = estimateVariance(X,Y),
              sigma.square.init = as.numeric(var(Y)),
              gamma.init = rep(1,ncol(X)), eta.init=rep(1,ncol(X)),
              threshold_eta = 10000,maxstep = 1000,margin=0.001,verbosity = 1,
              compute.likelihood = TRUE,method = "FixedPoint"){
        if (method == "L1"){
                result<-Reweighted_l1.ARD(X,Y,sigma.square=sigma.square,gamma.init=gamma.init,maxstep=maxstep,margin=margin,verbosity = verbosity,compute.likelihood=compute.likelihood)
                return(list(coefficients = result$coef,hyper_param = 1/result$gamma,sigma.square = sigma.square,n_iters = result$n_iter,likelihood.score = result$likelihood.score))
        }
        if (method == "FixedPoint"){
                result<-FixedPoint.ARD(X,Y,sigma.square.init = sigma.square.init, eta.init=eta.init, maxstep= maxstep,margin=margin,verbosity=verbosity,compute.likelihood = compute.likelihood)
                return(list(coefficients = result$mu_vec,hyper_param = result$eta,sigma.square = result$sigma.square,n_iters = result$n_iter,likelihood.score = result$likelihood.score))
        }
        if (method == "EM"){
                result<-EM.ARD(X,Y,sigma.square.init = sigma.square.init, eta.init=eta.init, maxstep= maxstep,margin=margin,verbosity=verbosity,compute.likelihood = compute.likelihood)
                return(list(coefficients = result$mu_vec,hyper_param = result$eta,sigma.square = result$sigma.square,n_iters = result$n_iter,likelihood.score = result$likelihood.score))
        }
}

