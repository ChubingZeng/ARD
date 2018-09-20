## 1. EM algorithm
EM.ARD<-function(X,Y,sigma.square.init = as.numeric(var(Y)),eta.init=rep(1,ncol(X)), maxstep=300,margin=0.001,verbosity=0,compute.likelihood = FALSE){
        n=nrow(X)
        p = ncol(X)
        # step 0: initialize eta
        eta = eta.init
        sigma.square = sigma.square.init
        k = 1
        #delta = Inf
        likelihood.score <- c()

        while(k < maxstep){
                # step 1: E step, compute diag(Sigma) and mu
                Rinv <- backsolve(chol(crossprod(X)/sigma.square + diag(c(eta))), diag(1, p))
                diagSigma <- rowSums(Rinv^2)
                mu_vec <- (Rinv %*% (crossprod(Rinv, crossprod(X,Y))))/sigma.square

                # Compute likelihood
                if (compute.likelihood == TRUE){
                        likelihood.score = c(likelihood.score,compute_likelihood(X,Y,eta,sigma.square))
                }

                if(k > 1){
                        if(sum(abs(mu_vec.old - mu_vec)) < margin ){
                                break
                        }
                }
                mu_vec.old <- mu_vec
                # step 2: M step, update eta vector
                eta = 1/(mu_vec^2 + diagSigma)
                sigma.square = (sum((Y - X%*%mu_vec)^2) + sigma.square*sum(1-eta*diagSigma))/n

                if(verbosity == 0){
                        print(paste("The ",k,"th iterations"))
                }

                k <- k+1
        }
        return(list(eta = eta,mu_vec = mu_vec,sigma.square = sigma.square,likelihood.score = likelihood.score,n_iter = k-1))
}

## 2. Fixed Point update algorithm
FixedPoint.ARD<-function(X,Y,sigma.square.init,eta.init,threshold_eta = 10000, maxstep=300,margin=0.001,verbosity=0,compute.likelihood = FALSE){
        n=nrow(X)
        p = ncol(X)
        # step 0: initialize eta
        eta = eta.init
        sigma.square = sigma.square.init
        k = 1
        likelihood.score <- c()
        #delta = Inf
        keep_eta <- 1:p
        mu_vec = rep(0,p)
        while(k < maxstep){
                # compute diag(Sigma) and mu
                Rinv <- backsolve(chol(crossprod(X[,keep_eta])/sigma.square + diag(c(eta[keep_eta]))), diag(1, ncol(X[,keep_eta])))
                diagSigma <- rowSums(Rinv^2)
                mu_vec[keep_eta] = (Rinv %*% (crossprod(Rinv, crossprod(X[,keep_eta],Y))))/sigma.square

                # Compute likelihood
                if (compute.likelihood == TRUE){
                        likelihood.score = c(likelihood.score,compute_likelihood(X,Y,eta,sigma.square))
                }

                if(k > 1){
                        if(sum(abs(mu_vec.old - mu_vec)) < margin ){
                                break
                        }
                }
                mu_vec.old <- mu_vec

                # step 2: update eta vector
                gamma = 1-eta[keep_eta]*diagSigma
                eta[keep_eta] = gamma/mu_vec[keep_eta]^2
                sigma.square = sum((Y - X%*%mu_vec)^2)/(n - sum(gamma))


                # prune the coefficients with a precision over a threshold
                keep_eta = which(eta < threshold_eta)
                mu_vec[-keep_eta] <- 0

                if(verbosity == 0){
                        print(paste("The ",k,"th iterations"))
                }

                k <- k+1
        }
        return(list(eta = eta,mu_vec = mu_vec,sigma.square=sigma.square,likelihood.score = likelihood.score,n_iter=k-1))
}

## 3. reweighted l1 algorithm
Reweighted_l1.ARD <- function(X,Y,sigma.square,gamma.init,maxstep,margin,verbosity,compute.likelihood){
        n = nrow(X);p=ncol(X)
        ## Initialize
        gamma.old = gamma.init
        k = 1
        likelihood.score <- c()
        while(k < maxstep){
                Sigma_y = sigma.square * diag(n) + X %*% diag(c(gamma.old)) %*% t(X)
                theta = diag(t(X) %*% solve(Sigma_y,X))

                ## update gamma given theta
                # deltaHat <- Variable(p)
                # obj <- sum((Y - X %*% deltaHat)^2) + 2*sigma.square*p_norm(sqrt(theta)*deltaHat, 1)
                # prob <- Problem(Minimize(obj))
                # result <- solve(prob)
                # delta.est <- result$getValue(deltaHat)
                # gamma.new <- 1/sqrt(theta)*abs(delta.est)
                #
                coef(glmnet(X,Y,alpha = 1, lambda = sigma.square/n*sum(sqrt(theta))/p), penalty.factor = sqrt(theta))
                delta.est <-  coef(glmnet(X,Y,alpha = 1, lambda = sigma.square/n*sum(sqrt(theta))/p), penalty.factor = sqrt(theta),intercept = F)[-1]
                gamma.new <- 1/sqrt(theta)*abs(delta.est)

                if(sum(abs(gamma.new - gamma.old)) < margin ){
                        break
                }
                gamma.old <- gamma.new

                # Compute likelihood
                if (compute.likelihood == TRUE){
                        likelihood.score = c(likelihood.score,compute_likelihood(X,Y,1/gamma.old,sigma.square))
                }

                if(verbosity == 0){
                        print(paste("The ",k,"th iterations"))
                }

                k <- k+1
        }
        coef = diag(c(gamma.old))%*%t(X)%*%solve(Sigma_y,Y)
        return(list(coef = coef,gamma = gamma.old,n_iter = k-1,likelihood.score = likelihood.score))
}

compute_likelihood <- function(X,Y,eta,sigma.square) {
        eta[eta>1e5] <- 1e5
        n = nrow(X)
        K = sigma.square * diag(n) + X %*% diag(c(eta)) %*% t(X)
        logdetK = determinant(K)$modulus[1]
        part1 = t(Y) %*% solve(K,Y)
        normapprox = 1/2 * (part1 + logdetK)
        return(as.numeric(normapprox))
}
