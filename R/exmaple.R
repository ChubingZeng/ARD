library(ARD)
library(classo)
library(selectiveInference)
set.seed(9)
n = 100
p = 600

x <- matrix(rnorm(2*n*p,0,1),nrow=2*n,ncol=p)
betas=rnorm(n = p, s = 2)
y <- x%*%betas + rnorm(2*n,0,1)
x_train = x[1:n,]
y_train = y[1:n]
x_test = x[(n+1):(2*n),]
y_test = y[(n+1):(2*n)]
ARD.fit.FP = ARD(x_train,y_train,method = "FixedPoint")
ARD.fit.EM = ARD(x_train,y_train,method = "EM")
ARD.fit.L1 = ARD(x_train,y_train,method = "L1")
plot(ARD.fit.FP$hyper_param,ARD.fit.EM$hyper_param)
plot(ARD.fit.FP$hyper_param,ARD.fit.L1$hyper_param)

1-get_mse(x_test%*%ARD.fit.FP$coef,y_test)/var(y_test)
1-get_mse(x_test%*%ARD.fit.EM$coef,y_test)/var(y_test)
1-get_mse(x_test%*%ARD.fit.L1$coef,y_test)/var(y_test)

plot.ARD(ARD.fit.EM)[[1]]
plot.ARD(ARD.fit.FP)[[1]]
plot.ARD(ARD.fit.L1)[[1]]
