printmain <- TRUE
printsub <- FALSE
sigma <- 1
stoptol <- 1e-4
X <- matrix(runif(3000*50), 3000, 50)
fixnumber <- ncol(X)
solution <- MARS(X, stoptol, Lambdapath = seq(1, 0.5, -0.01), maxiter = 10, stopmethod = "bigs", printmain, printsub, sigma, fixnumber, maxlambdacheck = TRUE)

