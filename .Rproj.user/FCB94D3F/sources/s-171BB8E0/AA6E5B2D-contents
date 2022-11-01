#' @title An example for running MARS
printmain <- TRUE
printsub <- FALSE
sigma <- 1
stoptol <- 1e-4
p <- 3000
n <- 50
X <- matrix(runif(p * n), p, n)
fixnumber <- ncol(X)
solution <- MARS(X, stoptol, Lambdapath = seq(1, 0.5, -0.01), maxiter = 10, stopmethod = "bigs", printmain, printsub, sigma, fixnumber, maxlambdacheck = TRUE)

