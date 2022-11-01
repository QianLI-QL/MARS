#' @title An example for running SSNAL or iADMM/eADMM
printmain <- TRUE
printsub <- FALSE
sigma <- 1
stoptol <- 1e-4
p <- 3000
n <- 50
X <- matrix(runif(p * n), p, n)
fixnumber <- ncol(X)
calmethod = "SSNAL" # you may choose another one, such as "iADMM" or "eADMM"
# note that if calmethod is "iADMM" or "eADMM", then maxiter sould be set larger, such as 1000, 3000, or even larger.
solution <- PMEAS(X, stoptol, Lambdapath = seq(1, 0.5, -0.01), calmethod, maxiter = 10, stopmethod = "bigs", printmain, printsub, sigma, fixnumber, maxlambdacheck = TRUE)

