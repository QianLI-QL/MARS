#' @title Estimation of a sequence of precision matrices (SSNAL, iADMM, eADMM)
#' @description  Under the high-dimensional setting, this function is designed
#'   for generating a solution path of the precision matrices by a second-order
#'   algorithm with the l1-penalized D-trace loss.
#' @details Input a p times n sample matrix and a vector of tuning parameters to
#'   obtain a path of estimated precision matrices, where p is the sample
#'   dimension and n is the sample size. Note that our algorithm is designed for
#'   the case that p is much greater than n.
#' @param X A p times n sample matrix.
#' @param stoptol A small value to control the stop tolerance.
#' @param Lambdapath A vector of tuning parameters.
#' @param calmethod A string for choosing the computation method from "SSNAL",
#'   "iADMM", "eADMM".
#' @param maxiter A integer is set to be the maximum iteration number for
#'   the main loop.
#' @param stopmethod A string to control the stop method, could be "bigs" or
#'   "fix". If it is "bigs", then the algorithm will stop when the number of
#'   non-zero off-diagonal components exceeds n times fixnumber. If it is "fix",
#'   then the algorithm will only calculate the estimator with first fixnumber
#'   biggest tuning parameters.
#' @param printmain A boolean for controlling whether to print the main results
#'   during the calculation or not.
#' @param printsub A boolean for controlling whether to print the inner problem
#'   results during the calculation or not.
#' @param sigma The initial value of parameter sigma.
#' @param fixnumber As described in the previous.
#' @param maxlambdacheck A boolean to decide whether to narrow the lambda from
#'   its upper bound path or not.
#' @return A list of solutions: Omegapath, Lambdapath, and timepath, which are
#'   the estimated solutions, the tunning parameters, and the computation time
#'   respectively.
#' @references Qian Li, Binyan Jiang, Defeng Sun. MARS: A second-order reduction
#'   algorithm for high-dimensional sparse precision matrices estimation.
#' @export
PMEAS <- function(X = NULL, stoptol = 1e-04, Lambdapath = seq(1,0.5,-0.1), calmethod = c("SSNAL", "iADMM", "eADMM"), maxiter = 10, stopmethod = "fix", printmain = FALSE, printsub = FALSE, sigma = 1, fixnumber = 10, maxlambdacheck = TRUE) {

    # check the inputs
    if (is.null(X)) {
        stop("Error: Doesn't find X!")
    }else if (!is.matrix(X)) {
        stop("Error: Input X must be matrix form!")
    }

    if (!all(Lambdapath > 0)) {
        stop("Error: Lambdapath must all be positive!")
    }

    if (length(stoptol) != 1) {
        stop("Error: stoptol must be a single value!")
    }

    if (stoptol <= 0) {
        stop("Error: stoptol must be positive!")
    } else if (stoptol >= 1) {
        stop("Error: stoptol must be less than 1!")
    }

    if ( maxiter <= 0) {
        stop("Error: maxiter must be a positive integer")
    }

    if (match(calmethod, c("SSNAL", "iADMM", "eADMM") )== 0) {
        stop("Error: calmethod must be one of {\"SSNAL\", \"iADMM\", \"eADMM\"}!")
    }


    if (match(stopmethod, c("fix", "bigs")) == 0) {
        stop("Error: stopmethod must be one of {\"fix\", \"bigs\"}!")
    }
    if (stopmethod == "fix") {
        if (!is.integer(fixnumber)) {
            stop("Error: fixnumber must be integer")
        }else if (fixnumber <= 0) {
            stop("Error: fixnumber must be positive")
        }
    }

    if (!is.logical(printmain)) {
        stop("Error: printmain must be \"TRUE\" or \"FALSE\"!")
    }

    if (!is.logical(printsub)) {
        stop("Error: printsub must be \"TRUE\" or \"FALSE\"!")
    }


    if (!is.double(sigma) || sigma < 0 || length(sigma) != 1){
        stop("Error: sigma must be a positive single double-precision value")
    }

    if (!is.logical(maxlambdacheck)) {
        stop("Error: maxlambdacheck must be \"TRUE\" or \"FALSE\"!")
    }


    pathsolution <- PMEASc(X, stoptol, Lambdapath, calmethod, stopmethod, maxiter, printmain, printsub, sigma, fixnumber, maxlambdacheck);


    return(pathsolution)
}
