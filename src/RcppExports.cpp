// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// MARSc
Rcpp::List MARSc(arma::mat X, double stoptol, arma::vec Lambdapath, std::string stopmethod, unsigned int maxiter, bool printyes, bool printyessub, double sigma, int numlam, bool maxlambdacheck);
RcppExport SEXP _MARS_MARSc(SEXP XSEXP, SEXP stoptolSEXP, SEXP LambdapathSEXP, SEXP stopmethodSEXP, SEXP maxiterSEXP, SEXP printyesSEXP, SEXP printyessubSEXP, SEXP sigmaSEXP, SEXP numlamSEXP, SEXP maxlambdacheckSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type stoptol(stoptolSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Lambdapath(LambdapathSEXP);
    Rcpp::traits::input_parameter< std::string >::type stopmethod(stopmethodSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< bool >::type printyes(printyesSEXP);
    Rcpp::traits::input_parameter< bool >::type printyessub(printyessubSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type numlam(numlamSEXP);
    Rcpp::traits::input_parameter< bool >::type maxlambdacheck(maxlambdacheckSEXP);
    rcpp_result_gen = Rcpp::wrap(MARSc(X, stoptol, Lambdapath, stopmethod, maxiter, printyes, printyessub, sigma, numlam, maxlambdacheck));
    return rcpp_result_gen;
END_RCPP
}
// PMEASc
Rcpp::List PMEASc(arma::mat X, double stoptol, arma::vec Lambdapath, std::string calmethod, std::string stopmethod, unsigned int maxiter, bool printyes, bool printyessub, double sigma, int numlam, bool maxlambdacheck);
RcppExport SEXP _MARS_PMEASc(SEXP XSEXP, SEXP stoptolSEXP, SEXP LambdapathSEXP, SEXP calmethodSEXP, SEXP stopmethodSEXP, SEXP maxiterSEXP, SEXP printyesSEXP, SEXP printyessubSEXP, SEXP sigmaSEXP, SEXP numlamSEXP, SEXP maxlambdacheckSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type stoptol(stoptolSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Lambdapath(LambdapathSEXP);
    Rcpp::traits::input_parameter< std::string >::type calmethod(calmethodSEXP);
    Rcpp::traits::input_parameter< std::string >::type stopmethod(stopmethodSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< bool >::type printyes(printyesSEXP);
    Rcpp::traits::input_parameter< bool >::type printyessub(printyessubSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type numlam(numlamSEXP);
    Rcpp::traits::input_parameter< bool >::type maxlambdacheck(maxlambdacheckSEXP);
    rcpp_result_gen = Rcpp::wrap(PMEASc(X, stoptol, Lambdapath, calmethod, stopmethod, maxiter, printyes, printyessub, sigma, numlam, maxlambdacheck));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MARS_MARSc", (DL_FUNC) &_MARS_MARSc, 10},
    {"_MARS_PMEASc", (DL_FUNC) &_MARS_PMEASc, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_MARS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
