#ifndef PMEASC_H
#define PMEASC_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <iostream>

using namespace arma;

#ifndef DOUBLE_EPS
#define DOUBLE_EPS 1E-15;
#endif


Rcpp::List PMEASc(arma::mat X, double stoptol, arma::vec Lambdapath, std::string calmethod,std::string stopmethod, unsigned int maxiter, bool printyes, bool printyessub, double sigma, int numlam, bool maxlambdacheck);

void PMESSNALmainc(arma::mat A, double lambda, double stoptol, int maxiter, arma::mat &Omega, arma::mat &Y, arma::mat &Z, int p, int n, double sigma, bool printyessub, double &primobj,\
                 double &dualobj, double &gap, double &primfeas, double &dualfeas, double &eta, int &nnzOmega, arma::mat Ip);

void PMESSNCGc(arma::mat Omega, arma::mat &Y, arma::mat &YAT, arma::mat &Z, arma::mat &Ztmp, arma::mat A, arma::mat Ip, int &subbreakyes, double lambda, double sigma, int maxitersub, double Stolconst,\
      double stoptol, int p, int n);

void PMECG(arma::mat res, double tolCG, int maxiterCG, arma::mat A, arma::mat U, int p, int n, double sigma, arma::mat &direction, int &solveok, std::vector<double> &err);

void findstep(arma::mat GradPsiY, double steptol, double stepop, double sigma, arma::mat &Ztmp, arma::mat &Z, arma::mat &Y, arma::mat A, arma::mat direction, double &PsiY, double lambda, double &alp, int p);


void PMEiADMMc(arma::mat A, double lambda, double stoptol, int maxiter, arma::mat& Omega, arma::mat& Y, arma::mat& Z, int p, int n, double sigma, bool printyessub, double& primobj, \
    double& dualobj, double& gap, double& primfeas, double& dualfeas, double& eta, int& nnzOmega, arma::mat Ip);

void PMEeADMMc(arma::mat S, arma::mat A, arma::vec eigenS, arma::mat eigenU, double lambda, double stoptol, int maxiter, arma::mat& Omega, arma::mat& Y, arma::mat& Z, int p, int n, double sigma, bool printyessub, double& primobj, \
    double& dualobj, double& gap, double& primfeas, double& dualfeas, double& eta, int& nnzOmega,arma::mat Ip);

void ADMMCG(arma::mat& Y, arma::mat res, arma::mat A, double tolCG, int maxiterCG, double sigma);

int modM(int a, int b);

/* ************************************ */
arma::mat findAM(arma::mat X);

arma::vec findAandeigen(arma::mat X, arma::mat& A, arma::mat &U);

double oneoffnorm(arma::mat X);

arma::mat proxBmainM(arma::mat XX, double lambda);

void updatesigmaM(int &primwin, int &dualwin, double &sigma, int iter, int subbreakyes);


arma::mat partgradient(arma::mat Zin, double lambda);

arma::mat OperatorVp(arma::mat Ag, arma::mat U);

double findmaxlambda(arma::mat S);


template <class T>
arma::vec vectoraccess(std::vector<T> a, int b, int c);

//arma::vec vectoraccess(std::vector<double> a, int b, int c);

//arma::mat vectoraccess(std::vector<int> a, int b, int c);

#endif//
