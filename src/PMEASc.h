#ifndef PMEC_H
#define PMEC_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <iostream>

using namespace arma;

#ifndef DOUBLE_EPS
#define DOUBLE_EPS 1E-15;
#endif


Rcpp::List PMEASc(mat X, double stoptol, vec Lambdapath, const std::string calmethod, const std::string stopmethod, const unsigned int maxiter, const bool printyes, const bool printyessub, double sigma, int numlam, bool maxlambdacheck);

void PMESSNALmainc(mat A, double lambda, double stoptol, int maxiter, mat &Omega, mat &Y, mat &Z, int p, int n, double sigma, bool printyessub, double &primobj,\
                 double &dualobj, double &gap, double &primfeas, double &dualfeas, double &eta, int &nnzOmega, mat Ip);

void PMESSNCGc(mat Omega, mat &Y, mat &YAT, mat &Z, mat &Ztmp, mat A, mat Ip, int &subbreakyes, double lambda, double sigma, int maxitersub, double Stolconst,\
      double stoptol, int p, int n);

void PMECG(mat res, double tolCG, int maxiterCG, mat A, mat U, int p, int n, double sigma, mat &direction, int &solveok, std::vector<double> &err);

void findstep(mat GradPsiY, double steptol, double stepop, double sigma, mat &Ztmp, mat &Z, mat &Y, mat A, mat direction, double &PsiY, double lambda, double &alp, int p);


void PMEiADMMc(mat A, double lambda, double stoptol, int maxiter, mat& Omega, mat& Y, mat& Z, int p, int n, double sigma, bool printyessub, double& primobj, \
    double& dualobj, double& gap, double& primfeas, double& dualfeas, double& eta, int& nnzOmega, mat Ip);

void PMEeADMMc(mat S, mat A, vec eigenS, mat eigenU, double lambda, double stoptol, int maxiter, mat& Omega, mat& Y, mat& Z, int p, int n, double sigma, bool printyessub, double& primobj, \
    double& dualobj, double& gap, double& primfeas, double& dualfeas, double& eta, int& nnzOmega,mat Ip);

void ADMMCG(mat& Y, mat res, mat A, double tolCG, int maxiterCG, double sigma);

int modM(int a, int b);

/* ************************************ */
mat findAM(mat X);

vec findAandeigen(mat X, mat& A, mat &U);

double oneoffnorm(mat X);

mat proxBmainM(mat XX, double lambda);

void updatesigmaM(int &primwin, int &dualwin, double &sigma, int iter, int subbreakyes);


mat partgradient(mat Zin, double lambda);

mat OperatorVp(mat A, mat g, mat U);

double findmaxlambda(mat S);

template <class T>
arma::vec vectoraccess(std::vector<T> a, int b, int c);

//vec vectoraccess(std::vector<double> a, int b, int c);

//vec vectoraccess(std::vector<int> a, int b, int c);

#endif//