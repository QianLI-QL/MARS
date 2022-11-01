#ifndef MARSS_H
#define MARSS_H

#include "armadillo"
#include <math.h>
#include <iostream>

using namespace arma;
using namespace std;

#ifndef DOUBLE_EPS
#define DOUBLE_EPS 1E-15;
#endif

Rcpp::List MARSc(arma::mat X, double stoptol, arma::vec Lambdapath, std::string stopmethod,  unsigned int maxiter, bool printyes,\
                  bool printyessub, double sigma, int numlam, bool maxlambdacheck);

void PMEASmainc(arma::mat A, double lambda, double stoptol, int maxiter, arma::uvec Index, arma::umat subindex, arma::mat& Omega, arma::mat& Y, int p, int n, double sigma, bool printyessub, double& primobj,\
           double& dualobj, double& gap, double& primfeas, double& dualfeas, double& eta, int& nnzOmega);

void PMEASSSNCGc(arma::mat &Y, arma::vec &z, arma::vec &ztmp, arma::vec &SY, int &subbreakyes, arma::mat A, arma::vec x, double lambda, double sigma, int maxitersub, double Stolconst, double stoptol, \
    int p, int n, arma::uvec Index, arma::umat subindex, arma::vec a, arma::vec b, arma::vec c, arma::vec d);

void PMEASCG(arma::mat res, double tolCG, int maxiterCG, arma::mat A, arma::umat subindex, arma::vec u, int p, int n, arma::vec a, double sigma, arma::mat &direction, int &solveok, std::vector<double> &err);


void findstep(arma::mat GradPsiY, double steptol, double stepop, double sigma, arma::mat direction, arma::mat A, arma::mat &Y, arma::vec &ztmp, arma::vec &z, double &PsiY, arma::uvec Index, arma::vec a, arma::vec d, int p, arma::umat subindex, double lambda, double &alp);

arma::mat findA( arma::mat X);

double findmaxlambda(arma::mat S, arma::vec& vecdiagS);

arma::mat proxBmain(arma::mat XX, double lambda);

arma::uvec findnewJ(arma::mat residual, int k);

arma::vec operatorSY(arma::mat Y, arma::mat A, arma::umat subindex);

//vec dotmult(vec x, vec y);

//mat dotmult(mat X, mat Y);

arma::mat operatorInvLA(arma::vec x, arma::mat A, arma::umat subindex);

arma::vec prox_b(arma::vec z, double lambda, arma::vec d);

arma::vec partgradient(arma::vec zin, arma::vec c, double lambda);

void updatesigma(int &primwin, int &dualwin, double &sigma, int iter, int subbreakyes);

int mod(int a, int n);

arma::vec vecOmega(arma::mat Omega, arma::umat index);

void findcd(arma::vec& c, arma::vec& d, arma::umat subindex);
#endif

