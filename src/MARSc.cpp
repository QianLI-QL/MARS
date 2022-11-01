// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
//#define ARMA_64BIT_WORD 1

#include <math.h>
#include <iostream>
#include "MARSc.h"

using namespace std;
using namespace arma;
using namespace Rcpp;

#ifndef DOUBLE_EPS
#define DOUBLE_EPS 1E-15;
#endif


/* main function:
 * using adaptive sieving strategy find a solution path
 */
//' @export
//'
// [[Rcpp::export]]
Rcpp::List MARSc(arma::mat X, double stoptol, arma::vec Lambdapath, std::string stopmethod,  unsigned int maxiter,  bool printyes,  bool printyessub, double sigma, int numlam, bool maxlambdacheck){
  // set some given information
  stoptol = max(stoptol, 1e-6);
  int p = X.n_rows;
  int n = X.n_cols;
  arma::mat A = findA(X);
  double primobj=0, dualobj=0, gap=0, primfeas=0, dualfeas=0, eta=0, etaorg=0;
  int nnzOmega = 0;
  cout.precision(4);
  arma::vec vecdiagS(p);


  // treat the lambda path
  double maxlambda = findmaxlambda(A, vecdiagS);
  if (printyes){
      cout << "\n Max lambda = " << fixed << maxlambda << endl;
  }

  if (maxlambdacheck){
    arma::uvec clamplambdaindex = find(Lambdapath < maxlambda);
    Lambdapath = Lambdapath(clamplambdaindex);
    if (!Lambdapath.is_sorted("descend")){
      Lambdapath = sort(Lambdapath,"descend");
    }
  }

  int pathlength = Lambdapath.n_elem;
  field<sp_mat> Omegapath(pathlength,1);
  arma::vec timepath(pathlength);

  // set initial values
  arma::mat Omega(p, p);
  arma::mat Y = A;
  arma::uvec Index(p);
  for (int i = 0; i < p; i++){
    Omega(i, i) = 1 / vecdiagS(i);
    Y.row(i) /= vecdiagS(i);
    Index(i) = i * p + i;
  }

  arma::umat subindex = arma::ind2sub(size(Omega), Index);
  wall_clock pathtime;
  pathtime.tic();

  // main iteration
  for (int iterpath = 0; iterpath < pathlength; iterpath++)
  {
    wall_clock itertime;
    itertime.tic();
    double lambda = Lambdapath(iterpath);
    double timeiter;
    if (printyes) {
      cout << "*********************************************************************" << endl;
      cout << "lambda = " << fixed << lambda << endl;
    }
    if (iterpath == 0){
      PMEASmainc(A, lambda, stoptol, maxiter, Index, subindex, Omega, Y, p, n, sigma, printyessub, primobj,\
                 dualobj, gap, primfeas, dualfeas, eta, nnzOmega);
      if (printyessub){
        cout << "subproblem results printing end." << endl;
        cout << "---------------------------------------" << endl;
      }
    }

    arma::mat gradP = 0.5 * ( (Omega * A) * A.t() + A * (A.t() * Omega) )- eye(p, p);
    arma::mat residual = gradP + proxBmain(Omega - gradP, lambda);
    etaorg = norm(residual, "fro") / (1 + norm(gradP, "fro") + norm(Omega, "fro"));
    int numAS = 0;
    if (printyes){
      if (etaorg < stoptol){
        cout << "etaorg = " << scientific << etaorg << ", which already < stoptol, solution is as same as the previous one " << endl;
      }else{
        cout << " Iter|   etaorg   |  primobjsub   dualobjsub     gapsub   | primfeassub  dualfeassub  etasub   | time   nnzoff" << endl;
      }
      if (iterpath == 0){
        timeiter = itertime.toc();
        numAS = 1;
        cout << "   " << dec << numAS << " | " << scientific <<  etaorg << " | " << primobj << "  " << dualobj << "  " << gap \
             << " |  " << primfeas << "  " << dualfeas << "  " << eta << " | " << fixed << timeiter << "   " << dec << nnzOmega << endl;
      }
    }

    while (etaorg > stoptol && numAS < 10){
      numAS++;
      int nonzeronumberIndex = Index.n_elem;
      arma::mat abstriures = abs(trimatu(residual));
      arma::uvec newJ = find(abstriures > 1e-5);
      arma::uvec newIndex = join_cols(Index, newJ);
      newIndex = unique(newIndex);
      int nonzeronumbernewIndex = newIndex.n_elem;
      if (nonzeronumbernewIndex > 5 * nonzeronumberIndex && nonzeronumberIndex != 0){
        arma::uvec sortindex = sort_index(abstriures, "descend");
        newJ = sortindex(span(0, nonzeronumberIndex));
        Index = join_cols(Index, newJ);
        Index = unique(Index);
      }else{
        Index = newIndex;
      }
      subindex = arma::ind2sub(size(Omega), Index);
      PMEASmainc(A, lambda, stoptol, maxiter, Index, subindex, Omega, Y, p, n, sigma, printyessub, primobj,\
                 dualobj, gap, primfeas, dualfeas, eta, nnzOmega);
      if (printyessub){
      cout << " subproblem results printing end." << endl;
      cout << " ---------------------------------------" << endl;
      }
      gradP = 0.5 * ( (Omega * A) * A.t() + A * (A.t() * Omega) )- eye(p, p);
      residual = gradP + proxBmain(Omega - gradP, lambda);
      etaorg = norm(residual, "fro") / (1 + norm(gradP, "fro") + norm(Omega, "fro"));
      if (printyes){
        if (printyessub){
           cout << " Iter|   etaorg   | primobjsub  dualobjsub    gapsub   | primfeassub  dualfeassub  etasub   | time   nnzoff" << endl;
        }
        timeiter = itertime.toc();
        cout << "   " << dec << numAS << " | " << scientific <<  etaorg << " | " << primobj << "  " << dualobj << "  " << gap \
             << " |  " << primfeas << "  " << dualfeas << "  " << eta << " | " << fixed << timeiter << "   " << dec << nnzOmega << endl;
      }
    }// end of while
    timeiter = itertime.toc();
    Omegapath(iterpath) = sp_mat(Omega);
    timepath(iterpath) = timeiter;

    if (printyes) {
      cout << " total time in this calculation process: " << fixed << timeiter << "(s)" << endl;
    }
    if (stopmethod == "bigs"){
      if (nnzOmega > (numlam * n) || numAS > 4){
        if (printyes){
          cout << " lambda small enough" << endl;
        }

        if (pathlength > 1 && iterpath < pathlength - 1){
          if (iterpath == pathlength - 2){
            Omegapath = Omegapath.rows(0, iterpath);
            Lambdapath.shed_rows(pathlength - 1, pathlength - 1);
            timepath.shed_rows(pathlength - 1, pathlength - 1);
          }else{
            Omegapath = Omegapath.rows(0, iterpath);
            Lambdapath.shed_rows((iterpath + 1), (pathlength - 1));
            timepath.shed_rows((iterpath + 1), (pathlength - 1));
          }
        }
        break;
      }
    }
    if(stopmethod == "fix" && iterpath == numlam - 1){
      if (pathlength > 1 && iterpath < pathlength - 1){
        if (iterpath == pathlength - 2){
          Omegapath = Omegapath.rows(0, iterpath);
          Lambdapath.shed_rows(pathlength - 1, pathlength - 1);
          timepath.shed_rows(pathlength - 1, pathlength - 1);
        }else{
          Omegapath = Omegapath.rows(0, iterpath);
          Lambdapath.shed_rows((iterpath + 1), (pathlength - 1));
          timepath.shed_rows((iterpath + 1), (pathlength - 1));
        }
      }
      break;
    }


    // end of stopmethod checking

  }// end of for loop

  if (printyes){
    double timepath = pathtime.toc();
    cout << "\n total time for calculating solution path: " << fixed << timepath << "(s)" << endl;
    cout << "\n *********************************************************************" << endl;
  }

  return Rcpp::List::create(Named("Omegapath") = Omegapath,
                            Named("Lambdapath") = Lambdapath,
                            Named("timepath") = timepath);
}


/* ***********************************
* main function
*/
void PMEASmainc(arma::mat A, double lambda, double stoptol, int maxiter, arma::uvec Index, arma::umat subindex, arma::mat &Omega, arma::mat &Y, int p, int n, double sigma, bool printyessub, double &primobj,\
                double  &dualobj, double &gap, double &primfeas, double &dualfeas, double &eta, int &nnzOmega){
  // preparation

  int lengthIndex = Index.n_elem;
  int maxitersub = 10;
  double Stolconst = 0.5;
  int primwin = 0, dualwin = 0;
  int iter = 0;
  bool breakyes = false;
  arma::vec c, d;
  findcd(c, d, subindex);
  arma::vec b = d + c, a = 0.25 * d + c;

  // set initial value
  arma::vec vomega = vecOmega(Omega, subindex);
  arma::vec x = vomega % b;
  arma::vec SY = 0.5 * operatorSY(Y, A, subindex);
  arma::vec SYc = SY - c;
  arma::mat OA = operatorInvLA(x, A, subindex);

  primfeas = norm(Y - OA, "fro") / (1 + norm(Y, "fro"));
  dualfeas = norm(SYc, "fro") / (1 + sqrt(p));
  primobj = 0.5 * pow(norm(OA, "fro"), 2) - accu((vomega % c)) + lambda * norm((vomega % d), 1);
  dualobj = - 0.5 * pow(norm(Y, "fro"), 2);
  gap = (primobj - dualobj) / (1 + fabs(primobj) + fabs(dualobj));
  eta = norm(SYc + prox_b(x - SYc, lambda, d), "fro") / (1 + norm(SYc, "fro") + norm(x, "fro"));

  if ((max(primfeas, dualfeas) < 500 * max(1e-6, stoptol)) && (eta < stoptol)){
    breakyes = true;
  }

  if (printyessub){
    cout << " ---------------------------------------" << endl;
    cout << " subproblem results printing start:" << endl;
    cout << " iter|  objprimal    objdual      gap    |  primfeas    dualfeas      eta       sigma  " << endl;
    cout << "   " << dec << iter << " | " << scientific << primobj << "  " << dualobj << " "<< gap << " | " \
         << primfeas << " " << dualfeas << " " << eta << " " << sigma << endl;
  }

  // start the main iteration of main function
  while (iter < maxiter && breakyes == false){
    iter++;
    if (dualfeas < 1e-3) {
      maxitersub = max(maxitersub, 30);
    } else if (dualfeas < 1e-1){
      maxitersub = max(maxitersub, 20);
    }

    int subbreakyes = 0;
    arma::vec z(lengthIndex), ztmp(lengthIndex);
    PMEASSSNCGc(Y, z, ztmp, SY, subbreakyes, A, x, lambda, sigma, maxitersub, Stolconst, stoptol, \
    p, n, Index, subindex, a, b, c, d);
    SYc = SY - c;
    x = sigma * ztmp;
    vomega = (x % a);

    OA = operatorInvLA(x, A, subindex);
    primfeas = norm(Y - OA, "fro")/(1 + norm(Y, "fro"));
    dualfeas = norm(SYc + z, "fro")/(1 + sqrt(p));
    primobj = 0.5 * pow(norm(OA, "fro"), 2) - accu((vomega % c)) + lambda * norm((vomega % d), 1);
    dualobj = - 0.5 * pow(norm(Y, "fro"), 2);
    gap = (primobj - dualobj)/(1 + fabs(primobj) + fabs(dualobj));
    eta = norm(SYc + prox_b(x - SYc, lambda, d), "fro") / (1 + norm(SYc, "fro") + norm(x, "fro"));

    if (max(primfeas, dualfeas) < 500 * max(1e-6, stoptol) && eta < stoptol){
       breakyes = true;
    }

    if (primfeas < dualfeas){
      primwin = primwin + 1;
    }
    else{
      dualwin = dualwin + 1;
    }

    updatesigma(primwin, dualwin, sigma, iter, subbreakyes);
    if (printyessub){
      cout << "   " << dec << iter << " | " << scientific << primobj << "  " << dualobj << " "<< gap << " | " \
           << primfeas << " " << dualfeas << " " << eta << " " << sigma << endl;
    }
  } // end of while
  nnzOmega = 0;
  Omega = sp_mat(zeros(p, p));
  for (int oidx = 0; oidx < lengthIndex; oidx++){
    int orow = subindex(0, oidx);
    int ocol = subindex(1, oidx);
    double elemo = vomega(oidx);
    if (abs(elemo) > (DOUBLE_EPS )){
      nnzOmega++;
      if (orow == ocol){
        Omega(orow, ocol) = elemo;
      }else{
        Omega(orow, ocol) = elemo;
        Omega(ocol, orow) = elemo;
      }
    }
  }
  nnzOmega = 2 * (nnzOmega - p);
  return;
}


/*
* main part of using semismooth Newton
*/
void PMEASSSNCGc(arma::mat &Y, arma::vec &z, arma::vec &ztmp, arma::vec &SY, int &subbreakyes, arma::mat A, arma::vec x, double lambda, double sigma, int maxitersub, double Stolconst, double stoptol, \
    int p, int n, arma::uvec Index, arma::umat subindex, arma::vec a, arma::vec b, arma::vec c, arma::vec d){
  int maxiterCG = 500;
  arma::vec zin = x / sigma - SY + c;
  z = prox_b(zin, lambda, d);
  ztmp = zin - z;
  double PsiY = 0.5 * pow(norm(Y, "fro"), 2) + 0.5 * sigma * pow(norm(ztmp, "fro"), 2);

  vector<double> SSNCGsubprimfeas, SSNCGsubdualfeas; // SSNCGPsiY; // clear it at the end;
  //vector<double> SSNCGCGiter, SSNCGsolveok;
  for (int itersub = 0; itersub < maxitersub; itersub++){
    arma::mat GradPsiY = sigma * operatorInvLA(ztmp, A, subindex) - Y;
    double subprimfeas = norm(GradPsiY, "fro") / (1 + norm(Y, "fro"));
    double subdualfeas = norm(SY + z - c, "fro") / (1 + sqrt(p));
    double primratio = 0, dualratio = 0;

    double tolsubconst;
    if (max(subprimfeas, subdualfeas) < stoptol) {
      tolsubconst = 0.9;
    }else{
      tolsubconst = 0.05;
    }

    double stoptolsub = max(min(1.0, Stolconst * subdualfeas), stoptol * tolsubconst);
    SSNCGsubprimfeas.push_back(subprimfeas);
    SSNCGsubdualfeas.push_back(subdualfeas);
    //SSNCGPsiY.push_back(PsiY);

    if (subprimfeas < stoptolsub && itersub > 0)
    {
      subbreakyes = -1;
      break;
    }

    if (subdualfeas > 1e-3 || itersub <= 5)
    {
      maxiterCG = min(maxiterCG, 200);
    }else if (subdualfeas > 1e-4)
    {
      maxiterCG = min(maxiterCG, 300);
    }else if (subdualfeas > 1e-5)
    {
      maxiterCG = min(maxiterCG, 400);
    }else if (subdualfeas > 5e-6)
    {
      maxiterCG = min(maxiterCG, 500);
    }

    if (itersub > 0)
    {
      primratio = subprimfeas/SSNCGsubprimfeas[itersub - 1];
      dualratio = subdualfeas/SSNCGsubdualfeas[itersub - 1];
    }

    arma::mat res = GradPsiY;
    double tolCG = min(5e-3, 0.1 * norm(res, "fro"));
    double tolCGconst = 1;
    if(itersub > 0 && (subprimfeas > 0.1 * SSNCGsubprimfeas[0] || primratio > 0.5)){
      tolCGconst = 0.5;
    }
    if(dualratio > 1.1){
      tolCGconst = 0.5;
    }
    tolCG = tolCGconst * tolCG;

    arma::mat direction = zeros(p, n);
    int solveok = 1;
    vector<double> err;
    arma::vec u = partgradient(ztmp + z, c, lambda); // this may put into CG function
    PMEASCG(res, tolCG, maxiterCG, A, subindex, u, p, n, a, sigma, direction, solveok, err);
    //SSNCGsolveok.push_back(solveok);
    //SSNCGCGiter.push_back(err.size() - 1);
    double steptol = 1e-5;
    double stepop;
    if (itersub < 2 || (itersub <= 2 && subdualfeas > 1e-4)){
      stepop = 1;
    }else
    {
      stepop = 2;
    }
    double alp = 1;
    findstep(GradPsiY, steptol, stepop, sigma, direction, A, Y, ztmp, z, PsiY, Index, a, d, p, subindex, lambda, alp);
    SY = 0.5 * operatorSY(Y, A, subindex);
    if ((alp) < (DOUBLE_EPS)){
      subbreakyes = 11;
    }
    err.clear();


  }// end of for loop


  return;
}


/* CG */
void PMEASCG(arma::mat res, double tolCG, int maxiterCG, arma::mat A, arma::umat subindex, arma::vec u, int p, int n, arma::vec a, double sigma, arma::mat &direction, int &solveok, vector<double> &err){
  int stagnatecheck = 20;
  err.push_back(norm(res, "fro"));
  arma::mat g = res, Vg;
  double rz1 = pow(err[0], 2);
  double rz2 = 1;

  double beta, alpha, denom, residual;
  arma::vec ldua;

  for (int iterCG = 0; iterCG < maxiterCG; iterCG++){
    if (iterCG > 0){
      beta = rz1 / rz2;
      g = res + beta * g;
    }
    ldua = 0.5 * operatorSY(g, A, subindex);
    ldua = ldua % u;
    Vg = g + sigma * operatorInvLA(ldua, A, subindex);
    denom = accu(g % Vg);
    if ((fabs(denom)) < (DOUBLE_EPS)){
      solveok = 2;
      break;
    }else
    {
      alpha = rz1 / denom;
      direction = direction + alpha * g;
      res = res - alpha * Vg;
    }
    //z = res;
    residual = norm(res, "fro");
    err.push_back(residual);
    if (residual < tolCG) {break;}
    rz2 = rz1;
    rz1 = accu(res % res);
    if (iterCG > stagnatecheck){
      arma::vec ratio(10);
      for (int i = 0; i < 10; i++)
      {
        ratio(i) = err[iterCG - i + 1] / err[iterCG - i];
      }
      if (min(ratio) > 0.997 && max(ratio) < 1.003){
        solveok = -1;
        break;
      }
    }

  } // end of for loop
  return;
}







/* ****************************
 * find step length
 */
void findstep(arma::mat GradPsiY, double steptol, double stepop, double sigma, arma::mat direction, arma::mat A, arma::mat &Y, arma::vec &ztmp, arma::vec &z, double &PsiY, arma::uvec Index, arma::vec a, arma::vec d, int p, arma::umat subindex, double lambda, double &alp){
  int maxiterstep = (int) (ceil(log(1 / (steptol + DOUBLE_EPS)) / log(2)));
  double c1 = 1e-4, c2 = 0.9;
  double change0 = accu(-GradPsiY % direction);
  if (change0 >= 0) {
    //cout << "warning" << endl;
    return;
  }
  double alpconst = 0.5, LB = 0, UB = 1, PsiYold = PsiY, change1;
  arma::mat Yold = Y, gradstep1;
  arma::vec zold = z, ztmpold = ztmp, zin;
  double gLB = change0, gUB = gLB;
  for (int iterstep = 0; iterstep < maxiterstep; iterstep++){
    if (iterstep > 0){
      alp = alpconst * (LB + UB);
    }
    Y = Yold + alp * direction;
    zin = ztmpold + zold - 0.5 * alp * operatorSY(direction, A, subindex);
    z = prox_b(zin, lambda, d);
    ztmp = zin - z;
    PsiY = 0.5 * pow(norm(Y, "fro"), 2) + 0.5 * sigma * pow(norm(ztmp, "fro"), 2);
    gradstep1 = Y - sigma * operatorInvLA(ztmp, A, subindex);
    change1 = accu(gradstep1 % (alp * direction));
    if (iterstep == 0){
      gUB = change1;
      if (sign(gLB) * sign(gUB) > 0){
        break;
      }
    }
    if(fabs(change1) < c2 * fabs(change0) && PsiY - PsiYold - c1 * alp * change0 <= 1e-8 / max((double)1, fabs(PsiYold))){
      if (stepop == 1 || (stepop == 2 && fabs(change1) < steptol))
      {
        if (iterstep == maxiterstep){
          cout << "\n Warning maximum iter in finding step length reached" << endl;
        }
        break;
      }
    }
    if (sign(change1) * sign(gUB) < 0){
      LB = alp;
      gLB = change1;
    }else if (sign(change1) * sign(gLB) < 0){
      UB = alp;
      gUB = change1;
    }
  }// end of for loop;
}




/* ****************************************************************************************************** */

/* From given data X to find the matrix A */

arma::mat findA(arma::mat X){
  int p = X.n_rows;
  int n = X.n_cols;
  double rowmeans;
  arma::mat centralizedX = X, U, V;
  arma::mat A = zeros(p, n);
  arma::vec s;

  // centralize
  for(int i = 0; i < p; i++){
    rowmeans = mean(X.row(i));
    centralizedX.row(i) -= rowmeans;
  }
  svd_econ(U, s, V, centralizedX,"left");
  double sqrtn1 = sqrt(n - 1);
  for (int j = 0; j < n; j++){
    A.col(j) = U.col(j) * s(j) / sqrtn1;
  }
  return A;
}


// find the maximum lambda

double findmaxlambda(arma::mat A, arma::vec& vecdiagS){
  arma::mat S = A * A.t();
  int p = A.n_rows;
  double elementS, maxlambda = 0;
  for(int i = 0; i < p; i++){
    vecdiagS(i) = S(i, i);
    for(int j = i + 1; j < p; j++){
      elementS = 0.5 * abs(S(i,j)/S(i,i) + S(i,j)/S(j,j));
      if (elementS > maxlambda){
        maxlambda = elementS;
      }
    }
  }
  return maxlambda;
}

// original indicator function (Checked)

arma::mat proxBmain(arma::mat XX, double lambda){
  int p = XX.n_rows;
  arma::mat  XXP = zeros(p, p);
  double elem;

  for (int i = 0; i < (p - 1); i++){
    for (int j = i + 1; j < p; j++){
      elem = max(-lambda, min(lambda, XX(i,j)));
      XXP(i,j) = elem;
      XXP(j,i) = elem;
    }
  }
  return XXP;
}


// find the largerst k elements of a given matrix

arma::uvec findnewJ(arma::mat residual, int k){
  arma::uvec newJ;
  residual = trimatu(residual, 1);
  newJ = sort_index(residual, "descend");
  newJ = newJ.in_range(span(0,k));
  return newJ;
}

// calculate operator SY (checked)
arma::vec operatorSY(arma::mat Y, arma::mat A, arma::umat subindex){
  int i, j;
  //double v1, v2;

  int t = subindex.n_cols;
  //int n = A.n_cols;
  arma::vec v = zeros(t);
  for (int k = 0; k < t; k++){
    i = subindex(0, k);
    j = subindex(1, k);
    /*v1 = 0;
    v2 = 0;
    for (int l = 0; l < n; l++) {
      v1 += Y(i,l) * A(j,l);
      v2 += A(i,l) * Y(j,l);
    }*/
    double v1 = sum(Y.row(i) % A.row(j));
    double v2 = sum(A.row(i) % Y.row(j));
    v(k) = v1 + v2;
  }
  return v;
}



// calculate the generalized inverse of L (checked)

arma::mat operatorInvLA(arma::vec x, arma::mat A, arma::umat subindex){
  int t = subindex.n_cols;
  int n = A.n_cols;
  int p = A.n_rows;
  arma::mat OA = zeros(p, n);
  double omek;
  int i, j, k;

  for (k = 0; k < t; k++){
    i = subindex(0, k);
    j = subindex(1, k);
    omek = x(k);
    if (i == j){
      OA.row(i) += omek * A.row(i);
    }else{
      OA.row(i) += 0.5 * omek * A.row(j);
      OA.row(j) += 0.5 * omek * A.row(i);
    }
  }
  return OA;
}


// prox_delta_b
arma::vec prox_b(arma::vec z, double lambda, arma::vec d){
  int t = z.n_elem;
  arma::vec v = zeros(t);
  for (int i = 0; i < t; i++){
      if ((abs(d(i))) > (DOUBLE_EPS)){
      v(i) = max(-lambda, min(lambda, z(i)));
    }
  }
  return v;
}


// update sigma
void updatesigma(int &primwin, int &dualwin, double &sigma, int iter, int subbreakyes){
  double minsig = 1e-4;
  double maxsig = 1e7;
  int sigiter;
  if (iter < 10){
    sigiter = 2;
  }else if (iter < 200)
  {
    sigiter = 3;
  }else if (iter < 500){
    sigiter = 10;
  }else{
    sigiter = 20;
  }
  //double dual = (double)dualwin, prim = (double)primwin;
  if (mod(iter, sigiter) == 0 && subbreakyes < 0){
    double mult = 5;
    if ((double)primwin > max(1.0, 1.2 * dualwin)){
      primwin = 0;
      sigma = min(maxsig, mult * sigma);
    }else if ((double)dualwin > max(1.0, 1.2 * primwin)){
      dualwin = 0;
      sigma = max(minsig, sigma / mult);
    }
  }
  if (subbreakyes >= 0 && iter > 10){
    sigma = max(minsig, sigma / 2.5);
  }
  return;
}


// calculate hessian
arma::vec partgradient(arma::vec zin, arma::vec c, double lambda){
  int t = zin.n_elem;
  arma::vec u = ones(t);
  for (int k = 0; k < t; k++){
      if (((fabs(c[k])) <= (DOUBLE_EPS)) && (fabs(zin[k]) <= lambda)){
      u[k] = 0;
    }
    // this function may improve by change the way of updating 1
  }
  return u;
}


int mod(int a, int n)
{
  return a - (int)(floor(a / n) * n);
}



arma::vec vecOmega(arma::mat Omega, arma::umat subindex) {
  int t = subindex.n_cols;
  arma::vec vomega(t);
  int i, j, k;

  for (k = 0; k < t; k++){
    i = subindex(0, k);
    j = subindex(1, k);
    vomega(k) = Omega(i, j);
  }
  return(vomega);
}



void findcd(arma::vec& c, arma::vec& d, arma::umat subindex) {
  int t = subindex.n_cols;
  c = zeros(t);
  d = ones(t);
  d *= 2;
  for(int i = 0; i < t; i++){
    if(subindex(0, i) == subindex(1, i)){
      c(i) = 1;
      d(i) = 0;
    }
  }
  return;
}
