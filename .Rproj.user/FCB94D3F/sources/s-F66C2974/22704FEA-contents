// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>

#include <math.h>
#include <iostream>
#include "PMEASc.h"


using namespace std;
using namespace arma;
using namespace Rcpp;

#ifndef DOUBLE_EPS
#define DOUBLE_EPS 1E-15;
#endif

// [[Rcpp::export]]
Rcpp::List PMEASc(arma::mat X, double stoptol, arma::vec Lambdapath,  std::string calmethod,  std::string stopmethod,  unsigned int maxiter,  bool printyes,  bool printyessub, double sigma, int numlam, bool maxlambdacheck){
  // set some given information
  stoptol = max(stoptol, 1e-6);
  int p = X.n_rows;
  int n = X.n_cols;
  mat A(p, n), S(p, p);
  vec eigenS(n);
  mat eigenU(p, n);

  if (calmethod != "eADMM") {
      A = findAM(X);
  }
  else {

      eigenS = findAandeigen(X, A, eigenU);
  }
  S = A * A.t();

  double primobj=0, dualobj=0, gap=0, primfeas=0, dualfeas=0, eta=0;
  int nnzOmega = 0;
  mat gradP, residual, abstriures;
  mat Ip = arma::eye(p, p);
  cout.precision(4);

  // deal the lambda path
  double maxlambda = findmaxlambda(S);
  if (printyes){
      cout << "\n Max lambda = " << fixed << maxlambda << endl;
      if (calmethod == "iADMM") {
          cout << " method = iADMM starts:" << endl;
      } else if (calmethod == "SSNAL") {
        cout << " method = SSNAL starts:" << endl;
      } else if (calmethod == "eADMM") {
        cout << " method = eADMM starts:" << endl;
      }
  }
  if (maxlambdacheck){
    uvec clamplambdaindex = find(Lambdapath < maxlambda);
    Lambdapath = Lambdapath(clamplambdaindex);
    if (!Lambdapath.is_sorted("descend")){
      Lambdapath = sort(Lambdapath,"descend");
    }
  }

  int pathlength = Lambdapath.n_elem;
  field<mat> Omegapath(pathlength,1);
  //vec nnzOmegapath(pathlength);
  vec timepath(pathlength);
  //vec nnsmallOmegapath(pathlength);

  // set initial values
  //mat Omega = mat(zeros(p, p));
  mat Omega = zeros(p, p);
  mat Y;
  if (calmethod != "eADMM") {
      Y = A;
      for (int i = 0; i < p; i++) {
          Omega(i, i) = 1 / S(i, i);
          Y.row(i) /= S(i, i);
      }
  }
  else {
      for (int i = 0; i < p; i++) {
          Omega(i, i) = 1 / S(i, i);
      }
      Y = Omega;
  }


  mat Z = zeros(p, p);

  wall_clock pathtime;
  pathtime.tic();
  if (printyes) {
      cout << "*********************************************************************" << endl;
      cout << " Iter|   lambda   |  primobj   dualobj     gap   | primfeas  dualfeas  eta   | time   nnzoff" << endl;
  }

  // main iteration
  for (int iterpath = 0; iterpath < pathlength; iterpath++)
  {
    wall_clock itertime;
    itertime.tic();
    double lambda = Lambdapath(iterpath);
    double timeiter;

    if (calmethod == "SSNAL") {
        PMESSNALmainc(A, lambda, stoptol, maxiter, Omega, Y, Z, p, n, sigma, printyessub, primobj, \
            dualobj, gap, primfeas, dualfeas, eta, nnzOmega, Ip);
    }
    else if (calmethod == "eADMM") {
        string admmmethod = "direct";
        PMEeADMMc(S, A, eigenS, eigenU, lambda, stoptol, maxiter, Omega, Y, Z, p, n, sigma, printyessub, primobj, \
            dualobj, gap, primfeas, dualfeas, eta, nnzOmega, Ip);
    }
    else if (calmethod == "iADMM") {
        string admmmethod = "inexact";
        PMEiADMMc(A, lambda, stoptol, maxiter, Omega, Y, Z, p, n, sigma, printyessub, primobj, \
            dualobj, gap, primfeas, dualfeas, eta, nnzOmega, Ip);
    }
    timeiter = itertime.toc();
    if (printyessub) {
        cout << " subproblem results printing end." << endl;
        cout << " ---------------------------------------" << endl;
    }
    Omegapath(iterpath) = sp_mat(Omega);
    //nnzOmegapath(iterpath) = nnzOmega;
    timepath(iterpath) = timeiter;
    //nnsmallOmegapath(iterpath) = accu(abs(Omega) > 1e-5) - p;
    if (printyes) {

      cout << "   " << dec << iterpath << " | "  << fixed << lambda << " | " << scientific << primobj << "  " << dualobj << "  " << gap \
             << " |  " << primfeas << "  " << dualfeas << "  " << eta << " | " << fixed << timeiter << "   " << dec << nnzOmega << endl;
    }

    if (stopmethod == "plot"){
      if (eta > 1e-5 || gap > 1e-5){
        if (pathlength > 1 && iterpath < pathlength - 1){
          if (iterpath == pathlength - 2){
            Omegapath = Omegapath.rows(0, iterpath);
            Lambdapath.shed_rows(pathlength - 1, pathlength - 1);
            //nnzOmegapath.shed_rows(pathlength - 1, pathlength - 1);
            //nnsmallOmegapath.shed_rows(pathlength - 1, pathlength - 1);
            timepath.shed_rows(pathlength - 1, pathlength - 1);
          }else{
            Omegapath = Omegapath.rows(0, iterpath);
            Lambdapath.shed_rows((iterpath + 1), (pathlength - 1));
            //nnzOmegapath.shed_rows((iterpath + 1), (pathlength - 1));
            //nnsmallOmegapath.shed_rows((iterpath + 1), (pathlength - 1));
            timepath.shed_rows((iterpath + 1), (pathlength - 1));
          }
        }
        break;
      }
    }

    if (stopmethod == "bigs"){
      if (nnzOmega > (numlam * n)){
        if (printyes){
          cout << " lambda small enough" << endl;
        }

        if (pathlength > 1 && iterpath < pathlength - 1){
          if (iterpath == pathlength - 2){
            Omegapath = Omegapath.rows(0, iterpath);
            Lambdapath.shed_rows(pathlength - 1, pathlength - 1);
            //nnzOmegapath.shed_rows(pathlength - 1, pathlength - 1);
            //nnsmallOmegapath.shed_rows(pathlength - 1, pathlength - 1);
            timepath.shed_rows(pathlength - 1, pathlength - 1);
          }else{
            Omegapath = Omegapath.rows(0, iterpath);
            Lambdapath.shed_rows((iterpath + 1), (pathlength - 1));
            //nnzOmegapath.shed_rows((iterpath + 1), (pathlength - 1));
            //nnsmallOmegapath.shed_rows((iterpath + 1), (pathlength - 1));
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
          //nnzOmegapath.shed_rows(pathlength - 1, pathlength - 1);
          //nnsmallOmegapath.shed_rows(pathlength - 1, pathlength - 1);
          timepath.shed_rows(pathlength - 1, pathlength - 1);
        }else{
          Omegapath = Omegapath.rows(0, iterpath);
          Lambdapath.shed_rows((iterpath + 1), (pathlength - 1));
          //nnzOmegapath.shed_rows((iterpath + 1), (pathlength - 1));
          //nnsmallOmegapath.shed_rows((iterpath + 1), (pathlength - 1));
          timepath.shed_rows((iterpath + 1), (pathlength - 1));
        }
      }


      break;
    }
  }// end of for loop

  if (printyes){
    double timepathwhole = pathtime.toc();
    cout << " total time for calculating solution path: " << fixed << timepathwhole << "(s)" << endl;
    cout << "*********************************************************************" << endl;
  }

  return Rcpp::List::create(Named("Omegapath") = Omegapath,
                            Named("Lambdapath") = Lambdapath,
                            Named("timepath") = timepath);
}// end of PMEc


// main function of SSNAL
void PMESSNALmainc(mat A, double lambda, double stoptol, int maxiter, mat &Omega, mat &Y, mat &Z, int p, int n, double sigma, bool printyessub, double &primobj,\
                   double &dualobj, double &gap, double &primfeas, double &dualfeas, double &eta, int &nnzOmega, mat Ip)
{
  // preparation
  int maxitersub = 10;
  double Stolconst = 0.5;
  int primwin = 0, dualwin = 0;
  bool breakyes = false;
  int iter = 0;


  mat YAT = Y * A.t();
  mat OmegaA = Omega * A;
  mat OmegaS = OmegaA * A.t();
  mat gradh = 0.5 * (OmegaS + OmegaS.t()) - Ip;


  primfeas = norm(Y - OmegaA, "fro") / (1 + norm(Y, "fro"));
  dualfeas = norm(0.5 * (YAT + YAT.t()) + Z - Ip, "fro") / (1 + sqrt(p));
  primobj = 0.5 * pow(norm(OmegaA, "fro"), 2) - accu(Omega.diag()) + lambda * oneoffnorm(Omega);
  dualobj = - 0.5 * pow(norm(Y, "fro"), 2);
  gap = (primobj - dualobj) / (1 + fabs(primobj) + fabs(dualobj));
  eta = norm(gradh + proxBmainM(Omega - gradh, lambda), "fro") / (1 + norm(gradh, "fro") + norm(Omega, "fro"));

  if ((std::max(primfeas, dualfeas) < 100 * std::max(1e-6, stoptol)) && (eta < stoptol)){
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
  while (iter < maxiter && breakyes == false) {
    iter++;
    if (dualfeas < 1e-3) {
      maxitersub = max(maxitersub, 30);
    } else if (dualfeas < 1e-1){
      maxitersub = max(maxitersub, 20);
    }

    int subbreakyes = 0;
    mat Ztmp = zeros(p, p);
    PMESSNCGc(Omega, Y, YAT, Z, Ztmp, A, Ip, subbreakyes, lambda, sigma, maxitersub, Stolconst,\
              stoptol, p, n);
    if (subbreakyes < 0){
      Stolconst = std::max(Stolconst / 1.06, 1e-3);
    }

    Omega = sigma * Ztmp;
    OmegaA = Omega * A;
    OmegaS = OmegaA * A.t();
    gradh = 0.5 * (OmegaS + OmegaS.t()) - Ip;

    primfeas = norm(Y - OmegaA, "fro") / (1 + norm(Y, "fro"));
    dualfeas = norm(0.5 * (YAT + YAT.t()) + Z - Ip, "fro") / (1 + sqrt(p));
    primobj = 0.5 * pow(norm(OmegaA, "fro"), 2) - accu(Omega.diag()) + lambda * oneoffnorm(Omega);
    dualobj = - 0.5 * pow(norm(Y, "fro"), 2);
    gap = (primobj - dualobj) / (1 + fabs(primobj) + fabs(dualobj));
    eta = norm(gradh + proxBmainM(Omega - gradh, lambda), "fro") / (1 + norm(gradh, "fro") + norm(Omega, "fro"));

    if (max(primfeas, dualfeas) < 500 * max(1e-6, stoptol) && eta < stoptol){
      breakyes = true;
    }

    if (primfeas < dualfeas){
      primwin = primwin + 1;
    }
    else{
      dualwin = dualwin + 1;
    }

    updatesigmaM(primwin, dualwin, sigma, iter, subbreakyes);
    if (printyessub){
      cout << "   " << dec << iter << " | " << scientific << primobj << "  " << dualobj << " "<< gap << " | " \
           << primfeas << " " << dualfeas << " " << eta << " " << sigma << endl;
    }


  }// end of while loop
  nnzOmega = accu(Omega != 0) - p;
  return;
}







void PMESSNCGc(mat Omega, mat &Y, mat &YAT, mat &Z, mat &Ztmp, mat A, mat Ip, int &subbreakyes, double lambda, double sigma, int maxitersub, double Stolconst,\
               double stoptol, int p, int n)
{
  int maxiterCG = 500;
  mat Zin = Omega / sigma - 0.5 * (YAT + YAT.t()) + Ip;
  Z = proxBmainM(Zin, lambda);
  Ztmp = Zin - Z;
  mat GradPsiY(p, n);
  double PsiY = 0.5 * pow(norm(Y, "fro"), 2) + 0.5 * sigma * pow(norm(Ztmp, "fro"), 2);
  vector<double> SSNCGsubprimfeas, SSNCGsubdualfeas, SSNCGPsiY;
  vector<int> SSNCGsolveok, SSNCGCGiter;
  mat res;

  for (int itersub = 0; itersub < maxitersub; itersub++) {
    GradPsiY = sigma * (Ztmp * A) - Y;
    double subprimfeas = norm(GradPsiY, "fro") / (1 + norm(Y, "fro"));
    double subdualfeas = norm(0.5 * (YAT + YAT.t())+ Z - Ip, "fro") / (1 + sqrt(p));
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
    SSNCGPsiY.push_back(PsiY);

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

    res = GradPsiY;
    double tolCG = min(5e-3, 0.1 * norm(res, "fro"));
    double tolCGconst = 1;
    if(itersub > 0 && (subprimfeas > 0.1 * SSNCGsubprimfeas[0] || primratio > 0.5)){
      tolCGconst = 0.5;
    }
    if(dualratio > 1.1){
      tolCGconst = 0.5;
    }
    tolCG = tolCGconst * tolCG;

    mat direction = zeros(p, n);
    vector<double> err;
    mat U = partgradient(Ztmp + Z, lambda);
    int solveok = 1;
    PMECG(res, tolCG, maxiterCG, A, U, p, n, sigma, direction, solveok, err);
    SSNCGsolveok.push_back(solveok);
    SSNCGCGiter.push_back(err.size()- 1);


    double steptol = 1e-5;
    double stepop;
    if (itersub < 2 || (itersub <= 2 && subdualfeas > 1e-4)){
      stepop = 1;
    }else
    {
      stepop = 2;
    }
    double alp = 1;

    findstep(GradPsiY, steptol, stepop, sigma, Ztmp, Z, Y, A, direction, PsiY, lambda, alp, p);

    YAT = Y * A.t();
    if ((alp) < (DOUBLE_EPS)){
      subbreakyes = 11;
    }

    if (itersub > 1000) {
      vec CGiter = vectoraccess(SSNCGCGiter, itersub - 3, itersub);
      vec tmp = vectoraccess(SSNCGsubprimfeas, itersub - 3, itersub);
      vec CGsolveok = vectoraccess(SSNCGsolveok, itersub - 3, itersub);
      double ratio = min(tmp) / max(tmp);
      if (all(CGsolveok <= -1) && ratio > 0.9 && min(CGiter) == max(CGiter) && max(tmp) < 5 * stoptol){
        subbreakyes = 1;
      }
      double const3 = 0.7;
      double priminf1half = min(vectoraccess(SSNCGsubprimfeas, 0, floor(itersub * const3)));
      double priminf2half = min(vectoraccess(SSNCGsubprimfeas, floor(itersub * const3), itersub));
      double priminfbest = min(vectoraccess(SSNCGsubprimfeas, 0, itersub));
      double priminfratio = SSNCGsubprimfeas[itersub] / SSNCGsubprimfeas[itersub - 1];
      uvec stagnateidx = find(vectoraccess(SSNCGsolveok, 0, itersub) <= -1);
      int stagnatecount = stagnateidx.n_elem;
      if (itersub >= 10 && all(vectoraccess(SSNCGsolveok, itersub - 3, itersub) == -1) && priminfbest < 1e-2 && subdualfeas < 1e-3) {
        tmp = vectoraccess(SSNCGsubprimfeas, std::max(0, itersub - 7), itersub);
        ratio = min(tmp) / max(tmp);
        if (ratio > 0.5) {
          subbreakyes = 2;
        }
      }
      if (itersub >= 15 && priminf1half < std::min(2e-3, priminf2half) && subdualfeas < 0.8 * SSNCGsubdualfeas[0] && subdualfeas < 1e-3 && stagnatecount >= 3) {
        subbreakyes = 3;
      }
      if (itersub >= 15 && priminfratio < 0.1 && subprimfeas < 0.8 * priminf1half && subdualfeas < min(1e-3, 2 * subprimfeas) && (subprimfeas < 1e-3 || (subdualfeas < 1e-5 && subprimfeas < 5e-3)) && stagnatecount >=3) {
        subbreakyes = 4;
      }
      if (itersub > 10 && subdualfeas > 5 * min(vectoraccess(SSNCGsubdualfeas, 0, itersub)) && subprimfeas > 2 * min(vectoraccess(SSNCGsubdualfeas, 0, itersub))) {
        subbreakyes = 5;
      }


      if (itersub > 20){
        vec dualinfratioall= vectoraccess(SSNCGsubdualfeas, 1, itersub) / vectoraccess(SSNCGsubdualfeas, 0, itersub - 1);
        uvec idx = find(dualinfratioall > 1);
        if (idx.n_elem > 3) {
          double dualinfincrement = mean(dualinfratioall(idx));
          if(dualinfincrement > 1.25) {
            subbreakyes = 6;
          }
        }
      }
      if (subbreakyes > 0) {
        break;
      }




    }


    err.clear();

  }


  return;
}




//CG
void PMECG(mat res, double tolCG, int maxiterCG, mat A, mat U, int p, int n, double sigma, mat &direction, int &solveok, vector<double> &err) {

  int stagnatecheck = 20;
  err.push_back(norm(res, "fro"));
  mat g = res, Vg;
  double rz1 = pow(err[0], 2);
  double rz2 = 1;

  double beta, alpha, denom, residual;
  //mat Ag;

  for (int iterCG = 0; iterCG < maxiterCG; iterCG++) {
    if (iterCG > 0){
      beta = rz1 / rz2;
      g = res + beta * g;
    }
    //Ag = g * A.t();
    Vg = g + 0.5 * sigma * OperatorVp(A, g, U) * A;
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
    residual = norm(res, "fro");
    err.push_back(residual);
    if (residual < tolCG) {break;}
    rz2 = rz1;
    rz1 = accu(res % res);
    if (iterCG > stagnatecheck){
      vec ratio(10);
      for (int i = 0; i < 10; i++)
      {
        ratio(i) = err[iterCG - i + 1] / err[iterCG - i];
      }
      if (min(ratio) > 0.997 && max(ratio) < 1.003){
        solveok = -1;
        break;
      }
    }


  }

  return;
}



// find step length

void findstep(mat GradPsiY, double steptol, double stepop, double sigma, mat &Ztmp, mat &Z, mat &Y, mat A, mat direction, double &PsiY, double lambda, double &alp, int p) {
  int maxiterstep = (int) (ceil(log(1 / (steptol + DOUBLE_EPS)) / log(2)));
  double c1 = 1e-4, c2 = 0.9;
  double change0 = accu(-GradPsiY % direction);
  // double change0 = accu((Y - sigma * Ztmp * A) % direction);
  if (change0 >= 0) {
    //cout << "warning" << endl;
    return;
  }
  double alpconst = 0.5, LB = 0, UB = 1, PsiYold = PsiY, change1;
  mat Yold = Y, gradstep1;
  mat Zold = Z;
  mat Ztmpold = Ztmp, Zin;
  double gLB = change0, gUB = 0;
  for (int iterstep = 0; iterstep < maxiterstep; iterstep++){
    if (iterstep > 0){
      alp = alpconst * (LB + UB);
    }
    Y = Yold + alp * direction;
    Zin = Ztmpold + Zold - 0.5 * alp * (direction * A.t() + A * direction.t());
    Z = proxBmainM(Zin, lambda);
    Ztmp = Zin - Z;
    PsiY = 0.5 * pow(norm(Y, "fro"), 2) + 0.5 * sigma * pow(norm(Ztmp, "fro"), 2);
    gradstep1 = Y - sigma * (Ztmp * A);
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
        /*if (iterstep == maxiterstep - 1){
          cout << "\n Warning maximum iter in finding step length reached" << endl;
        }*/
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







void PMEiADMMc(mat A, double lambda, double stoptol, int maxiter, mat &Omega, mat &Y, mat &Z, int p, int n, double sigma, bool printyessub, double &primobj,\
            double &dualobj, double &gap, double &primfeas, double &dualfeas, double &eta, int &nnzOmega, mat Ip)
{
    double gamma = 1.618;
    int maxiterCG = 200;
    double sigamupdateconst = 1.125;
    double sigmamax = 1e6;
    double sigmamin = 1e-4;
    mat YAT = Y * A.t();
    mat OmegaA = Omega * A;
    mat partgradD = 0.5 * (YAT + YAT.t());
    primfeas = norm(Y - Omega * A, "fro") / (1 + norm(Y, "fro"));
    dualfeas = norm(partgradD + Z - Ip, "fro") / (1 + sqrt(p));
    double maxfeas = std::max(primfeas, dualfeas);
    primobj = 0.5 * pow(norm(OmegaA, "fro"), 2) - accu(Omega.diag()) + lambda * oneoffnorm(Omega);
    dualobj = - 0.5 * pow(norm(Y, "fro"), 2);
    gap = (primobj - dualobj) / (1 + fabs(primobj) + fabs(dualobj));
    mat OmegaS = OmegaA * A.t();
    mat gradh = 0.5 * (OmegaS + OmegaS.t()) - Ip;
    eta = norm(gradh + proxBmainM(Omega - gradh, lambda), "fro") / (1 + norm(gradh, "fro") + norm(Omega, "fro"));

    bool breakyes = false;
    if (maxfeas < 500 * stoptol && eta < stoptol) {
        breakyes = true;
    }
    int primwin = 0, dualwin = 0;
    int iter = 0;
    if (printyessub) {
        cout << " ---------------------------------------" << endl;
        cout << " subproblem results printing start:" << endl;
        cout << " iter |    primobj  dualobj     gap   |   primfeas      dualfeas        eta     |   sigma   " << endl;
        cout << "   " << dec << iter << " | " << scientific << primobj << "  " << dualobj << " "<< gap << " | " \
         << primfeas << " " << dualfeas << " " << eta << " " << sigma << endl;
    }
    if (breakyes) {
        return;
    }

    for (iter = 1; iter <= maxiter; iter++) {

            double stoptolCG = std::max(0.9 * stoptol, std::min(1 / pow((double)iter, 1.1), 0.9 * maxfeas));
            mat res = OmegaA + sigma * (Ip - Z - partgradD) * A - Y;
            //mat res = sigma * (Omega / sigma + Ip - Z) * A - sigma * 0.5 * (YAT + A * Y.t()) * A - Y;
            if (dualfeas > 1e-3 || iter <= 5) {
                maxiterCG = std::max(maxiterCG, 200);
            }else if (dualfeas > 1e-4) {
                maxiterCG = std::max(maxiterCG, 300);
            }else if (dualfeas > 1e-5) {
                maxiterCG = std::max(maxiterCG, 400);
            }else if (dualfeas > 5e-6) {
                maxiterCG = std::max(maxiterCG, 500);
            }

            ADMMCG(Y, res, A, stoptolCG, maxiterCG, sigma);
            YAT = Y * A.t();
            partgradD = 0.5 * (YAT + YAT.t());
            mat Zin = Omega / sigma - partgradD + Ip;
            Z = proxBmainM(Zin, lambda);
            Omega = Omega - sigma * gamma * (partgradD + Z - Ip);
            OmegaA = Omega * A;


            primfeas = norm(Y - Omega * A, "fro") / (1 + norm(Y, "fro"));
            dualfeas = norm(partgradD + Z - Ip, "fro") / (1 + sqrt(p));
            maxfeas = std::max(primfeas, dualfeas);
            primobj = 0.5 * pow(norm(OmegaA, "fro"), 2) - accu(Omega.diag()) + lambda * oneoffnorm(Omega);
            dualobj = - 0.5 * pow(norm(Y, "fro"), 2);
            gap = (primobj - dualobj) / (1 + fabs(primobj) + fabs(dualobj));
            OmegaS = OmegaA * A.t();
            gradh = 0.5 * (OmegaS + OmegaS.t()) - Ip;
            eta = norm(gradh + proxBmainM(Omega - gradh, lambda), "fro") / (1 + norm(gradh, "fro") + norm(Omega, "fro"));
            if (maxfeas < 500 * stoptol && eta < stoptol) {
                breakyes = true;
            }



            if (printyessub){
                int printiter;
                if (iter <= 200) {
                    printiter = 20;
                } else if (iter <= 2000) {
                    printiter = 100;
                } else if (iter <= 10000) {
                    printiter = 200;
                } else {
                    printiter = 500;
                }
                if ((modM(iter, printiter) == 0 || iter == maxiter) || breakyes || iter < 15) {
                    cout << "   " << dec << iter << " | " << scientific << primobj << "  " << dualobj << " "<< gap << " | " \
                    << primfeas << " " << dualfeas << " " << eta << " " << sigma << endl;
                }
            }


            if (breakyes) {
                break;
            }
            double feasratio = primfeas / dualfeas;
            if (feasratio < 1) {
                primwin += 1;
            } else {
                dualwin += 1;
            }
            int sigmaupdateiter;
            if (iter < 30) {
                sigmaupdateiter = 3;
            }else if (iter < 60) {
                sigmaupdateiter = 6;
            }else if (iter < 120) {
                sigmaupdateiter = 12;
            }else if (iter < 250) {
                sigmaupdateiter = 25;
            }else if (iter < 500) {
                sigmaupdateiter = 50;
            }else {
                sigmaupdateiter = 100;
            }

            if (modM(iter, sigmaupdateiter) == 0) {
                if (iter < 2500) {
                    if ((double)primwin > std::max(1.0, 1.2 * dualwin)) {
                        primwin = 0;
                        sigma = std::min(sigmamax, sigma * sigamupdateconst);
                    }else if ((double)dualwin > std::max(1.0, 1.2 * primwin)) {
                        dualwin = 0;
                        sigma = std::max(sigmamin, sigma / sigamupdateconst);
                    }
                }
            }


    }//end of for loop
    nnzOmega = accu(Omega != 0) - p;
    return;
}




// dual admm with exact solve the subproblem
void PMEeADMMc(mat S, mat A, vec eigenS, mat eigenU, double lambda, double stoptol, int maxiter, mat& Omega, mat& Y, mat& Z, int p, int n, double sigma, bool printyessub, double& primobj, \
               double& dualobj, double& gap, double& primfeas, double& dualfeas, double& eta, int& nnzOmega, mat Ip) {
  double gamma = 1.618;
  double sigamupdateconst = 1.125;
  double sigmamax = 1e6;
  double sigmamin = 1e-4;

  mat YS = (Y * A) * A.t();
  mat partgradD = 0.5 * (YS + YS.t());
  primfeas = norm(Y - Omega, "fro") / (1 + norm(Y, "fro"));
  dualfeas = norm(partgradD + Z - Ip, "fro") / (1 + sqrt(p));
  double maxfeas = std::max(primfeas, dualfeas);
  primobj = 0.5 * pow(norm(Omega * A, "fro"), 2) - accu(Omega.diag()) + lambda * oneoffnorm(Omega);
  dualobj = -0.5 * pow(norm(Y * A, "fro"), 2);
  gap = (primobj - dualobj) / (1 + fabs(primobj) + fabs(dualobj));
  mat OmegaS = (Omega * A) * A.t();
  mat gradh = 0.5 * (OmegaS + OmegaS.t()) - Ip;
  eta = norm(gradh + proxBmainM(Omega - gradh, lambda), "fro") / (1 + norm(gradh, "fro") + norm(Omega, "fro"));

  bool breakyes = false;
  if (maxfeas < 500 * stoptol && eta < stoptol) {
    breakyes = true;
  }
  int primwin = 0, dualwin = 0;
  int iter = 0;
  if (printyessub) {
    cout << " ---------------------------------------" << endl;
    cout << " subproblem results printing start:" << endl;
    cout << " iter |    primobj  dualobj     gap   |   primfeas      dualfeas        eta     |   sigma   " << endl;
    cout << "   " << dec << iter << " | " << scientific << primobj << "  " << dualobj << " " << gap << " | " \
         << primfeas << " " << dualfeas << " " << eta << " " << sigma << endl;
  }
  if (breakyes) {
    return;
  }
  bool sigmachange = true;
  mat Lam1(n, n), Lam2(n, n);
  for (iter = 1; iter <= maxiter; iter++) {
    if (sigmachange) {
      double rau = 2 / sigma;
      vec eigenl(n);
      for (int i = 0; i < n; i++) {
        double eigeni = eigenS(i);
        eigenl[i] = eigeni / (eigeni + rau);
        for (int j = 0; j < n; j++) {
          double eigenj = eigenS(j);
          double deno = (eigeni + rau) * (eigenj + rau) * (eigeni + eigenj + rau);
          Lam2(i, j) = eigeni * eigenj * (eigeni + eigenj + 2 * rau) / deno;
        }
      }
      Lam1 = diagmat(eigenl);
      sigmachange = false;
    }

    mat C = Omega / sigma + Ip - Z;
    mat CP1 = (C * eigenU) * (Lam1 * eigenU.t());
    mat P2 = eigenU * (Lam2 % (eigenU.t() * C * eigenU)) * eigenU.t();
    Y = sigma * (C - CP1 - CP1.t() + P2);
    Y = 0.5 * (Y + Y.t());
    YS = (Y * A) * A.t();
    partgradD = 0.5 * (YS + YS.t());
    mat Zin = Omega / sigma - partgradD + Ip;
    Z = proxBmainM(Zin, lambda);
    Omega = Omega - sigma * gamma * (partgradD + Z - Ip);


    primfeas = norm(Y - Omega, "fro") / (1 + norm(Y, "fro"));
    dualfeas = norm(partgradD + Z - Ip, "fro") / (1 + sqrt(p));
    maxfeas = std::max(primfeas, dualfeas);
    primobj = 0.5 * pow(norm(Omega * A, "fro"), 2) - accu(Omega.diag()) + lambda * oneoffnorm(Omega);
    dualobj = -0.5 * pow(norm(Y * A, "fro"), 2);
    gap = (primobj - dualobj) / (1 + fabs(primobj) + fabs(dualobj));
    OmegaS = Omega * S;
    gradh = 0.5 * (OmegaS + OmegaS.t()) - Ip;
    eta = norm(gradh + proxBmainM(Omega - gradh, lambda), "fro") / (1 + norm(gradh, "fro") + norm(Omega, "fro"));
    if (maxfeas < 500 * stoptol && eta < stoptol) {
      breakyes = true;
    }

    if (printyessub) {
      int printiter;
      if (iter <= 200) {
        printiter = 20;
      }
      else if (iter <= 2000) {
        printiter = 100;
      }
      else if (iter <= 10000) {
        printiter = 200;
      }
      else {
        printiter = 500;
      }
      if ((modM(iter, printiter) == 0 || iter == maxiter) || breakyes || iter < 20) {
        cout << "   " << dec << iter << " | " << scientific << primobj << "  " << dualobj << " " << gap << " | " \
             << primfeas << " " << dualfeas << " " << eta << " " << sigma << endl;
      }
    }


    if (breakyes) {
      break;
    }
    double feasratio = primfeas / dualfeas;
    if (feasratio < 1) {
      primwin += 1;
    }
    else {
      dualwin += 1;
    }
    int sigmaupdateiter;
    if (iter < 30) {
      sigmaupdateiter = 3;
    }
    else if (iter < 60) {
      sigmaupdateiter = 6;
    }
    else if (iter < 120) {
      sigmaupdateiter = 12;
    }
    else if (iter < 250) {
      sigmaupdateiter = 25;
    }
    else if (iter < 500) {
      sigmaupdateiter = 50;
    }
    else {
      sigmaupdateiter = 100;
    }

    if (modM(iter, sigmaupdateiter) == 0) {
      if (iter < 2500) {
        if ((double)primwin > std::max(1.0, 1.2 * dualwin)) {
          primwin = 0;
          sigma = std::min(sigmamax, sigma * sigamupdateconst);
          sigmachange = true;
        }
        else if ((double)dualwin > std::max(1.0, 1.2 * primwin)) {
          dualwin = 0;
          sigma = std::max(sigmamin, sigma / sigamupdateconst);
          sigmachange = true;
        }
      }
    }



  }// end of for loop

  nnzOmega = accu(mat(Omega) != 0) - p;
  return;
}






void ADMMCG(mat &Y, mat res, mat A, double tolCG, int maxiterCG, double sigma) {
  //int p = Y.n_rows;
  //int n = Y.n_cols;
  int stagnatecheck = 20;
  std::vector<double> err;
  err.push_back(norm(res, "fro"));
  mat g = res, Vg;
  double rz1 = pow(err[0], 2);
  double rz2 = 1;
  //cout << "In CG" << endl;
  double beta, alpha, denom, residual;
  mat Agt;

  for (int iterCG = 0; iterCG < maxiterCG; iterCG++) {
    if (iterCG > 0){
      beta = rz1 / rz2;
      g = res + beta * g;
    }
    Agt = g * A.t();
    Vg = g + 0.5 * sigma * (Agt.t() + Agt) * A;
    denom = accu(g % Vg);
    if ((fabs(denom)) < (DOUBLE_EPS)){
      //solveok = 2;
      break;
    } else
    {
      alpha = rz1 / denom;
      Y = Y + alpha * g;
      res = res - alpha * Vg;
    }
    residual = norm(res, "fro");
    // cout << residual << endl;
    err.push_back(residual);
    if (residual < tolCG) {break;}
    rz2 = rz1;
    rz1 = accu(res % res);
    if (iterCG > stagnatecheck) {
      vec ratio(10);
      for (int i = 0; i < 10; i++)
      {
        ratio(i) = err[static_cast<std::vector<double, std::allocator<double>>::size_type>(iterCG) - i + 1] / err[static_cast<std::vector<double, std::allocator<double>>::size_type>(iterCG) - i];
      }
      if (min(ratio) > 0.997 && max(ratio) < 1.003) {
        //solveok = -1;
        break;
      }
    }
  }

  return;

}
















int modM(int a, int n)
{
  return a - (int)(floor(a / n) * n);
}




















































/* ****************************************************************************************************** */

/* From given data X to find the matrix A */

mat findAM(mat X){
  int p = X.n_rows;
  int n = X.n_cols;
  double rowmeans;
  mat centralizedX = X, U, V;
  mat A = zeros(p, n);
  vec s;

  // centralize
  for(int i = 0; i < p; i++){
    rowmeans = mean(X.row(i));
    for(int j = 0; j < n; j++){
      centralizedX(i,j) = centralizedX(i, j) - rowmeans;
    }
  }
  svd_econ(U, s, V, centralizedX,"left");
  double sqrtn1 = sqrt(n - 1);
  for (int j = 0; j < n; j++){
    A.col(j) = U.col(j) * s(j) / sqrtn1;
  }
  return A;
}

/* From given data X to find the matrix A */

vec findAandeigen(mat X, mat &A, mat &U) {
  int p = X.n_rows;
  int n = X.n_cols;
  double rowmeans;
  mat centralizedX = X, V;
  //mat A = zeros(p, n);
  vec s;

  // centralize
  for (int i = 0; i < p; i++) {
    rowmeans = mean(X.row(i));
    for (int j = 0; j < n; j++) {
      centralizedX(i, j) = centralizedX(i, j) - rowmeans;
    }
  }
  svd_econ(U, s, V, centralizedX, "left");
  double sqrtn1 = sqrt(n - 1);
  for (int j = 0; j < n; j++) {
    A.col(j) = U.col(j) * s(j) / sqrtn1;
  }
  s = (s % s) / (n - 1);
  return s;
}


// calculate the diagonal-off one norm of a symmetric matrix

double oneoffnorm(mat X) {
  double oneoffabs = 0;
  int p = X.n_rows;
  for (int i = 0; i < (p - 1); i++)
  {
    for(int j = i + 1; j < p; j++)
    {
      oneoffabs += fabs(X(i, j));
    }
  }
  oneoffabs *= 2;
  return oneoffabs;
}


// original indicator function delta_{B_\lambda}

mat proxBmainM(mat XX, double lambda) {
  int p = XX.n_rows;
  mat  XXP = zeros(p, p);
  double elem;

  for (int i = 0; i < (p - 1); i++){
    for (int j = i + 1; j < p; j++){
      elem = max(-lambda, min(lambda, XX(i,j)));
      XXP(i,j) = elem;
      XXP(j,i) = elem;
    }
  }
  return mat(XXP);
}



// update sigma
void updatesigmaM(int &primwin, int &dualwin, double &sigma, int iter, int subbreakyes){
  double minsig = 1e-4;
  double maxsig = 1e7;
  int sigiter;
  if (iter < 10){
    sigiter = 3;
  }else if (iter < 200)
  {
    sigiter = 5;
  }else if (iter < 500){
    sigiter = 10;
  }else{
    sigiter = 20;
  }
  //double dual = (double)dualwin, prim = (double)primwin;
  if (modM(iter, sigiter) == 0 && subbreakyes < 0){
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
    sigma = max(minsig, 2 * sigma / 5);
  }
  return;
}


// calculate hessian

mat partgradient(mat Zin, double lambda) {
  int p = Zin.n_rows;
  mat U = ones(p, p);
  for (int i = 0; i < (p - 1); i++) {
    for (int j = i + 1; j < p; j++) {
      //cout << i << " , " << j << endl;
      if (fabs(Zin(i, j)) <= lambda){
        U(i, j) = 0;
        U(j, i) = 0;
      }
    }
  }
  return U;
}


mat OperatorVp(mat A, mat g, mat U) {

  int k,j;
  mat AgU = diagmat(U);
  sp_mat supperU = sp_mat(trimatu(U));

  sp_mat::const_iterator i = supperU.begin();
  sp_mat::const_iterator i_end = supperU.end();

  for (; i != i_end; ++i)
  {
    k = i.row();
    j = i.col();
    // Ag = g * A.t()
    //AgU(k,j) = Ag(k,j) + Ag(j,k);
    AgU(k,j) = dot(g.row(k),A.row(j)) + dot(g.row(j),A.row(k));
    AgU(j,k) = AgU(k,j);
  }
  return AgU;
}


// find the maximum lambda

double findmaxlambda(mat S){
  int p = S.n_rows;
  double elementS, maxlambda = 0;
  for(int i = 0; i < p; i++){
    for(int j = i + 1; j < p; j++){
      elementS = 0.5 * abs(S(i,j)/S(i,i) + S(i,j)/S(j,j));
      if (elementS > maxlambda){
        maxlambda = elementS;
      }
    }
  }
  return maxlambda;
}



template <class T>
arma::vec vectoraccess(std::vector<T> a, int b, int c){
  int lengthvec = c - b + 1;
  vec partvec(lengthvec);
  for(int i = 0; i < lengthvec; i++){
    partvec(i) = a[i];
  }
  return partvec;
}
