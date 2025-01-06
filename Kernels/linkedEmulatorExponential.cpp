#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cKSingle(double test, double train, double gammaVal){

	double cor = exp(-(abs(test - train)/(gammaVal)));
    return cor;
}

// [[Rcpp::export]]
double xiSingle(double mu, double sigma, double gammaVal, double wT_ik){

	double mu_A = mu - pow(sigma,2)/gammaVal;
	double mu_B = mu + pow(sigma,2)/gammaVal;
	double xi = exp((pow(sigma,2) + 2*gammaVal * (wT_ik - mu))/(2*pow(gammaVal,2))) * R::pnorm((mu_A - wT_ik)/(sigma), 0.0, 1.0, 1, 0) + exp((pow(sigma,2) - 2*gammaVal * (wT_ik - mu))/(2*pow(gammaVal,2))) * R::pnorm((wT_ik - mu_B)/(sigma), 0.0, 1.0, 1, 0);

	return xi;
}

// [[Rcpp::export]]
double psiSingle(double mu, double sigma, double gammaVal, double wT_jk){

	double mu_A = mu - pow(sigma,2)/gammaVal;
	double mu_B = mu + pow(sigma,2)/gammaVal;
	double psi = exp((pow(sigma,2) + 2*gammaVal * (wT_jk - mu))/(2*pow(gammaVal,2))) * (mu_A * R::pnorm((mu_A - wT_jk)/(sigma), 0.0, 1.0, 1, 0) + (sigma/sqrt(2*M_PI)) * exp(-(pow(wT_jk - mu_A,2)/(2*pow(sigma,2))))) - exp((pow(sigma,2) - 2*gammaVal * (wT_jk - mu))/(2*pow(gammaVal,2))) * (mu_B * R::pnorm((wT_jk - mu_B)/(sigma), 0.0, 1.0, 1, 0) - (sigma/sqrt(2*M_PI)) * exp(-(pow(wT_jk - mu_B,2)/(2*pow(sigma,2)))));

	return psi;
}


// [[Rcpp::export]]
NumericVector iFun(NumericVector test, NumericMatrix trainW, NumericMatrix trainZ, NumericVector gamma, NumericVector mu, NumericVector sigma)
{
  int d = trainW.ncol();
  int p = trainZ.ncol();
  int m = trainW.nrow();
  NumericVector iVec(m);
  for (int i = 0; i < m; i++) {
    iVec[i] = 1.0;
  }
  for (int i = 0; i < m; i++) {
	for (int k = 0; k < p; k++) {
    double cVal = cKSingle(test[d+k], trainZ(i,k), gamma[d+k]);
	  iVec[i] = iVec[i]*cVal;
    }
    for (int k = 0; k < d; k++) {
    double xiVal = xiSingle(mu[k], sigma[k], gamma[k], trainW(i,k));
	  iVec[i] = iVec[i]*xiVal;
    }
  }
  
  return iVec;
}

// [[Rcpp::export]]
NumericMatrix bFun(NumericVector test, NumericMatrix trainW, NumericMatrix trainZ, NumericVector gamma, NumericVector mu, NumericVector sigma)
{
  int d = trainW.ncol();
  int p = trainZ.ncol();
  int m = trainW.nrow();
  NumericMatrix bMat(d,m);
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < m; j++) {
      bMat(i,j) = 1.0;
    }
  }
  for (int l = 0; l < d; l++) {
	  for (int j = 0; j < m; j++) {
      bMat(l,j) = bMat(l,j)*psiSingle(mu[l], sigma[l], gamma[l], trainW(j,l));
      for (int k = 0; k < d; k++) {
        if (k != l) {
          double xiVal = xiSingle(mu[k], sigma[k], gamma[k], trainW(j,k));
          bMat(l,j) = bMat(l,j)*xiVal;
        }
      }
	  for (int k = 0; k < p; k++) {
        double cVal = cKSingle(test[d+k], trainZ(j,k), gamma[d+k]);
        bMat(l,j) = bMat(l,j)*cVal;
      }
    }
  }

  return bMat;
}

// [[Rcpp::export]]
double hZeta(double mu, double sigma, double gammaVal, double x1, double x2) {

  double mu_C = mu - (2*pow(sigma,2))/gammaVal;
  double mu_D = mu + (2*pow(sigma,2))/gammaVal;
  double zetaV = exp((2*(pow(sigma,2)) + gammaVal * (x1 + x2 - 2*mu))/(pow(gammaVal,2))) * R::pnorm((mu_C - x2)/(sigma), 0.0, 1.0, 1, 0) + exp(-((x2-x1)/gammaVal))*(R::pnorm((x2-mu)/(sigma), 0.0, 1.0, 1, 0) - R::pnorm((x1-mu)/(sigma), 0.0, 1.0, 1, 0)) + exp((2*(pow(sigma,2)) - gammaVal * (x1 + x2 - 2*mu))/(pow(gammaVal,2))) * R::pnorm((x1 - mu_D)/(sigma), 0.0, 1.0, 1, 0);
  
  return zetaV;
}

// [[Rcpp::export]]
double zetaSingle(double mu, double sigma, double gammaVal, double x1, double x2)
{
  double zetaVal;
  if (x2 >= x1) {
    zetaVal = hZeta(mu, sigma, gammaVal, x1, x2);
  } else {
    zetaVal = hZeta(mu, sigma, gammaVal, x2, x1);
  }
  return zetaVal;
}

// [[Rcpp::export]]
NumericVector jFun(NumericVector test, NumericMatrix trainW, NumericMatrix trainZ, NumericVector gamma, NumericVector mu, NumericVector sigma)
{
  int d = trainW.ncol();
  int p = trainZ.ncol();
  int m = trainW.nrow();
  NumericMatrix jMat(m,m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      jMat(i,j) = 1.0;
    }
  }
  for (int i = 0; i < m; i++) {
	for (int j = 0; j < m; j++) {
		for (int k = 0; k < p; k++) {
		  double cVali = cKSingle(test[d+k], trainZ(i,k), gamma[d+k]);
		  double cValj = cKSingle(test[d+k], trainZ(j,k), gamma[d+k]);
		  jMat(i,j) = jMat(i,j)*cVali*cValj;
		}
		for (int k = 0; k < d; k++) {
		  double zetaVal = zetaSingle(mu[k], sigma[k], gamma[k], trainW(i,k), trainW(j,k));
		  jMat(i,j) = jMat(i,j)*zetaVal;
		}
		//jMat(j,i) = jMat(i,j); //copy across diagonal
	}
  }
  
  return jMat;
}

/*
// [[Rcpp::export]]
double gCFunction(double x1, double x2, double c) {
  double r = sqrt(pow(x1 - x2, 2)); // Euclidean distance calculation - might need to adjust for more dimensions
  
  double wght = 0;
  double absr = fabs(r);
  
  if (0 <= absr && absr <= c) {
    wght = -0.25*pow(absr/c,5) + 0.5*pow(absr/c,4) + 0.625*pow(absr/c,3) - (5.0/3.0) * pow(absr/c,2) + 1;
  }
  else if (c <= absr && absr <= 2 * c) {
    wght = (1.0/12.0)*pow(absr/c,5) - 0.5*pow(absr/c,4) + 0.625*pow(absr/c,3) + (5.0/3.0)*pow(absr/c,2) - 5*(absr/c) + 4 - (2.0/3.0)*(c/absr);
  }
  
  return wght;
}

// [[Rcpp::export]]
double gCFunction2D(NumericVector x1, NumericVector x2, double c) {
  double r = sqrt(pow(x1[0] - x2[0], 2) + pow(x1[1] - x2[1], 2)); // Euclidean distance calculation - might need to adjust for more dimensions
  
  double wght = 0;
  double absr = fabs(r);
  
  if (0 <= absr && absr <= c) {
    wght = -0.25*pow(absr/c,5) + 0.5*pow(absr/c,4) + 0.625*pow(absr/c,3) - (5.0/3.0) * pow(absr/c,2) + 1;
  }
  else if (c <= absr && absr <= 2 * c) {
    wght = (1.0/12.0)*pow(absr/c,5) - 0.5*pow(absr/c,4) + 0.625*pow(absr/c,3) + (5.0/3.0)*pow(absr/c,2) - 5*(absr/c) + 4 - (2.0/3.0)*(c/absr);
  }
  
  return wght;
}*/