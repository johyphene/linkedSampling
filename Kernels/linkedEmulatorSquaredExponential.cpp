#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cKSingle(double test, double train, double gammaVal){

	double cor = exp(-(pow(test - train, 2)/(pow(gammaVal,2))));
    return cor;
}

// [[Rcpp::export]]
double xiSingle(double mu, double sigma, double gammaVal, double wT_ik){

	double xi = (1/(sqrt(1 + (2*(pow(sigma,2)))/(pow(gammaVal,2))))) * exp(-(pow(mu - wT_ik, 2))/(2*pow(sigma,2) + pow(gammaVal,2))); 

	return xi;
}

// [[Rcpp::export]]
double psiSingle(double mu, double sigma, double gammaVal, double wT_jk){

  double psi = (1/(sqrt(1 + (2*(pow(sigma,2)))/(pow(gammaVal,2))))) * exp(-(pow(mu - wT_jk, 2))/(2*pow(sigma,2) + pow(gammaVal,2))) * (2*pow(sigma,2)*wT_jk + pow(gammaVal,2)*mu)/(2*pow(sigma,2) + pow(gammaVal,2)); 

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
double zetaSingle(double mu, double sigma, double gammaVal, double x1, double x2)
{
  double zetaVal = (1/(sqrt(1 + (4*(pow(sigma,2)))/(pow(gammaVal,2))))) * exp(-(pow(((x1 + x2)/2) - mu, 2))/(2*pow(sigma,2) + (pow(gammaVal,2)/2)) - (pow(x1 - x2,2))/(2*pow(gammaVal,2))); 

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