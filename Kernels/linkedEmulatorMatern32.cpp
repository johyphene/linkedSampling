#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cKSingle(double test, double train, double gammaVal){
	
	double cor = (1 + (sqrt(3)*abs(test - train))/gammaVal)*exp(-(sqrt(3)*abs(test - train))/gammaVal);
	
	return cor;
}

// [[Rcpp::export]]
double xiSingle(double mu, double sigma, double gammaVal, double wT_ik){

	double mu_A = mu - (sqrt(3)*pow(sigma,2))/gammaVal;
	double mu_B = mu + (sqrt(3)*pow(sigma,2))/gammaVal;
	NumericVector E_1 = {1 - (sqrt(3)*wT_ik)/(gammaVal), sqrt(3)/gammaVal};
	NumericVector E_2 = {1 + (sqrt(3)*wT_ik)/(gammaVal), sqrt(3)/gammaVal};
	NumericVector Lambda_11 = {1, mu_A};
	NumericVector Lambda_12 = {0, 1};
	NumericVector Lambda_21 = {1, -mu_B};
	NumericVector Lambda_22 = {0, 1};
	
	double xi = exp((3*pow(sigma,2) + 2*sqrt(3)*gammaVal*(wT_ik - mu))/(2*pow(gammaVal,2)))*((sum(E_1*Lambda_11)*R::pnorm((mu_A - wT_ik)/(sigma), 0.0, 1.0, 1, 0)) + (sum(E_1*Lambda_12)*(sigma/sqrt(2*M_PI))*exp(-(pow(wT_ik - mu_A,2)/(2*pow(sigma,2)))))) + exp((3*pow(sigma,2) - 2*sqrt(3)*gammaVal*(wT_ik - mu))/(2*pow(gammaVal,2)))*((sum(E_2*Lambda_21)*R::pnorm((wT_ik - mu_B)/(sigma), 0.0, 1.0, 1, 0)) + (sum(E_2*Lambda_22)*(sigma/sqrt(2*M_PI))*exp(-(pow(wT_ik - mu_B,2)/(2*pow(sigma,2))))));

	return xi;
}

// [[Rcpp::export]]
double psiSingle(double mu, double sigma, double gammaVal, double wT_jk){

  double mu_A = mu - (sqrt(3)*pow(sigma,2))/gammaVal;
  double mu_B = mu + (sqrt(3)*pow(sigma,2))/gammaVal;
  NumericVector E_1 = {1 - (sqrt(3)*wT_jk)/gammaVal, sqrt(3)/gammaVal};
  NumericVector E_2 = {1 + (sqrt(3)*wT_jk)/gammaVal, sqrt(3)/gammaVal};
  NumericVector Lambda_61 = {mu_A, (pow(mu_A,2) + pow(sigma,2))};
  NumericVector Lambda_62 = {1, (mu_A + wT_jk)};
  NumericVector Lambda_71 = {-mu_B, (pow(mu_B,2) + pow(sigma,2))};
  NumericVector Lambda_72 = {1, (-mu_B - wT_jk)};
    
  double psi = exp((3*pow(sigma,2) + 2*sqrt(3)*gammaVal*(wT_jk - mu))/(2*pow(gammaVal,2)))*((sum(E_1*Lambda_61)*R::pnorm((mu_A - wT_jk)/(sigma), 0.0, 1.0, 1, 0)) + (sum(E_1*Lambda_62)*(sigma/sqrt(2*M_PI))*exp(-(pow(wT_jk - mu_A,2)/(2*pow(sigma,2)))))) - exp((3*pow(sigma,2) - 2*sqrt(3)*gammaVal*(wT_jk - mu))/(2*pow(gammaVal,2)))*((sum(E_2*Lambda_71)*R::pnorm((wT_jk - mu_B)/(sigma), 0.0, 1.0, 1, 0)) + (sum(E_2*Lambda_72)*(sigma/sqrt(2*M_PI))*exp(-(pow(wT_jk - mu_B,2)/(2*pow(sigma,2))))));
  
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
  
  double mu_C = mu - ((2*sqrt(3)*pow(sigma,2))/gammaVal);
  double mu_D = mu + ((2*sqrt(3)*pow(sigma,2))/gammaVal);

  NumericVector Lambda_31 = {1, mu_C, pow(mu_C,2) + pow(sigma,2)};
  NumericVector Lambda_32 = {0, 1, mu_C + x2};
  NumericVector Lambda_41 = {1, mu, pow(mu,2) + pow(sigma,2)};
  NumericVector Lambda_42 = {0, 1, mu + x1};
  NumericVector Lambda_43 = {0, 1, mu + x2};
  NumericVector Lambda_51 = {1, -mu_D, pow(mu_D,2) + pow(sigma,2)};
  NumericVector Lambda_52 = {0, 1, -mu_D - x1};
  
  NumericVector E_3 = {1 + ((3*x1*x2 - sqrt(3)*gammaVal*(x1+x2))/(pow(gammaVal,2))), (2*sqrt(3)*gammaVal-3*(x1 + x2))/(pow(gammaVal,2)), 3/(pow(gammaVal,2))};
  NumericVector E_4 = {1 + ((sqrt(3)*gammaVal*(x2 - x1) - 3*x1*x2)/(pow(gammaVal,2))), (3*(x1 + x2))/(pow(gammaVal,2)), -3/(pow(gammaVal,2))};
  NumericVector E_5 = {1 + ((3*x1*x2 + sqrt(3)*gammaVal*(x1+x2))/(pow(gammaVal,2))), (2*sqrt(3)*gammaVal+3*(x1 + x2))/(pow(gammaVal,2)), 3/(pow(gammaVal,2))};

  double zetaV = exp((6*pow(sigma, 2) + sqrt(3)*gammaVal*(x1 + x2 - 2*mu))/(pow(gammaVal,2)))*(sum(E_3*Lambda_31)*R::pnorm((mu_C - x2)/(sigma), 0.0, 1.0, 1, 0) + sum(E_3*Lambda_32)*((sigma)/(sqrt(2*M_PI)))*exp(-pow((x2 - mu_C),2)/(2*pow(sigma,2)))) + exp(-((sqrt(3)*(x2 - x1))/(gammaVal)))*(sum(E_4*Lambda_41)*(R::pnorm((x2 - mu)/(sigma), 0.0, 1.0, 1, 0) - R::pnorm((x1 - mu)/(sigma), 0.0, 1.0, 1, 0)) + sum(E_4*Lambda_42)*((sigma)/(sqrt(2*M_PI)))*exp(-pow((x1 - mu),2)/(2*pow(sigma,2))) - sum(E_4*Lambda_43)*((sigma)/(sqrt(2*M_PI)))*exp(-pow((x2 - mu),2)/(2*pow(sigma,2)))) + exp((6*pow(sigma,2) - sqrt(3)*gammaVal*(x1 + x2 - 2*mu))/(pow(gammaVal,2)))*(sum(E_5*Lambda_51)*R::pnorm((x1 - mu_D)/(sigma), 0.0, 1.0, 1, 0) + sum(E_5*Lambda_52)*((sigma)/(sqrt(2*M_PI)))*exp(-pow((x1 - mu_D),2)/(2*pow(sigma,2))));

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