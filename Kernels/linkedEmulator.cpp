#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cKSingle(double test, double train, double gammaVal){

	double cor = (1 + (sqrt(5)*abs(test - train))/gammaVal + (5*pow(test - train, 2))/(3*pow(gammaVal,2)))*exp(-(sqrt(5)*abs(test - train))/gammaVal);

    return cor;
}

// [[Rcpp::export]]
double xiSingle(double mu, double sigma, double gammaVal, double wT_ik){

	double mu_A = mu - (sqrt(5)*pow(sigma,2))/gammaVal;
	double mu_B = mu + (sqrt(5)*pow(sigma,2))/gammaVal;
	NumericVector E_1 = {1 - (sqrt(5)*wT_ik)/(gammaVal) + (5*pow(wT_ik,2))/(3*pow(gammaVal,2)), sqrt(5)/gammaVal - (10*wT_ik)/(3*pow(gammaVal,2)), 5/(3*pow(gammaVal,2))};
	NumericVector E_2 = {1 + (sqrt(5)*wT_ik)/(gammaVal) + (5*pow(wT_ik,2))/(3*pow(gammaVal,2)), sqrt(5)/gammaVal + (10*wT_ik)/(3*pow(gammaVal,2)), 5/(3*pow(gammaVal,2))};
	NumericVector Lambda_11 = {1, mu_A, pow(mu_A,2) + pow(sigma,2)};
	NumericVector Lambda_12 = {0, 1, mu_A + wT_ik};
	NumericVector Lambda_21 = {1, -mu_B, pow(mu_B,2) + pow(sigma,2)};
	NumericVector Lambda_22 = {0, 1, -mu_B - wT_ik};
	double xi = exp((5*pow(sigma,2) + 2*sqrt(5)*gammaVal*(wT_ik - mu))/(2*pow(gammaVal,2)))*((sum(E_1*Lambda_11)*R::pnorm((mu_A - wT_ik)/(sigma), 0.0, 1.0, 1, 0)) + (sum(E_1*Lambda_12)*(sigma/sqrt(2*M_PI))*exp(-(pow(wT_ik - mu_A,2)/(2*pow(sigma,2)))))) + exp((5*pow(sigma,2) - 2*sqrt(5)*gammaVal*(wT_ik - mu))/(2*pow(gammaVal,2)))*((sum(E_2*Lambda_21)*R::pnorm((wT_ik - mu_B)/(sigma), 0.0, 1.0, 1, 0)) + (sum(E_2*Lambda_22)*(sigma/sqrt(2*M_PI))*exp(-(pow(wT_ik - mu_B,2)/(2*pow(sigma,2))))));

	return xi;
}

// [[Rcpp::export]]
double psiSingle(double mu, double sigma, double gammaVal, double wT_jk){

  double mu_A = mu - (sqrt(5)*pow(sigma,2))/gammaVal;
  double mu_B = mu + (sqrt(5)*pow(sigma,2))/gammaVal;
  NumericVector E_1 = {1 - (sqrt(5)*wT_jk)/gammaVal + (5*pow(wT_jk,2))/(3*pow(gammaVal,2)), sqrt(5)/gammaVal - (10*wT_jk)/(3*pow(gammaVal,2)), 5/(3*pow(gammaVal,2))};
  NumericVector E_2 = {1 + (sqrt(5)*wT_jk)/gammaVal + (5*pow(wT_jk,2))/(3*pow(gammaVal,2)), sqrt(5)/gammaVal + (10*wT_jk)/(3*pow(gammaVal,2)), 5/(3*pow(gammaVal,2))};
  NumericVector Lambda_61 = {mu_A, (pow(mu_A,2) + pow(sigma,2)), (pow(mu_A,3) + 3*mu_A*pow(sigma,2))};
  NumericVector Lambda_62 = {1, (mu_A + wT_jk), (pow(mu_A,2) + 2*pow(sigma,2) + pow(wT_jk,2) + mu_A*wT_jk)};
  NumericVector Lambda_71 = {-mu_B, (pow(mu_B,2) + pow(sigma,2)), (-pow(mu_B,3) - 3*mu_B*pow(sigma,2))};
  NumericVector Lambda_72 = {1, (-mu_B - wT_jk), (pow(mu_B,2) + 2*pow(sigma,2) + pow(wT_jk,2) + mu_B*wT_jk)};
  
  double psi = exp((5*pow(sigma,2) + 2*sqrt(5)*gammaVal*(wT_jk - mu))/(2*pow(gammaVal,2)))*((sum(E_1*Lambda_61)*R::pnorm((mu_A - wT_jk)/(sigma), 0.0, 1.0, 1, 0)) + (sum(E_1*Lambda_62)*(sigma/sqrt(2*M_PI))*exp(-(pow(wT_jk - mu_A,2)/(2*pow(sigma,2)))))) - exp((5*pow(sigma,2) - 2*sqrt(5)*gammaVal*(wT_jk - mu))/(2*pow(gammaVal,2)))*((sum(E_2*Lambda_71)*R::pnorm((wT_jk - mu_B)/(sigma), 0.0, 1.0, 1, 0)) + (sum(E_2*Lambda_72)*(sigma/sqrt(2*M_PI))*exp(-(pow(wT_jk - mu_B,2)/(2*pow(sigma,2))))));
  
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
  double mu_C = mu - ((2*sqrt(5)*pow(sigma,2))/gammaVal);
  double mu_D = mu + ((2*sqrt(5)*pow(sigma,2))/gammaVal);

  NumericVector Lambda_31 = {1, mu_C, pow(mu_C,2) + pow(sigma,2), pow(mu_C,3) + 3*pow(sigma,2)*mu_C, pow(mu_C,4) + 6*pow(sigma,2)*pow(mu_C,2) + 3*pow(sigma,4)};
  NumericVector Lambda_32 = {0, 1, mu_C + x2, pow(mu_C,2) + 2*pow(sigma,2) + pow(x2,2) + mu_C*x2, pow(mu_C,3) + pow(x2,3) + x2*pow(mu_C,2) + mu_C*pow(x2,2) + 3*pow(sigma,2)*x2 + 5*pow(sigma,2)*mu_C};
  NumericVector Lambda_41 = {1, mu, pow(mu,2) + pow(sigma,2), pow(mu,3) + 3*pow(sigma,2)*mu, pow(mu,4) + 6*pow(sigma,2)*pow(mu,2) + 3*pow(sigma,4)};
  NumericVector Lambda_42 = {0, 1, mu + x1, pow(mu,2) + 2*pow(sigma,2) + pow(x1,2) + mu*x1, pow(mu,3) + pow(x1,3) + x1*pow(mu,2) + mu*pow(x1,2) + 3*pow(sigma,2)*x1 + 5*pow(sigma,2)*mu};
  NumericVector Lambda_43 = {0, 1, mu + x2, pow(mu,2) + 2*pow(sigma,2) + pow(x2,2) + mu*x2, pow(mu,3) + pow(x2,3) + x2*pow(mu,2) + mu*pow(x2,2) + 3*pow(sigma,2)*x2 + 5*pow(sigma,2)*mu};
  NumericVector Lambda_51 = {1, -mu_D, pow(mu_D,2) + pow(sigma,2), -pow(mu_D,3) - 3*pow(sigma,2)*mu_D, pow(mu_D,4) + 6*pow(sigma,2)*pow(mu_D,2) + 3*pow(sigma,4)};
  NumericVector Lambda_52 = {0, 1, -mu_D - x1, pow(mu_D,2) + 2*pow(sigma,2) + pow(x1,2) + mu_D*x1, -pow(mu_D,3) - pow(x1,3) - x1*pow(mu_D,2) - mu_D*pow(x1,2) - 3*pow(sigma,2)*x1 - 5*pow(sigma,2)*mu_D};

  double E_30 = 1 + ((25*pow(x1,2)*pow(x2,2) - 3*sqrt(5)*(3*pow(gammaVal,3) + 5*gammaVal*x1*x2)*(x1 + x2) + (15*pow(gammaVal,2))*(pow(x1,2) + pow(x2,2) + 3*x1*x2))/(9*pow(gammaVal,4)));
  double E_31 = ((18*sqrt(5)*pow(gammaVal,3) + 15*sqrt(5)*gammaVal*(pow(x1,2) + pow(x2,2)) - (75*pow(gammaVal,2) + 50*x1*x2)*(x1 + x2) + 60*sqrt(5)*gammaVal*x1*x2)/(9*pow(gammaVal,4)));
  double E_32 = ((5*(5*pow(x1,2) + 5*pow(x2,2) + 15*pow(gammaVal,2) - 9*sqrt(5)*gammaVal*(x1 + x2) + 20*x1*x2))/(9*pow(gammaVal,4)));
  double E_33 = ((10*(3*sqrt(5)*gammaVal - 5*x1 - 5*x2))/(9*pow(gammaVal,4)));
  double E_34 = (25)/(9*pow(gammaVal,4));
  NumericVector E_3 = NumericVector::create(E_30, E_31, E_32, E_33, E_34);

  double E_40 = 1 + ((25*pow(x1,2)*pow(x2,2) + 3*sqrt(5)*(3*pow(gammaVal,3) - 5*gammaVal*x1*x2)*(x2 - x1) + (15*pow(gammaVal,2))*(pow(x1,2) + pow(x2,2) - 3*x1*x2))/(9*pow(gammaVal,4)));
  double E_41 = ((5*(3*sqrt(5)*gammaVal*(pow(x2,2) - pow(x1,2)) + (3*pow(gammaVal,2))*(x1 + x2) - 10*x1*x2*(x1 + x2)))/(9*pow(gammaVal,4)));
  double E_42 = ((5*(5*pow(x1,2) + 5*pow(x2,2) - 3*pow(gammaVal,2) - 3*sqrt(5)*gammaVal*(x2 - x1) + 20*x1*x2))/(9*pow(gammaVal,4)));
  double E_43 = -((50*(x1 + x2))/(9*pow(gammaVal,4)));
  double E_44 = (25)/(9*pow(gammaVal,4));
  NumericVector E_4 = NumericVector::create(E_40, E_41, E_42, E_43, E_44);

  double E_50 = 1 + ((25*pow(x1,2)*pow(x2,2) + 3*sqrt(5)*(3*pow(gammaVal,3) + 5*gammaVal*x1*x2)*(x1 + x2) + (15*pow(gammaVal,2))*(pow(x1,2) + pow(x2,2) + 3*x1*x2))/(9*pow(gammaVal,4)));
  double E_51 = ((18*sqrt(5)*pow(gammaVal,3) + 15*sqrt(5)*gammaVal*(pow(x1,2) + pow(x2,2)) + (75*pow(gammaVal,2) + 50*x1*x2)*(x1 + x2) + 60*sqrt(5)*gammaVal*x1*x2)/(9*pow(gammaVal,4)));
  double E_52 = ((5*(5*pow(x1,2) + 5*pow(x2,2) + 15*pow(gammaVal,2) + 9*sqrt(5)*gammaVal*(x1 + x2) + 20*x1*x2))/(9*pow(gammaVal,4)));
  double E_53 = ((10*(3*sqrt(5)*gammaVal + 5*x1 + 5*x2))/(9*pow(gammaVal,4)));
  double E_54 = (25)/(9*pow(gammaVal,4));
  NumericVector E_5 = NumericVector::create(E_50, E_51, E_52, E_53, E_54);

  double zetaV = exp((10*pow(sigma, 2) + sqrt(5)*gammaVal*(x1 + x2 - 2*mu))/(pow(gammaVal,2)))*(sum(E_3*Lambda_31)*R::pnorm((mu_C - x2)/(sigma), 0.0, 1.0, 1, 0) + sum(E_3*Lambda_32)*((sigma)/(sqrt(2*M_PI)))*exp(-pow((x2 - mu_C),2)/(2*pow(sigma,2)))) + exp(-((sqrt(5)*(x2 - x1))/(gammaVal)))*(sum(E_4*Lambda_41)*(R::pnorm((x2 - mu)/(sigma), 0.0, 1.0, 1, 0) - R::pnorm((x1 - mu)/(sigma), 0.0, 1.0, 1, 0)) + sum(E_4*Lambda_42)*((sigma)/(sqrt(2*M_PI)))*exp(-pow((x1 - mu),2)/(2*pow(sigma,2))) - sum(E_4*Lambda_43)*((sigma)/(sqrt(2*M_PI)))*exp(-pow((x2 - mu),2)/(2*pow(sigma,2)))) + exp((10*pow(sigma,2) - sqrt(5)*gammaVal*(x1 + x2 - 2*mu))/(pow(gammaVal,2)))*(sum(E_5*Lambda_51)*R::pnorm((x1 - mu_D)/(sigma), 0.0, 1.0, 1, 0) + sum(E_5*Lambda_52)*((sigma)/(sqrt(2*M_PI)))*exp(-pow((x1 - mu_D),2)/(2*pow(sigma,2))));

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

// [[Rcpp::export]]
NumericVector jFun2(NumericVector test, NumericMatrix trainW, NumericVector gamma, NumericVector mu, NumericVector sigma)
{
  int d = trainW.ncol();
  int m = trainW.nrow();
  NumericMatrix jMat(m,m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      jMat(i,j) = 1.0;
    }
  }
  for (int i = 0; i < m; i++) {
	for (int j = 0; j < m; j++) {
		for (int k = 0; k < d; k++) {
		  double zetaVal = zetaSingle(mu[k], sigma[k], gamma[k], trainW(i,k), trainW(j,k));
		  jMat(i,j) = jMat(i,j)*zetaVal;
		}
		jMat(j,i) = jMat(i,j); //copy across diagonal
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