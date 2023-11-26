##BURGERS LINKED##
burgers <- function(xl, xr, tb, te, M, N, alf, bet, D)
{
  l = alf
  r = 0*t
  h = (xr-xl)/M
  k=(te-tb)/N
  m = M+1
  n = N
  sigma=D*k/(h*h)
  index = 1
  
  # for ii = xl:h:xr
  vals <- seq(from = xl, to = xr, by = h)
  for (ii in vals) {
    if (ii <= 0.25) {
      temp[index] = alf
    }
    if (ii > 0.25 && ii < 0.75) {
      temp[index] = -2*alf*ii + 1.5*alf
    }
    if (ii >= 0.75) {
      temp[index] = 0
    }
    index = index + 1
  }
  w[,1] = t(temp)
  w1 = w
  # for j=1:n
  for (j in 1:n) {
    # for it=1:3
    for (it in 1:3) {
      DF1 = zeros(m)
      DF2 = zeros(m)
      DF1 = diag(1+2*sigma, m)+rbind(zeros(1,m), cbind(diag(-sigma, m-1), zeros((m-1),1)))
      DF1 = DF1+rbind(cbind(zeros((m-1),1), diag(-sigma, m-1)), zeros(1,m))
      # DF2 = diag([0;k*w1(3:m)/(2*h);0])-diag([0;k*w1(1:m-2)/(2*h);0]);
      DF2 = diag(k*w1[3:m]*(1/(2*h)))
      DF2 = rbind(zeros(1, (m-2)), DF2, zeros(1, (m-2)))
      DF2 = cbind(zeros(m,1), DF2, zeros(m,1))
      DF2 = DF2 - cbind(zeros(m,1), rbind(zeros(1, (m-2)), diag(k*w1[1:(m-2)]*(1/(2*h))), zeros(1, (m-2))), zeros(m,1))
      # DF2 = DF2+diag([0;k*w1(2:m-1)/(2*h)],1)-diag([k*w1(2:m-1)/(2*h);0],-1);
      DF22 = diag(k*w1[2:m-1]*(1/(2*h)))
      DF22 = cbind(zeros((m-2),2), DF22)
      DF22 = rbind(zeros(1,m), DF22, zeros(1,m))
      DF23 = diag(k*w1[2:m-1]*(1/(2*h)))
      DF23 = cbind(DF23, zeros((m-2),2))
      DF23 = rbind(zeros(1,m), DF23, zeros(1,m))
      DF2 = DF2 + DF22 - DF23
      DF = DF1+DF2
      FMat = -w[,j]+(DF1+DF2*0.5)%*%w1
      DF[1,] = rbind(1, zeros(1, (m-1))) #[1 zeros(1,m-1)]; 
      DF[m,] = rbind(zeros(1, (m-1)),1) #[zeros(1, m-1) 1];
      FMat[1] = w1[1]-l[j] 
      FMat[m] = w1[m]-r[j] #look here if there are issues
      w1 = w1-solve(DF, FMat)
    }
    w[,j+1]=w1
  }
}

mu2 <- function(b1, b2, b3, w, z)
{
  return (b3+as.vector(b1)*w + as.vector(b2)*z)
}

#Matern 5/2
cKSingle <- function(test, train, gammaVal)
{
  cor <- (1 + (sqrt(5)*abs(test - train))/gammaVal + (5*(test - train)^2)/(3*gammaVal^2))*exp(-(sqrt(5)*abs(test - train))/gammaVal)
  return (cor)
}

cFunction <- function(X)
{
  cMat <- matrix(1, nrow = length(X[,1,1]), ncol = length(X[1,,1]))
  for (i in 1:length(X[,1,1])) {
    for (j in 1:length(X[1,,1])) {
      for (k in 1:length(X[1,1,])) {
        cMat[i, j] <- cMat[i,j]*X[i,j,k]
      }
    }
  }
  return (cMat)
}

cKSingle <- function(test, train, gammaVal)
{
  cor <- (1 + (sqrt(5)*abs(test - train))/gammaVal + (5*(test - train)^2)/(3*gammaVal^2))*exp(-(sqrt(5)*abs(test - train))/gammaVal)
  return (cor)
}

xiSingle <- function(mu, sigma, gammaVal, wT_ik)
{
  mu_A <- mu - (sqrt(5)*sigma^2)/gammaVal
  mu_B <- mu + (sqrt(5)*sigma^2)/gammaVal
  E_1 <- rbind(1 - (sqrt(5)*wT_ik)/(gammaVal) + (5*wT_ik^2)/(3*gammaVal^2), sqrt(5)/gammaVal - (10*wT_ik)/(3*gammaVal^2), 5/(3*gammaVal^2))
  E_2 <- rbind(1 + (sqrt(5)*wT_ik)/(gammaVal) + (5*wT_ik^2)/(3*gammaVal^2), sqrt(5)/gammaVal + (10*wT_ik)/(3*gammaVal^2), 5/(3*gammaVal^2))
  Lambda_11 <- rbind(1, mu_A, mu_A^2 + sigma^2)
  Lambda_12 <- rbind(0, 1, mu_A + wT_ik)
  Lambda_21 <- rbind(1, -mu_B, mu_B^2 + sigma^2)
  Lambda_22 <- rbind(0, 1, -mu_B - wT_ik)
  xi <- exp((5*sigma^2 + 2*sqrt(5)*gammaVal*(wT_ik - mu))/(2*gammaVal^2)) * ((t(E_1)%*%Lambda_11*pnorm((mu_A - wT_ik)/(sigma))) + (t(E_1)%*%Lambda_12*(sigma/sqrt(2*pi))*exp(-((wT_ik - mu_A)^2/(2*sigma^2))))) + exp((5*sigma^2 - 2*sqrt(5)*gammaVal*(wT_ik - mu))/(2*gammaVal^2)) * ((t(E_2)%*%Lambda_21*pnorm((wT_ik - mu_B)/(sigma))) + (t(E_2)%*%Lambda_22*(sigma/sqrt(2*pi))*exp(-((wT_ik - mu_B)^2/(2*sigma^2)))))
  
  # if ((wT_ik == mu) & (sigma == 0))
  # {
  #   xi <- exp((5*sigma^2 + 2*sqrt(5)*gammaVal*(wT_ik - mu))/(2*gammaVal^2)) * ((t(E_1)%*%Lambda_11*pnorm(0)) + (t(E_1)%*%Lambda_12*(sigma/sqrt(2*pi))*exp(0))) + exp((5*sigma^2 - 2*sqrt(5)*gammaVal*(wT_ik - mu))/(2*gammaVal^2)) * ((t(E_2)%*%Lambda_21*pnorm(0)) + (t(E_2)%*%Lambda_22*(sigma/sqrt(2*pi))*exp(0)))
  # }
  return (xi)
}

psiSingle <- function(mu, sigma, gammaVal, wT_jk)
{
  mu_A <- mu - (sqrt(5)*sigma^2)/gammaVal
  mu_B <- mu + (sqrt(5)*sigma^2)/gammaVal
  E_1 <- rbind(1 - (sqrt(5)*wT_jk)/gammaVal + (5*wT_jk^2)/(3*gammaVal^2), sqrt(5)/gammaVal - (10*wT_jk)/(3*gammaVal^2), 5/(3*gammaVal^2))
  E_2 <- rbind(1 + (sqrt(5)*wT_jk)/gammaVal + (5*wT_jk^2)/(3*gammaVal^2), sqrt(5)/gammaVal + (10*wT_jk)/(3*gammaVal^2), 5/(3*gammaVal^2))
  Lambda_61 <- rbind(mu_A, (mu_A^2 + sigma^2), (mu_A^3 + 3*mu_A*sigma^2))
  Lambda_62 <- rbind(1, (mu_A + wT_jk), (mu_A^2 + 2*sigma^2 + wT_jk^2 + mu_A*wT_jk))
  Lambda_71 <- rbind(-mu_B, (mu_B^2 + sigma^2), (-mu_B^3 - 3*mu_B*sigma^2))
  Lambda_72 <- rbind(1, (-mu_B - wT_jk), (mu_B^2 + 2*sigma^2 + wT_jk^2 + mu_B*wT_jk))
  
  psi <- exp((5*sigma^2 + 2*sqrt(5)*gammaVal*(wT_jk - mu))/(2*gammaVal^2)) * ((t(E_1)%*%Lambda_61*pnorm((mu_A - wT_jk)/(sigma))) + (t(E_1)%*%Lambda_62*(sigma/sqrt(2*pi))*exp(-((wT_jk - mu_A)^2/(2*sigma^2))))) - exp((5*sigma^2 - 2*sqrt(5)*gammaVal*(wT_jk - mu))/(2*gammaVal^2)) * ((t(E_2)%*%Lambda_71*pnorm((wT_jk - mu_B)/(sigma))) + (t(E_2)%*%Lambda_72*(sigma/sqrt(2*pi))*exp(-((wT_jk - mu_B)^2/(2*sigma^2)))))
  
  # if ((wT_jk == mu) & (sigma == 0))
  # {
  #   psi <- exp((5*sigma^2 + 2*sqrt(5)*gammaVal*(wT_jk - mu))/(2*gammaVal^2)) * ((t(E_1)%*%Lambda_61*pnorm(0)) + (t(E_1)%*%Lambda_62*(sigma/sqrt(2*pi))*exp(0))) - exp((5*sigma^2 - 2*sqrt(5)*gammaVal*(wT_jk - mu))/(2*gammaVal^2)) * ((t(E_2)%*%Lambda_71*pnorm(0)) + (t(E_2)%*%Lambda_72*(sigma/sqrt(2*pi))*exp(0)))
  # }
  return (psi)
}


iFun <- function(test, trainW, trainZ, gamma, mu, sigma)
{
  d <- length(trainW[1,])
  p <- length(trainZ[1,])
  m <- length(trainW[,1])
  xis <- matrix(1,nrow = m, ncol = d)
  cs <- matrix(1,nrow = m, ncol = p)
  iVec <- c(rep(1, each=m))
  for (i in 1:m) {
    for (k in 1:p) {
      cs[i,k] = cKSingle(test[d+k], trainZ[i,k], gamma[d+k])
      #iVec[i] = iVec[i] * cVal
    }
    for (k in 1:d) {
      xis[i,k] = xiSingle(mu[k], sigma[k], gamma[k], trainW[i,k])
      #iVec[i] = iVec[i] * xiVal
    }
  }
  for (k in 1:p) {
    iVec <- iVec*cs[,k]
  }
  for (k in 1:d) {
    iVec <- iVec*xis[,k]
  }
  
  return (iVec)
}


bFun <- function(test, trainW, trainZ, gamma, mu, sigma)
{
  d <- length(trainW[1,])
  p <- length(trainZ[1,])
  m <- length(trainW[,1])
  xis <- array(1,dim=c(d,m,d))
  cs <- array(1,dim=c(d,m,p))
  bMat <- matrix(1,nrow = d, ncol = m)
  for (l in 1:d) {
    for (j in 1:m) {
      bMat[l,j] = bMat[l,j] * psiSingle(mu[l], sigma[l], gamma[l], trainW[j,l]) 
      for (k in 1:d) {
        if (k != l) {
          xis[l,j,k] = xiSingle(mu[k], sigma[k], gamma[k], trainW[j,k])
          #bMat[l,j] = bMat[l,j] * xiVal
        }
        
      }
      for (k in 1:p) {
        cs[l,j,k] = cKSingle(test[d+k], trainZ[j,k], gamma[d+k])
        #bMat[l,j] = bMat[l,j] * cVal
      }
    }
  }
  for (k in 1:d) {
    bMat <- bMat*xis[,,k]
  }
  for (k in 1:p) {
    bMat <- bMat*cs[,,k]
  }
  return (bMat)
}

hZeta <- function(mu, sigma, gammaVal, x1, x2)
{
  mu_C <- mu - ((2*sqrt(5)*sigma^2)/(gammaVal))
  mu_D <- mu + ((2*sqrt(5)*sigma^2)/(gammaVal))
  Lambda_31 <- rbind(1, mu_C, (mu_C^2 + sigma^2), (mu_C^3 + 3*(sigma^2)*mu_C), (mu_C^4 + 6*(sigma^2)*(mu_C^2) + 3*sigma^4))
  Lambda_32 <- rbind(0, 1, (mu_C + x2), (mu_C^2 + 2*sigma^2 + x2^2 + mu_C*x2), (mu_C^3 + x2^3 + x2*mu_C^2 + mu_C*x2^2 + 3*(sigma^2)*x2 + 5*(sigma^2)*mu_C))
  Lambda_41 <- rbind(1, mu, (mu^2 + sigma^2), (mu^3 + 3*(sigma^2)*mu), (mu^4 + 6*(sigma^2)*(mu^2) + 3*sigma^4))
  Lambda_42 <- rbind(0, 1, (mu + x1), (mu^2 + 2*sigma^2 + x1^2 + mu*x1), (mu^3 + x1^3 + x1*mu^2 + mu*x1^2 + 3*(sigma^2)*x1 + 5*(sigma^2)*mu))
  Lambda_43 <- rbind(0, 1, (mu + x2), (mu^2 + 2*sigma^2 + x2^2 + mu*x2), (mu^3 + x2^3 + x2*mu^2 + mu*x2^2 + 3*(sigma^2)*x2 + 5*(sigma^2)*mu))
  Lambda_51 <- rbind(1, -mu_D, (mu_D^2 + sigma^2), (-mu_D^3 - 3*(sigma^2)*mu_D), (mu_D^4 + 6*(sigma^2)*(mu_D^2) +3*sigma^4))
  Lambda_52 <- rbind(0, 1, (-mu_D - x1), (mu_D^2 + 2*sigma^2 + x1^2 + mu_D*x1), (-mu_D^3 - x1^3 - x1*mu_D^2 - mu_D*x1^2 - 3*(sigma^2)*x1) - 5*(sigma^2)*mu_D)
  E_30 <- 1 + ((25*(x1^2)*(x2^2) - 3*sqrt(5)*(3*gammaVal^3 + 5*gammaVal*x1*x2)*(x1 + x2) + (15*gammaVal^2)*(x1^2 + x2^2 +3*x1*x2))/(9*gammaVal^4))
  E_31 <- ((18*sqrt(5)*gammaVal^3 + 15*sqrt(5)*gammaVal*(x1^2 + x2^2) - (75*gammaVal^2 + 50*x1*x2)*(x1 + x2) + 60*sqrt(5)*gammaVal*x1*x2)/(9*gammaVal^4))
  E_32 <- ((5*(5*x1^2 + 5*x2^2 + 15*gammaVal^2 - 9*sqrt(5)*gammaVal*(x1 + x2) + 20*x1*x2))/(9*gammaVal^4))
  E_33 <- ((10*(3*sqrt(5)*gammaVal - 5*x1 - 5*x2))/(9*gammaVal^4))
  E_34 <- ((25)/(9*gammaVal^4))
  E_40 <- 1 + ((25*(x1^2)*(x2^2) + 3*sqrt(5)*(3*gammaVal^3 - 5*gammaVal*x1*x2)*(x2 - x1) + (15*gammaVal^2)*(x1^2 + x2^2 - 3*x1*x2))/(9*gammaVal^4))
  E_41 <- ((5*(3*sqrt(5)*gammaVal*(x2^2 - x1^2) + (3*gammaVal^2)*(x1 + x2) - 10*x1*x2*(x1 + x2)))/(9*gammaVal^4))
  E_42 <- ((5*(5*x1^2 + 5*x2^2 -3*gammaVal^2 - 3*sqrt(5)*gammaVal*(x2 - x1) + 20*x1*x2))/(9*gammaVal^4))
  E_43 <- -((50*(x1 + x2))/(9*gammaVal^4))
  E_44 <- ((25)/(9*gammaVal^4))
  E_50 <- 1 + ((25*(x1^2)*(x2^2) + 3*sqrt(5)*(3*gammaVal^3 + 5*gammaVal*x1*x2)*(x1 + x2) + (15*gammaVal^2)*(x1^2 + x2^2 + 3*x1*x2))/(9*gammaVal^4))
  E_51 <- ((18*sqrt(5)*gammaVal^3 + 15*sqrt(5)*gammaVal*(x1^2 + x2^2) + (75*gammaVal^2 + 50*x1*x2)*(x1 + x2) + 60*sqrt(5)*gammaVal*x1*x2)/(9*gammaVal^4))
  E_52 <- ((5*(5*x1^2 + 5*x2^2 + 15*gammaVal^2 + 9*sqrt(5)*gammaVal*(x1 + x2) + 20*x1*x2))/(9*gammaVal^4))
  E_53 <- ((10*(3*sqrt(5)*gammaVal + 5*x1 + 5*x2))/(9*gammaVal^4))
  E_54 <- ((25)/(9*gammaVal^4))
  E_3 <- rbind(E_30, E_31, E_32, E_33, E_34)
  E_4 <- rbind(E_40, E_41, E_42, E_43, E_44)
  E_5 <- rbind(E_50, E_51, E_52, E_53, E_54)
  
  zetaV <- exp((10*sigma^2 +sqrt(5)*gammaVal*(x1 + x2 - 2*mu))/(gammaVal^2)) * (t(E_3)%*%Lambda_31*pnorm((mu_C - x2)/(sigma)) + t(E_3)%*%Lambda_32*((sigma)/(sqrt(2*pi)))*exp(-(((x2 - mu_C)^2)/(2*sigma^2)))) + exp(-((sqrt(5)*(x2-x1))/(gammaVal)))*(t(E_4)%*%Lambda_41*(pnorm((x2 - mu)/(sigma))-pnorm((x1-mu)/(sigma))) + t(E_4)%*%Lambda_42*((sigma)/(sqrt(2*pi)))*exp(-((x1-mu)^2/(2*sigma^2))) - t(E_4)%*%Lambda_43*((sigma)/(sqrt(2*pi)))*exp(-((x2-mu)^2/(2*sigma^2)))) + exp((10*sigma^2 - sqrt(5)*gammaVal*(x1+x2-2*mu))/(gammaVal^2))*(t(E_5)%*%Lambda_51*pnorm((x1-mu_D)/(sigma)) + t(E_5)%*%Lambda_52*((sigma)/(sqrt(2*pi)))*exp(-((x1-mu_D)^2/(2*sigma^2))))
  
  # if (sigma == 0) {
  #   if (x1 == mu) {
  #     zetaV <- exp((10*sigma^2 +sqrt(5)*gammaVal*(x1 + x2 - 2*mu))/(gammaVal^2)) * (t(E_3)%*%Lambda_31*pnorm((mu_C - x2)/(sigma)) + t(E_3)%*%Lambda_32*((sigma)/(sqrt(2*pi)))*exp(-(((x2 - mu_C)^2)/(2*sigma^2)))) + exp(-((sqrt(5)*(x2-x1))/(gammaVal)))*(t(E_4)%*%Lambda_41*(pnorm((x2 - mu)/(sigma))-pnorm(0)) + t(E_4)%*%Lambda_42*((sigma)/(sqrt(2*pi)))*exp(0) - t(E_4)%*%Lambda_43*((sigma)/(sqrt(2*pi)))*exp(-((x2-mu)^2/(2*sigma^2)))) + exp((10*sigma^2 - sqrt(5)*gammaVal*(x1+x2-2*mu))/(gammaVal^2))*(t(E_5)%*%Lambda_51*pnorm(0) + t(E_5)%*%Lambda_52*((sigma)/(sqrt(2*pi)))*exp(0))
  #   }
  #   if (x2 == mu) {
  #     zetaV <- exp((10*sigma^2 +sqrt(5)*gammaVal*(x1 + x2 - 2*mu))/(gammaVal^2)) * (t(E_3)%*%Lambda_31*pnorm(0) + t(E_3)%*%Lambda_32*((sigma)/(sqrt(2*pi)))*exp(0)) + exp(-((sqrt(5)*(x2-x1))/(gammaVal)))*(t(E_4)%*%Lambda_41*(pnorm(0)-pnorm((x1-mu)/(sigma))) + t(E_4)%*%Lambda_42*((sigma)/(sqrt(2*pi)))*exp(-((x1-mu)^2/(2*sigma^2))) - t(E_4)%*%Lambda_43*((sigma)/(sqrt(2*pi)))*exp(0)) + exp((10*sigma^2 - sqrt(5)*gammaVal*(x1+x2-2*mu))/(gammaVal^2))*(t(E_5)%*%Lambda_51*pnorm((x1-mu_D)/(sigma)) + t(E_5)%*%Lambda_52*((sigma)/(sqrt(2*pi)))*exp(-((x1-mu_D)^2/(2*sigma^2))))
  #   }
  #   if ((x1 == mu) & (x2 == mu)) {
  #     zetaV <- exp((10*sigma^2 +sqrt(5)*gammaVal*(x1 + x2 - 2*mu))/(gammaVal^2)) * (t(E_3)%*%Lambda_31*pnorm(0) + t(E_3)%*%Lambda_32*((sigma)/(sqrt(2*pi)))*exp(0)) + exp(-((sqrt(5)*(x2-x1))/(gammaVal)))*(t(E_4)%*%Lambda_41*(0) + t(E_4)%*%Lambda_42*((sigma)/(sqrt(2*pi)))*exp(0) - t(E_4)%*%Lambda_43*((sigma)/(sqrt(2*pi)))*exp(0)) + exp((10*sigma^2 - sqrt(5)*gammaVal*(x1+x2-2*mu))/(gammaVal^2))*(t(E_5)%*%Lambda_51*pnorm(0) + t(E_5)%*%Lambda_52*((sigma)/(sqrt(2*pi)))*exp(0))
  #   }
  # }
  
  
  return (zetaV)
}

zetaSingle <- function(mu, sigma, gamma, x1, x2)
{
  if (x2 >= x1) {
    zetaVal = hZeta(mu, sigma, gamma, x1, x2)
  } else {
    zetaVal = hZeta(mu, sigma, gamma, x2, x1)
  }
  return (zetaVal)
}


jFun <- function(test, trainW, trainZ, gamma, mu, sigma)
{
  d <- length(trainW[1,])
  p <- length(trainZ[1,])
  m <- length(trainW[,1])
  cVi <- array(1,dim=c(m,m,p))
  cVj <- array(1,dim=c(m,m,p))
  zetas <- array(1,dim=c(m,m,d))
  jMat <- matrix(1,nrow = m, ncol = m)
  for (i in 1:m) {
    for (j in 1:m) {
      for (k in 1:p) {
        cVi[i,j,k] = cKSingle(test[d+k], trainZ[i,k], gamma[d+k]) #consider making this parallel
        cVj[i,j,k] = cKSingle(test[d+k], trainZ[j,k], gamma[d+k])
        #jMat[i,j] = jMat[i,j] * cVali * cValj
      }
      for (k in 1:d) {
        zetas[i,j,k] = zetaSingle(mu[k], sigma[k], gamma[k], trainW[i,k], trainW[j,k]) #also here
        #jMat[i,j] = jMat[i,j] * zetaVal
      }
    }
  }
  for (k in 1:p) {
    jMat <- jMat*cVi[,,k]*cVj[,,k]
  }
  for (k in 1:d) {
    jMat <- jMat*zetas[,,k]
  }
  return (jMat)
}

## Squared Exponential
# cKSingle <- function(test, train, gammaVal)
# {
#   cor <- exp(-((test - train)^2/(gammaVal^2)))
#   return (cor)
# }
# xiSingle <- function(mu, sigma, gammaVal, wT_ik)
# {
#   xi <- (1/(sqrt(1 + ((2*sigma^2)/(gammaVal^2)))))*exp(-((mu - wT_ik)^2/(2*sigma^2 + gammaVal^2)))
#   return (xi)
# }
# 
# psiSingle <- function(mu, sigma, gammaVal, wT_jk)
# {
#   psi <- (1/(sqrt(1 + ((2*sigma^2)/(gammaVal^2)) )))*exp(-((mu - wT_jk)^2/(2*sigma^2 + gammaVal^2)))*((2*(sigma^2)*wT_jk + (gammaVal^2)*mu)/(2*sigma^2 + gammaVal^2))
#   return (psi)
# }
# 
# zetaSingle <- function(mu, sigma, gammaVal, wT_ik, wT_jk)
# {
#   zeta <- (1/(sqrt(1 + (4*sigma^2)/(gammaVal^2))))*exp(-((((wT_ik + wT_jk)/2) - mu)^2/(((gammaVal^2)/2) + 2*sigma^2))-((wT_ik - wT_jk)^2/(2*gammaVal^2)))
#   return (zeta)
# }

##Generate design##
set.seed(9)
##First level of emulation##
library(readxl)
library(ggplot2)
library(lhs)
library(matlib)
library(plot3D)
library(pracma)
library(RobustGaSP)
library(svMisc)

num_inputs = 3
num_outputs = 1

#pull files from github repo to run this
training_input_raw <- read_excel("bTrainIn.xlsx",  col_names = FALSE, col_types = "numeric")
training_output_raw <- read_excel("bTrainOut.xlsx", col_names = FALSE, col_types = "numeric")

pts <- training_input_raw$...3[1:150]

training_input <- cbind(training_input_raw$...1, training_input_raw$...2, training_input_raw$...3, training_output_raw$...1)
#find closest indices
indices <- pts*0
lp <- length(pts)

for (i in 1:lp) {
  temp <- abs(training_input[,3] - pts[i])
  indices[i] <- which.min(temp)
}
numCurves <- 6
for (i in 1:numCurves) {

  for (j in 1:lp) {
    if ((i == 1) && (j == 1)) {
      trainData <- training_input[indices[1],]
    }
    else {
      temp <- training_input[((i-1)*150 + indices[j]),]
      trainData <- rbind(trainData, temp)
    }
  }
  if (i == 1) {
    trainInputs <- trainData[((i-1)*lp+1),1:2]
    trainOutputs <- trainData[,4]
  } else {
    trainInputs <- rbind(trainInputs, trainData[((i-1)*lp+1),1:2])
    trainOutputs <- rbind(trainOutputs, trainData[((i-1)*lp+1):(i*lp),4])
  }
  
}
xVals <- trainData[,3]
yVals <- trainData[,4]
trainData <- trainData[,1:2]
m <- dim(trainInputs)[1]
trendMat <- cbind(trainInputs,matrix(1,m,1))
gpModel <- ppgasp(design=trainInputs, response=trainOutputs, trend=trendMat)
fGammas <- 1/gpModel@beta_hat

#see files in github repo for these files
testing_input_raw1 <- read_excel("bTestOut.xlsx",  col_names = FALSE, col_types = "numeric")
testing_input_raw2 <- read_excel("bTestIn.xlsx",  col_names = FALSE, col_types = "numeric")

temp <- cbind(testing_input_raw1$...1, testing_input_raw2$...1)

test_all <- temp
gTestInputs <- test_all
numTest <- dim(test_all)[1]

testInputs <- cbind(0.947788063691687, 0.0439093869169055)

testTrendMat <- cbind(testInputs, 1)
pred <- predict.ppgasp(gpModel, testInputs, testing_trend = testTrendMat)

# ##Second level of emulation##
d = 1
p = 1
g = 1

wT <- rbind(yVals[15], yVals[80], yVals[135])
zT <- rbind(pts[15], pts[80], pts[135])
for (i in 2:numCurves) {
  wT <- rbind(wT, yVals[(i-1)*150 + 15], yVals[(i-1)*150 + 80], yVals[(i-1)*150 + 135])
  zT <- rbind(zT, pts[15], pts[80], pts[135])
}
yT <- wT
in_all <- cbind(wT, zT)
zT <- in_all[,2]
w_all <- in_all
g_all <- yT
m <- length(wT)
trendMat <- matrix(1,m,1)
gpModel2 <- ppgasp(design=w_all, response=as.matrix(g_all), trend = trendMat)

gTestInputs <- test_all
numTest <- dim(test_all)[1]
testTrendMat <- matrix(1,numTest,1)
pred2 <- predict(gpModel2, test_all, testing_trend = testTrendMat)

muVec <- pred$mean
sigVec <- pred$sd

in_all <- cbind(wT, zT)
yT <- wT
w_all <- in_all
g_all <- yT
m <- dim(in_all)[1]
#
eta <- gpModel2@nugget #nugget term
gSigVec <- gpModel2@sigma2_hat
gamma <- gpModel2@beta_hat
gamma <- 1/gamma

ndp=dim(w_all)[1]
numVars = dim(w_all)[2]
B=matrix(1, nrow = ndp, ncol = ndp)
for (j in 1:ndp) {
  for (i in 1:ndp) {
    for (k in 1:numVars) {
      B[i,j] = B[i,j] * cKSingle(w_all[i,k], w_all[j,k], gamma[k])
    }
  }
}
B=B+(0.00001)^2 * diag(ndp)
saveB <- B

c_k <- array(0, dim = c(m, m, (d+p)))
for (i in 1:m) {
  for (j in 1:m) {
    for (k in 1:(d+p)) {
      c_k[i,j,k] <- cKSingle(in_all[i,k], in_all[j,k], gamma[k])
    }
  }
}

R <- cFunction(c_k)
for (i in 1:m) {
  for (j in 1:m) {
    if (i == j) {
      R[i,j] <- R[i,j] + eta + 0.000001
    }
  }
  
}
invR <- solve(R)
HzT <- matrix(1,m,1)
hTilde <- matrix(1,m,1)
hats <- inv(t(hTilde) %*% invR %*% hTilde)%*%t(hTilde) %*% invR %*% yT
#normally you would use the commented bits here 
#but for the trend we are using we use the uncommented lines
# thetaHat <- hats[1:d,]
# betaHat <- hats[(d+1):(d+p+1),]
# thetaHat <- gpModel2@theta_hat[1:d]
# betaHat <- gpModel2@theta_hat[(d+1):(d+p+1)]
thetaGat <- 0
betaHat <- gpModel2@theta_hat
B <- matrix(0,nrow = d, ncol = m)
Q <- invR %*% hTilde %*% inv(t(hTilde) %*% invR %*% hTilde) %*% t(hTilde) %*% invR - invR
C <- inv(t(hTilde) %*% invR %*% hTilde)
muLs <- matrix(c(0), nrow = numTest, ncol = g)
sigma2Ls <- matrix(c(0), nrow = numTest, ncol = g)
saveC <- matrix(0, nrow = numTest, ncol = m)
for (i in 1:g) {
  A <- invR%*%(yT - wT * thetaHat - HzT %*% betaHat)
  for(k in 1:numTest){
    print(k)
    Sys.sleep(0.001)
    flush.console()
    hZ <- 1
    iVec <- iFun(gTestInputs[k,], as.matrix(wT), as.matrix(zT), gamma, muVec[,k], sigVec[k])
    B <- bFun(gTestInputs[k,], as.matrix(wT), as.matrix(zT), gamma, muVec[,k], sigVec[k])
    K <- cbind(t(B), iVec%*%t(hZ))
    J <- jFun(gTestInputs[k,], as.matrix(wT), as.matrix(zT), gamma, muVec[,k], sigVec[k])
    Omega <- sigVec[k]^2
    Omega <- as.matrix(Omega)
    P <- blkdiag(Omega, matrix(0,nrow = length(hZ), ncol = length(hZ)))
    G <- cbind(t(muVec[,k]), t(hZ))
    G <- t(G)
    mu_L <- (t(muVec[,k]) * thetaHat) + (t(hZ) %*% betaHat) + (t(iVec) %*% A)
    
    sigma2_L <- t(A) %*% (J - iVec%*%t(iVec)) %*% A + gSigVec[i]*(1 + eta + inv(t(hTilde)%*%invR%*%hTilde) + tr(Q%*%J) - 2*tr(C * t(hTilde) %*% invR %*% iVec))
    muLs[k, i] <- mu_L
    sigma2Ls[k, i] <- sigma2_L
    if(sigma2_L < 0) {
      sigma2Ls[k,i] = abs(sigma2_L)
    }
    saveC[k,] <- iVec
  }
}
N <- numTest
A=matrix(1, nrow = N, ncol = N)
C <- saveC
C2=matrix(1, nrow = N, ncol = ndp)

for (j in 1:N) {
  for (i in 1:N) {
    if (i == j) {
      A[i,j] = A[i,j] * sigma2Ls[i]
    }
    else if (i != j) {
      A[i,j] = A[i,j] * sqrt(sigma2Ls[i]) * sqrt(sigma2Ls[j]) * cKSingle(muLs[i], muLs[j], gamma[1]) *  cKSingle(test_all[i,2], test_all[j,2], gamma[2])
      # A[i,j] = A[i,j] * sqrt(sigma2Ls[i]) * sqrt(sigma2Ls[j]) * rhoCalc(x_i, x_j, fGammas, sigma_i, sigma_j) *  cKSingle(test_all[i,2], test_all[j,2], gamma[2])
    }
  }
  for (i in 1:ndp) {
    for (k in 1:numVars) {
      C2[j,i] = C2[j,i] * cKSingle(w_all[i,k], test_all[j,k], gamma[k])
    }
  }
}

R = A + (0.000001 * diag(N))
L = chol(R)
L = t(L)
samples <- matrix(0, nrow = numTest, ncol = 100)
beta1 <- thetaHat
beta2 <- 0
beta3 <- betaHat
for (i in 1:100) {
  u = rnorm(N)
  ###
  ysamp = as.vector(beta3) + C %*% (inv(saveB)%*%(g_all - as.vector(beta3))) + L%*%u
  
  samples[,i] <- ysamp
}
plotXs <- matrix(rep(pts, times=100), nrow = N)
plotRuns <- t(matrix(rep(1:100, times=N), ncol = N))
df <- data.frame(x = as.vector(plotXs), y = as.vector(samples), runId = as.vector(plotRuns))

g2 <- ggplot(df,aes(x = x, y = y)) + geom_line(aes(color=factor(runId)))
df2 <- data.frame(x = gTestInputs[,2], mu = muLs, lb = muLs - 1.96*sqrt(sigma2Ls), ub = muLs + 1.96*sqrt(sigma2Ls), trueVals = gTestInputs[,1], ppgaspPred = pred2$mean, plb = pred2$lower95, pub = pred2$upper95)
g2 <- g2 + geom_line(data = df2, aes(x = x, y = trueVals), lwd=1, color='blue')
g2 <- g2 + geom_line(data = df2, aes(x = x, y = mu), lwd=1)
g2 <- g2 + geom_line(data = df2, aes(x = x, y = lb), lwd=1, linetype = "dashed")
g2 <- g2 + geom_line(data = df2, aes(x = x, y = ub), lwd=1, linetype = "dashed")
g2 <- g2 + theme(legend.position = "none")
g2

g3 <- ggplot(df,aes(x = x, y = y)) + geom_line(aes(color=factor(runId)))
df22 <- data.frame(xVals = pts, predVals = t(pred$mean), l95 = t(pred$lower95), u95 = t(pred$upper95))
g3 <- g3 + geom_line(data = df22, aes(x = xVals, y = predVals), color = "black", lwd=1)
g3 <- g3 + geom_line(data = df22, aes(x = xVals, y = l95), color = "black", lwd=1, linetype = "dashed")
g3 <- g3 + geom_line(data = df22, aes(x = xVals, y = u95), color = "black", lwd=1, linetype = "dashed")
g3 <- g3 + theme(legend.position = "none")
g3

g22 <- ggplot(df,aes(x = x, y = y)) + geom_line(aes(color=factor(runId)))
g22 <- g22 + geom_line(data = df2, aes(x = x, y = ppgaspPred), lwd=1, color='black')
g22 <- g22 + geom_line(data = df2, aes(x = x, y = plb), lwd=1, linetype = "dashed", color='black')
g22 <- g22 + geom_line(data = df2, aes(x = x, y = pub), lwd=1, linetype = "dashed", color='black')
g22 <- g22 + theme(legend.position = "none")
g22
