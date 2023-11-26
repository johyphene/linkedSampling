f <- function(x,z)
{
  return(cos(3+5*x*sin(x))*exp(-z))
}

mu2 <- function(b1, b2, b3, w, z)
{
  return (b3+as.vector(b1)*w + as.vector(b2)*z)
}

#Matern 5/2 Calculations
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
    }
    for (k in 1:d) {
      xis[i,k] = xiSingle(mu[k], sigma[k], gamma[k], trainW[i,k])
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
        }
        
      }
      for (k in 1:p) {
        cs[l,j,k] = cKSingle(test[d+k], trainZ[j,k], gamma[d+k])
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
        cVi[i,j,k] = cKSingle(test[d+k], trainZ[i,k], gamma[d+k])
        cVj[i,j,k] = cKSingle(test[d+k], trainZ[j,k], gamma[d+k])
      }
      for (k in 1:d) {
        zetas[i,j,k] = zetaSingle(mu[k], sigma[k], gamma[k], trainW[i,k], trainW[j,k])
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

## Squared Exponential Calculations
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

#create test data
x <- matrix(seq(-2,2, by = 0.5))
z <- matrix(seq(0,3, by = 0.375))

testInputs <- x
testOutputs <- matrix(0, nrow = length(x), ncol = length(z))
testOutputs2 <- matrix(0, nrow = length(x)*length(z), ncol = 3)
for (i in 1:length(x)) {
  for (j in 1:length(z)) {
    testOutputs[i,j] = f(x[i],z[j])
    testOutputs2[((i-1)*length(x) + j),1] = f(x[i],z[j])
    testOutputs2[((i-1)*length(x) + j),2] = x[i]
    testOutputs2[((i-1)*length(x) + j),3] = z[j]
  }
}

#create train data
pts <- matrix(seq(-1.9, 1.9, by = .38))
pts <- pts + 0.025
z2 <- z

trainInputs <- pts
trainOutputs <- matrix(0, nrow = length(pts), ncol = length(z2))
trainOutputs2 <- matrix(0, nrow = length(pts)*length(z2), ncol = 3)
for (i in 1:length(pts)) {
  for (j in 1:length(z2)) {
    trainOutputs[i,j] = f(pts[i],z2[j])
    trainOutputs2[((i-1)*length(x) + j),1] = f(pts[i],z2[j])
    trainOutputs2[((i-1)*length(x) + j),2] = pts[i]
    trainOutputs2[((i-1)*length(x) + j),3] = z2[j]
  }
}

numCurves <- length(pts)

xVals <- pts
yVals <- trainOutputs
trainData <- trainInputs
m <- size(trainInputs)[1]
trendMat <- cbind(trainInputs, matrix(1,m,1))
gpModel <- ppgasp(design=trainInputs, response=trainOutputs, trend=trendMat)

fGamma <- 1/gpModel@beta_hat

numTest <- dim(testInputs)[1]

testTrendMat <- cbind(testInputs, 1)
pred <- predict.ppgasp(gpModel, testInputs, testing_trend = testTrendMat)

##Second level of emulation##
d = 1
p = 1
g = 1

wT <- matrix(t(trainOutputs))
zT <- matrix(z)

yT <- wT
in_all <- cbind(wT, as.vector(zT))
zT <- in_all[,2]
w_all <- in_all
g_all <- yT
m <- length(wT)
trendMat <- matrix(1,m,1)
gpModel2 <- ppgasp(design=w_all, response=as.matrix(g_all), trend = trendMat)

test_all <- cbind(matrix(t(testOutputs)),as.vector(z))
gTestInputs <- test_all
numTest <- size(test_all)[1]
testTrendMat <- matrix(1,numTest,1)
pred2 <- predict(gpModel2, test_all, testing_trend = testTrendMat)
muVec <- t(matrix(pred$mean))
sigVec <- t(matrix(pred$sd))

##linked emulation begins##
##Calculations follow process outlined in 
##Linked Gaussian Process Emulation for Systems of Computer Models Using
##Matern Kernels and Adaptive Design
##By Deyu Ming and Serge Guillas

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
      R[i,j] <- R[i,j] + eta + 0.000000001
    }
  }
  
}
invR <- solve(R)
HzT <- matrix(1,m,1)
hTilde <- matrix(1,m,1)
hats <- inv(t(hTilde) %*% invR %*% hTilde)%*%t(hTilde) %*% invR %*% yT
thetaHat <- 0
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

##begin spatially correlated sampling##
N <- numTest
A=matrix(1, nrow = N, ncol = N)
C <- saveC
C2=matrix(1, nrow = N, ncol = ndp)
numInputs <- length(testInputs)
for (j in 1:N) {
  for (i in 1:N) {
    if (i == j) {
      A[i,j] = A[i,j] * sigma2Ls[i]
    }
    else if (i != j) {
      if (i %% numInputs == 0)
        indexI = i/numInputs
      else
        indexI = floor(i/numInputs)+1
      if (j %% numInputs == 0)
        indexJ = j/numInputs
      else
        indexJ = floor(j/numInputs)+1
      
      rho_f = cKSingle(testInputs[indexI], testInputs[indexJ], fGamma[1])
      tempSig = sqrt(sigma2Ls[i] - 2*rho_f*sqrt(sigma2Ls[i])*sqrt(sigma2Ls[j]) + sigma2Ls[j])
      rho_eta = xiSingle(muLs[i], tempSig, gamma[1], muLs[j])
      A[i,j] = A[i,j] * sqrt(sigma2Ls[i]) * sqrt(sigma2Ls[j]) * rho_eta * cKSingle(test_all[i,2], test_all[j,2], gamma[2])
    }
  }
}

R = A + (0.0000000001 * diag(N))
L = chol(R)
L = t(L)
samples <- matrix(0, nrow = numTest, ncol = 100)
beta1 <- thetaHat
beta2 <- 0
beta3 <- betaHat
for (i in 1:100) {
  u = rnorm(N)
  #this calculation is a bit different because no trend is in the second ppgasp emulator
  ysamp = as.vector(beta3) + C %*% (inv(saveB)%*%(g_all - as.vector(beta3))) + L%*%u
  samples[,i] <- ysamp
}


##begin plotting##
plotX <- testOutputs2[,2] #x
plotX <- matrix(testOutputs2[,2], nrow = length(x), byrow = TRUE)
plotY <- testOutputs2[,3] #z
plotY <- matrix(testOutputs2[,3], nrow = length(x), byrow = TRUE)
plotZ <- testOutputs2[,1]  #y
plotZ <- matrix(muLs, nrow = length(x), byrow = TRUE)
ub <- muLs + 1.96*sqrt(sigma2Ls)
lb <- muLs - 1.96*sqrt(sigma2Ls)
plotUB <- matrix(ub, nrow = length(x), byrow = TRUE)
plotLB <- matrix(lb, nrow = length(x), byrow = TRUE)
plotSamp <- matrix(samples[,2], nrow = length(x), byrow = TRUE)


library(plotly)
fig1 <- plot_ly(showscale = FALSE)
fig1 <- fig1 %>% add_surface(x = ~plotX, y = ~plotY, z = ~plotZ)
fig1 <- fig1 %>% add_surface(x = ~plotX, y = ~plotY, z = ~plotLB, opacity = 0.5)
fig1 <- fig1 %>% add_surface(x = ~plotX, y = ~plotY, z = ~plotUB, opacity = 0.5)
fig1

plotSamp1 <- matrix(samples[,1], nrow = length(x), byrow = TRUE)

scene = list(aspectmode = "cube", camera = list(eye = list(x = -2.25, y = 0.75, z = 0.75)), xaxis = list(title = 'x'), yaxis = list(title = 'z'), zaxis = list(title = 'f(x,z)'))

fig2 <- plot_ly(showscale = FALSE)
fig2 <- fig2 %>% add_surface(x = ~plotX, y = ~plotY, z = ~plotSamp1)
fig2 <- fig2 %>% add_surface(x = ~plotX, y = ~plotY, z = ~plotLB, opacity = 0.45)
fig2 <- fig2 %>% add_surface(x = ~plotX, y = ~plotY, z = ~plotUB, opacity = 0.45)
fig2 <- fig2 %>% layout(scene = scene)
fig2
