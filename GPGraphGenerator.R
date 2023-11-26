library(lhs)
library(matlib)
library(plot3D)
library(pracma)
library(RobustGaSP)
library(svMisc)

f1 <- function(x)
{
  return (3 + 5*x*sin(x))
}

f <- function(x)
{
  return (2*x+cos(5*x))
}

mu <- function(b1, b2, x)
{
  return (b1+(b2)*x)
}

mu2 <- function(b1, b2, b3, w, z)
{
  return (b3+as.vector(b1)*w + as.vector(b2)*z)
}

ypred <- function(z, y, x, gammaVal, b1, b2)
{
  # (r(z))'*(B\(y-mu(x)))
  return (t(cKSingle(z, x, gammaVal)) %*% (inv(B)%*%(y - mu(b1, b2, x))))
}

g1 <- function(w1, z1)
{
  return (cos(w1) * (exp(-1*z1)))
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
  
  if ((wT_ik == mu) & (sigma == 0))
  {
    xi <- exp((5*sigma^2 + 2*sqrt(5)*gammaVal*(wT_ik - mu))/(2*gammaVal^2)) * ((t(E_1)%*%Lambda_11*pnorm(0)) + (t(E_1)%*%Lambda_12*(sigma/sqrt(2*pi))*exp(0))) + exp((5*sigma^2 - 2*sqrt(5)*gammaVal*(wT_ik - mu))/(2*gammaVal^2)) * ((t(E_2)%*%Lambda_21*pnorm(0)) + (t(E_2)%*%Lambda_22*(sigma/sqrt(2*pi))*exp(0)))
  }
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
  
  if ((wT_jk == mu) & (sigma == 0))
  {
    psi <- exp((5*sigma^2 + 2*sqrt(5)*gammaVal*(wT_jk - mu))/(2*gammaVal^2)) * ((t(E_1)%*%Lambda_61*pnorm(0)) + (t(E_1)%*%Lambda_62*(sigma/sqrt(2*pi))*exp(0))) - exp((5*sigma^2 - 2*sqrt(5)*gammaVal*(wT_jk - mu))/(2*gammaVal^2)) * ((t(E_2)%*%Lambda_71*pnorm(0)) + (t(E_2)%*%Lambda_72*(sigma/sqrt(2*pi))*exp(0)))
  }
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
  for (k in 1:p) {
    bMat <- bMat*xis[,,k]
  }
  for (k in 1:d) {
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
  
  if (sigma == 0) {
    if (x1 == mu) {
      zetaV <- exp((10*sigma^2 +sqrt(5)*gammaVal*(x1 + x2 - 2*mu))/(gammaVal^2)) * (t(E_3)%*%Lambda_31*pnorm((mu_C - x2)/(sigma)) + t(E_3)%*%Lambda_32*((sigma)/(sqrt(2*pi)))*exp(-(((x2 - mu_C)^2)/(2*sigma^2)))) + exp(-((sqrt(5)*(x2-x1))/(gammaVal)))*(t(E_4)%*%Lambda_41*(pnorm((x2 - mu)/(sigma))-pnorm(0)) + t(E_4)%*%Lambda_42*((sigma)/(sqrt(2*pi)))*exp(0) - t(E_4)%*%Lambda_43*((sigma)/(sqrt(2*pi)))*exp(-((x2-mu)^2/(2*sigma^2)))) + exp((10*sigma^2 - sqrt(5)*gammaVal*(x1+x2-2*mu))/(gammaVal^2))*(t(E_5)%*%Lambda_51*pnorm(0) + t(E_5)%*%Lambda_52*((sigma)/(sqrt(2*pi)))*exp(0))
    }
    if (x2 == mu) {
      zetaV <- exp((10*sigma^2 +sqrt(5)*gammaVal*(x1 + x2 - 2*mu))/(gammaVal^2)) * (t(E_3)%*%Lambda_31*pnorm(0) + t(E_3)%*%Lambda_32*((sigma)/(sqrt(2*pi)))*exp(0)) + exp(-((sqrt(5)*(x2-x1))/(gammaVal)))*(t(E_4)%*%Lambda_41*(pnorm(0)-pnorm((x1-mu)/(sigma))) + t(E_4)%*%Lambda_42*((sigma)/(sqrt(2*pi)))*exp(-((x1-mu)^2/(2*sigma^2))) - t(E_4)%*%Lambda_43*((sigma)/(sqrt(2*pi)))*exp(0)) + exp((10*sigma^2 - sqrt(5)*gammaVal*(x1+x2-2*mu))/(gammaVal^2))*(t(E_5)%*%Lambda_51*pnorm((x1-mu_D)/(sigma)) + t(E_5)%*%Lambda_52*((sigma)/(sqrt(2*pi)))*exp(-((x1-mu_D)^2/(2*sigma^2))))
    }
    if ((x1 == mu) & (x2 == mu)) {
      zetaV <- exp((10*sigma^2 +sqrt(5)*gammaVal*(x1 + x2 - 2*mu))/(gammaVal^2)) * (t(E_3)%*%Lambda_31*pnorm(0) + t(E_3)%*%Lambda_32*((sigma)/(sqrt(2*pi)))*exp(0)) + exp(-((sqrt(5)*(x2-x1))/(gammaVal)))*(t(E_4)%*%Lambda_41*(0) + t(E_4)%*%Lambda_42*((sigma)/(sqrt(2*pi)))*exp(0) - t(E_4)%*%Lambda_43*((sigma)/(sqrt(2*pi)))*exp(0)) + exp((10*sigma^2 - sqrt(5)*gammaVal*(x1+x2-2*mu))/(gammaVal^2))*(t(E_5)%*%Lambda_51*pnorm(0) + t(E_5)%*%Lambda_52*((sigma)/(sqrt(2*pi)))*exp(0))
    }
  }
  
  
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

# # Squared Exponential
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

set.seed(9)
m = 100
x = randomLHS(m,1)
x = 6*x - 3 #for [-3, 3]
y=f1(x)

trendMat <- cbind(x,matrix(1,m,1))
#use top line for squared exponential case
# gpModel <- ppgasp(design=x, response=as.matrix(y), trend = trendMat, kernel_type = 'pow_exp')
gpModel <- ppgasp(design=x, response=as.matrix(y), trend = trendMat)

beta1 <- gpModel@theta_hat[2] #intercept
beta2 <- gpModel@theta_hat[1] #x

gamma <- gpModel@beta_hat
gamma <- 1/gamma

ndp=length(x);
B=matrix(0, nrow = ndp, ncol = ndp)
for (j in 1:ndp) {
  for (i in 1:ndp) {
    B[i,j]=cKSingle(x[i], x[j], gamma);
  }
}
B=B+(0.00001)^2 * diag(ndp);

N=51
xx2 = seq(-3,3, length.out = N)
xx <- rep(xx2, times = N)
N = 2601
testTrendMat <- cbind(xx,matrix(1,N,1))
pred <- predict(gpModel, as.matrix(xx), testing_trend = testTrendMat)

ygp = pred$mean
v = pred$sd
N <- 51
w <- y
wTest <- ygp
z = randomLHS(m,1)
z = 6*z - 3 #for [-3, 3]
zT <- seq(-3,3, length.out = N)
zT2 <- zT
zTest <- t(zT)
for (i in 1:(N-1)) {
  zTest <- rbind(zTest, zT)
}
zTest <- as.vector(zTest)

w_all <- cbind(w, z)
g_all <- g1(w, z)
test_all <- cbind(wTest, zTest)

trendMat <- cbind(w_all,matrix(1,m,1))

#use top line for squared exponential case
# gpModel2 <- ppgasp(design=w_all, response=as.matrix(g_all), trend = trendMat, kernel_type = 'pow_exp')
gpModel2 <- ppgasp(design=w_all, response=as.matrix(g_all), trend = trendMat)

beta1 <- gpModel2@theta_hat[1] #w
beta2 <- gpModel2@theta_hat[2] #z
beta3 <- gpModel2@theta_hat[3] #intercept

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

d = 1
p = 1
g = 1
numTest <- dim(test_all)[1]
wT <- y
zT <- z
yT <- g_all
in_all <- w_all
gTestInputs <- test_all
testTrendMat <- cbind(test_all,matrix(1,numTest,1))
pred2 <- predict(gpModel2, test_all, testing_trend = testTrendMat)

muVec <- t(pred$mean)
sigVec <- pred$sd

eta <- gpModel2@nugget #nugget term
gSigVec <- gpModel2@sigma2_hat
gamma <- gpModel2@beta_hat
gamma <- 1/gamma
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
      R[i,j] <- R[i,j] + eta
    }
  }
  
}
HzT <- cbind(zT,matrix(1,m,1))
hTilde <- cbind(wT, HzT)
hats <- inv(t(hTilde) %*% inv(R) %*% hTilde)%*%t(hTilde) %*% inv(R) %*% yT
thetaHat <- hats[1:d,]
betaHat <- hats[(d+1):(d+p+1),]
thetaHat2 <- gpModel2@theta_hat[1:d]
betaHat2 <- gpModel2@theta_hat[(d+1):(d+p+1)]
B <- matrix(0,nrow = d, ncol = m)
Q <- inv(R) %*% hTilde %*% inv(t(hTilde) %*% inv(R) %*% hTilde) %*% t(hTilde) %*% inv(R) - inv(R)
C <- inv(t(hTilde) %*% inv(R) %*% hTilde)
muLs <- matrix(c(0), nrow = numTest, ncol = g)
sigma2Ls <- matrix(c(0), nrow = numTest, ncol = g)
saveC <- matrix(0, nrow = numTest, ncol = ndp)
for (i in 1:g) {
  A <- inv(R)%*%(yT - wT * thetaHat - HzT %*% betaHat)
  for(k in 1:numTest){
    print(k)
    Sys.sleep(0.01)
    flush.console()
    hZ <- c(gTestInputs[k,2],1)
    iVec <- iFun(gTestInputs[k,], as.matrix(wT), as.matrix(zT), gamma, muVec[,k], sigVec[k,])
    B <- bFun(gTestInputs[k,], as.matrix(wT), as.matrix(zT), gamma, muVec[,k], sigVec[k,])
    K <- cbind(t(B), iVec%*%t(hZ))
    J <- jFun(gTestInputs[k,], as.matrix(wT), as.matrix(zT), gamma, muVec[,k], sigVec[k,])
    Omega <- sigVec[k,]^2
    Omega <- as.matrix(Omega)
    P <- blkdiag(Omega, matrix(0,nrow = length(hZ), ncol = length(hZ)))
    G <- cbind(t(muVec[,k]), t(hZ))
    G <- t(G)
    mu_L <- (t(muVec[,k]) %*% thetaHat) + (t(hZ) %*% betaHat) + (t(iVec) %*% A)
    
    sigma2_L <- t(A) %*% (J - iVec%*%t(iVec)) %*% A + 2*t(thetaHat)%*%(B - muVec[,k] %*% t(iVec)) %*% A + tr(thetaHat %*% t(thetaHat) %*% Omega) + gSigVec[i] *(1 + eta + tr(Q%*%J) + t(G) %*% C %*% G + tr(C %*% P - 2*C %*% t(hTilde) %*% inv(R) %*% K))
    muLs[k, i] <- mu_L
    sigma2Ls[k, i] <- sigma2_L
    if(sigma2_L < 0) {
      sigma2Ls[k,i] = abs(sigma2_L)
    }
    saveC[k,] <- iVec
  }
}

B <- saveB
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

for (i in 1:100) {
  u = rnorm(N)
  
  ###
  ysamp = mu2(beta1, beta2, beta3, test_all[,1], test_all[,2]) + C %*% (inv(B)%*%(g_all - mu2(beta1, beta2, beta3, w_all[,1], w_all[,2]))) + L%*%u
  samples[,i] <- ysamp
}
vars <- c(0)
means <- c(0)
for (i in 1:numTest) {
  vars[i] <- var(samples[i,])
  means[i] <- mean(samples[i,])
}
#~~~#
#COMPOSITE EMULATOR#
compTrainVals <- cbind(x, z)
compTrendMat <- cbind(compTrainVals,matrix(1,m,1))
compY <- g1(f1(x), z)
compositeModel <- ppgasp(design=compTrainVals, response=compY, trend = compTrendMat)
compTestVals <- cbind(xx, zTest)
compTestTrendMat <- cbind(compTestVals,matrix(1,numTest,1))
compPred <- predict(compositeModel, compTestVals, testing_trend = compTestTrendMat)
#~~~#
trueOutputs <- g1(test_all[,1], test_all[,2])
tw <- matrix(xx, nrow = 51, byrow = TRUE)
tz <- matrix(test_all[,2], nrow = 51, byrow = TRUE)
ty <- matrix(ysamp, nrow = 51, byrow = TRUE)
tvars <- matrix(vars, nrow = 51, byrow = TRUE)
gpvar <- (pred2$sd)^2
muvar <- matrix(pred2$mean, nrow = 51, byrow = TRUE)
gpVars <- matrix(gpvar, nrow = 51, byrow = TRUE)
manualVars <- matrix(sigma2Ls[,1], nrow = 51, byrow = TRUE)

#looking at a slice of the surface
library(ggplot2)
startpt <- 1735
endpt <- 1785
plotXs <- matrix(rep(tw[1,], times=100), nrow = 51)
plotRuns <- t(matrix(rep(1:100, times=51), ncol = 51))
plotVars <- samples[startpt:endpt,]
trueOuts <- trueOutputs[startpt:endpt]

df <- data.frame(x = as.vector(plotXs), y = as.vector(plotVars), runId = as.vector(plotRuns))
g01 <- ggplot(df,aes(x = x, y = y)) + geom_line(aes(color=factor(runId)))
df2 <- data.frame(x = tw[1,], mu = compPred$mean[startpt:endpt], lb = compPred$lower95[startpt:endpt], ub = compPred$upper95[startpt:endpt], trueVals = trueOuts)
g01 <- g01 + geom_line(data = df2, aes(x = x, y = trueVals), lwd=1, color='blue')
g01 <- g01 + geom_line(data = df2, aes(x = x, y = mu), lwd=1)
g01 <- g01 + geom_line(data = df2, aes(x = x, y = lb), lwd=1, linetype = "dashed")
g01 <- g01 + geom_line(data = df2, aes(x = x, y = ub), lwd=1, linetype = "dashed")
g01 <- g01 + theme(legend.position = "none")
g01 <- g01 + xlab("x") + ylab("g")
g01


g2 <- ggplot(df,aes(x = x, y = y)) + geom_line(aes(color=factor(runId)))
df2 <- data.frame(x = tw[1,], mu = pred2$mean[startpt:endpt], lb = pred2$lower95[startpt:endpt], ub = pred2$upper95[startpt:endpt], trueVals = trueOuts)
g2 <- g2 + geom_line(data = df2, aes(x = x, y = trueVals), lwd=1, color='blue')
g2 <- g2 + geom_line(data = df2, aes(x = x, y = mu), lwd=1)
g2 <- g2 + geom_line(data = df2, aes(x = x, y = lb), lwd=1, linetype = "dashed")
g2 <- g2 + geom_line(data = df2, aes(x = x, y = ub), lwd=1, linetype = "dashed")
g2 <- g2 + theme(legend.position = "none")
g2 <- g2 + xlab("x") + ylab("g")
g2

g3 <- ggplot(df,aes(x = x, y = y)) + geom_line(aes(color=factor(runId)))
df2 <- data.frame(x = tw[1,], mu = muLs[startpt:endpt], lb = muLs[startpt:endpt] - 1.96*sqrt(sigma2Ls[startpt:endpt]), ub = muLs[startpt:endpt] + 1.96*sqrt(sigma2Ls[startpt:endpt]), trueVals = trueOuts)
g3 <- g3 + geom_line(data = df2, aes(x = x, y = trueVals), lwd=1, color='blue')
g3 <- g3 + geom_line(data = df2, aes(x = x, y = mu), lwd=1)
g3 <- g3 + geom_line(data = df2, aes(x = x, y = lb), lwd=1, linetype = "dashed")
g3 <- g3 + geom_line(data = df2, aes(x = x, y = ub), lwd=1, linetype = "dashed")
g3 <- g3 + theme(legend.position = "none")
g3 <- g3 + xlab("x") + ylab("g")
g3

g4 <- ggplot(df,aes(x = x, y = y)) + geom_line(aes(color=factor(runId)))
df2 <- data.frame(x = tw[1,], mu = means[startpt:endpt], lb = means[startpt:endpt] - 1.96*sqrt(vars[startpt:endpt]), ub = means[startpt:endpt] + 1.96*sqrt(vars[startpt:endpt]), trueVals = trueOuts)
g4 <- g4 + geom_line(data = df2, aes(x = x, y = trueVals), lwd=1, color='blue')
g4 <- g4 + geom_line(data = df2, aes(x = x, y = mu), lwd=1)
g4 <- g4 + geom_line(data = df2, aes(x = x, y = lb), lwd=1, linetype = "dashed")
g4 <- g4 + geom_line(data = df2, aes(x = x, y = ub), lwd=1, linetype = "dashed")
g4 <- g4 + theme(legend.position = "none")
g4 <- g4 + xlab("x") + ylab("g")
g4
