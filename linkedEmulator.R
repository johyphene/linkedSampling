fEmulator <- function(fLocs, fOutput, fLocsPred, kernel = "mat52") {
  m <- nrow(fLocs)
  trendMat <- cbind(fLocs, matrix(1,m,1)) #the trend always be this way for f

  if (is.integer(dim(fLocsPred)) == TRUE)
  {
    testTrendMat <- cbind(fLocsPred, 1) 
  } else {
    testTrendMat <- cbind(t(as.matrix(fLocsPred)), 1)
  }
  
  if (kernel == "exp") { #exponential currently does not work with the RobustGaSP package
    fModel <- ppgasp(design=fLocs, response=fOutput, trend=trendMat, kernel_type = 'pow_exp', alpha = 1.0) #ppgasp for f layer
    fPred <- predict.ppgasp(fModel, fLocsPred, testing_trend = testTrendMat, kernel_type = 'pow_exp', alpha = 1.0)
  }
  else if (kernel == "sq_exp" || kernel == "pow_exp") {
    fModel <- ppgasp(design=fLocs, response=fOutput, trend=trendMat, kernel_type = 'pow_exp') #ppgasp for f layer
    fPred <- predict.ppgasp(fModel, fLocsPred, testing_trend = testTrendMat, kernel_type = 'pow_exp')
  }
  else if (kernel == "mat32") {
    fModel <- ppgasp(design=fLocs, response=fOutput, trend=trendMat, kernel_type = 'matern_3_2') #ppgasp for f layer
    fPred <- predict.ppgasp(fModel, fLocsPred, testing_trend = testTrendMat, kernel_type = 'matern_3_2')
  }
  else {
    fModel <- ppgasp(design=fLocs, response=as.matrix(fOutput), trend=trendMat) #ppgasp for f layer
    fPred <- predict.ppgasp(fModel, fLocsPred, testing_trend = testTrendMat)
    }
  
  outList <- list("model" = fModel, "pred" = fPred)
  return(outList) 
}
###
gEmulator <- function(gLocs, gOutput, gLocsPred, trend, kernel = "mat52") {
  m <- nrow(gLocs)
  numTest <- nrow(gLocsPred)
  if(trend == 'constant') {
    trendMat <- matrix(1,m,1)
    testTrendMat <- matrix(1,numTest,1)
  }
  if(trend == 'linear') {
    trendMat <- cbind(gLocs, matrix(1,m,1))
    testTrendMat <- cbind(gLocsPred, matrix(1,numTest,1))
  }
  
  if (kernel == "exp") {
    gModel<- ppgasp(design=gLocs, response=gOutput, trend=trendMat, kernel_type = 'pow_exp', alpha = 1.0)
    gPred <- predict.ppgasp(gModel, gLocsPred, testing_trend = testTrendMat, kernel_type = 'pow_exp', alpha = 1.0)
  }
  else if (kernel == "sq_exp" || kernel == "pow_exp") {
    gModel<- ppgasp(design=gLocs, response=gOutput, trend=trendMat, kernel_type = 'pow_exp')
    gPred <- predict.ppgasp(gModel, gLocsPred, testing_trend = testTrendMat, kernel_type = 'pow_exp')
  }
  else if (kernel == "mat32") {
    gModel<- ppgasp(design=gLocs, response=gOutput, trend=trendMat, kernel_type = 'matern_3_2')
    gPred <- predict.ppgasp(gModel, gLocsPred, testing_trend = testTrendMat, kernel_type = 'matern_3_2')
  }
  else {
    gModel<- ppgasp(design=gLocs, response=gOutput, trend=trendMat)
    gPred <- predict.ppgasp(gModel, gLocsPred, testing_trend = testTrendMat)
  }
  outList <- list("model" = gModel, "pred" = gPred)
  return(outList) 
}
###
linkedSampling <- function(fLocsPred, gLocsPred, muLs, sigma2Ls, fGamma, gGamma, N, gPred){
  
  A=matrix(1, nrow = N, ncol = N)

  numFCols = ncol(fLocsPred)
  if (is.null(numFCols) == TRUE) {
    numFCols = length(fLocsPred)
    fLocsPred = matrix(fLocsPred, nrow = 1)
  }
  numGCols = ncol(gLocsPred)
  
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        A[i,j] = sigma2Ls[i]
      }
      else if (i != j) {
        rho_f = 1
        for (k in 1:numFCols)  {
          if (nrow(fLocsPred) > 1) {
            rho_f = rho_f * cKSingle(fLocsPred[i,k],fLocsPred[j,k],fGamma[k]) 
            } else {
            rho_f = rho_f * cKSingle(fLocsPred[k],fLocsPred[k],fGamma[k]) 
            }
        }
        tempSig = sqrt(sigma2Ls[i] - 2*rho_f*sqrt(sigma2Ls[i])*sqrt(sigma2Ls[j]) + sigma2Ls[j])
        rho_eta = 1
        for (k in 1:numFCols)  {
          rho_eta = rho_eta * xiSingle(muLs[i], tempSig, gGamma[k], muLs[j]) 
          }

        A[i,j] = sqrt(sigma2Ls[i]) * sqrt(sigma2Ls[j]) * rho_eta * rho_f
      }
    }
  }
  L = chol(A + eye(N)*0.0000000001) #regularizing the matrix so it is positive definite
  L = t(L)
  numTest <- dim(gLocsPred)[1]
  samples <- matrix(0, nrow = numTest, ncol = 100)
  for (i in 1:100) {
    u = rnorm(N)
    ysamp = gPred + L%*%as.matrix(u)
    samples[,i] <- ysamp
  }
  outList <- list("A" = A, "samples" = samples) 
  return(outList)
  
}
###
linkedEmulator <- function(fLocs, fOutput, fLocsPred, zLocs, gOutput, zLocsPred, trend, kernel = "mat52") {
  #load in correct .cpp file for kernel specified
  library(Rcpp)
  setwd("C:/") #this will be different for package I assume
  if (kernel == "exp") {
    sourceCpp('linkedEmulatorExponential.cpp') # file with C++ versions of the linked emulator functions
  } else if (kernel == "sq_exp" || kernel == "pow_exp") {
    sourceCpp('linkedEmulatorSquaredExponential.cpp') # file with C++ versions of the linked emulator functions
  } else if (kernel == "mat32") {
    sourceCpp('linkedEmulatorMatern32.cpp') # file with C++ versions of the linked emulator functions
  } else {
    sourceCpp('linkedEmulator.cpp') # file with C++ versions of the linked emulator functions
  }
  
  set.seed(9)
  N <- nrow(as.matrix(fLocsPred))
  numFInputs <- ncol(fLocs)
  numFOutputs <- ncol(fOutput)
  ygpVals <- matrix(0, nrow = N, ncol = numFOutputs)
  sdVals <- matrix(0, nrow = N, ncol = numFOutputs)
  gammas <- matrix(0, nrow = numFOutputs, ncol = numFInputs)
  fGamma <- rep(1, times=numFInputs)
  #first layer PPE - take in f^D inputs/outputs, take in f^* inputs (should not need outputs?)
  for (jj in 1:numFOutputs) {
    y <- fOutput[,jj]
    fVals = fEmulator(fLocs, y, fLocsPred, kernel)
    
    ygp = fVals$pred$mean
    v = fVals$pred$sd
    ygpVals[,jj] <- ygp
    sdVals[,jj] <- v
    gammas[jj,] <- 1/fVals$model@beta_hat
    for (jjj in 1:numFInputs) {
      fGamma[jjj] <- fGamma[jjj]*gammas[jj,jjj]
    }
  }

  fModel <- fVals$model
  fPred <- list(mean = ygpVals, sd = sdVals)
  
  if (is.null(dim(fOutput)) == TRUE) {
    gLocs <- cbind(matrix(t(fOutput), ncol = 1), zLocs)
  } else {
    gLocs <- cbind(fOutput, zLocs)
  }
  
  if (dim(as.matrix(fPred$mean))[1] == 1) {
    gLocsPred <- cbind(t(fPred$mean), zLocsPred)
  } else {
    gLocsPred <- cbind(fPred$mean, zLocsPred)
  }
  #g layer - all inputs and outputs are imported - might want to consider how that will look
  #train inputs are called gLocs, train outputs are called gOutput, test inputs are called fLocsPred
  gVals = gEmulator(gLocs, gOutput, gLocsPred, trend, kernel)
  gModel <- gVals$model
  gPred <- gVals$pred

  #based on trend call one of two functions for the linked emulator
  muVec = as.matrix(fPred$mean)

  sigVec = as.matrix(fPred$sd)

  d <- ncol(muVec)
  p <- ncol(gLocs) - d
  in_all <- gLocs
  wT <- gLocs[,1:d]
  zT <- gLocs[,(d+1):(d+p)]
  m <- dim(in_all)[1]
  eta <- gModel@nugget #nugget term
  gSigVec <- gModel@sigma2_hat
  gGamma <- 1/gModel@beta_hat

  c_k <- array(0, dim = c(m, m, (d+p)))
  for (i in 1:m) {
    for (j in 1:m) {
      for (k in 1:(d+p)) {
        c_k[i,j,k] <- cKSingle(in_all[i,k], in_all[j,k], gGamma[k])
      }
    }
  }
  for (k in 1:(d+p)) {
    if (k == 1)
      R <- c_k[,,k]    
    else 
      R <- R * c_k[,,k]
  }
  if (eta == 0) {
    eta = 0.000000000001
  }
  invR <- solve(R + diag(eta, nrow(R), ncol(R))) #regularize R matrix
  g = ncol(as.matrix(gOutput)) #number of outputs we are interested in
  print(paste0("g = ", g))
  numTest <- nrow(gLocsPred)
  if (g == 1) {
    yT <- gOutput
    if (trend == "constant") {
      #call 1D constant trend
      HzT <- matrix(1,m,1)
      hTilde <- matrix(1,m,1)
      hats <- solve(t(hTilde) %*% invR %*% hTilde)%*%t(hTilde) %*% invR %*% yT
      thetaHat <- rep(0, times = ncol(as.matrix(wT)))
      thetaHat <- as.matrix(thetaHat)
      betaHat <- as.matrix(hats)
      sigInvMat <- solve(t(hTilde)%*%invR%*%hTilde)
      B <- matrix(0,nrow = d, ncol = m)
      Q <- invR %*% hTilde %*% solve(t(hTilde) %*% invR %*% hTilde) %*% t(hTilde) %*% invR - invR
      C <- solve(t(hTilde) %*% invR %*% hTilde)
      muLs <- matrix(c(0), nrow = numTest, ncol = g)
      sigma2Ls <- matrix(c(0), nrow = numTest, ncol = g)
      for (i in 1:g) {
        A <- invR%*%(yT - wT %*% thetaHat - HzT %*% betaHat)
        for(k in 1:numTest){
          print(k)  
          Sys.sleep(0.001)  
          flush.console()  
          hZ <- 1
          
          iVec <- iFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          B <- bFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          K <- cbind(t(B), iVec%*%t(hZ))
          J <- jFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          tempOmega <- sigVec[k,]^2
          Omega <- matrix(diag(tempOmega), ncol = length(tempOmega))
          P <- blkdiag(Omega, matrix(0,nrow = length(hZ), ncol = length(hZ)))
          G <- cbind(t(muVec[k,]), t(hZ))
          G <- t(G)
          trQJ = 0
          for (bb in 1:ncol(Q)) {
            trQJ = trQJ + (Q[bb,] %*% J[,bb])
          }
          mu_L <- (t(muVec[k,]) %*% thetaHat) + (t(hZ) %*% betaHat) + (t(iVec) %*% A)
          sigma2_L <- t(A) %*% (J - iVec%*%t(iVec)) %*% A + gSigVec[i]*(1 + eta + sigInvMat + trQJ - 2*tr(C * t(hTilde) %*% invR %*% iVec))
          muLs[k, i] <- mu_L
          sigma2Ls[k, i] <- sigma2_L
          if(sigma2_L < 0) { 
            if(sigma2_L < -(10^-7)) {
              if(sigma2_L < -(10^-7)) {
                print(paste0("Negative variance at index ", k))
                print(paste0("Value of ", sigma2_L))
              }
            }
            sigma2Ls[k,i] = abs(sigma2_L)
          }
        }
      }
      samps = linkedSampling(fLocsPred, gLocsPred, muLs, sigma2Ls, fGamma, gGamma, numTest, gPred$mean)
    }
    if (trend == "linear") {
      #call 1D linear trend
      HzT <- cbind(zT,matrix(1,m,1))
      hTilde <- cbind(wT, HzT)
      hats <- solve(t(hTilde) %*% invR %*% hTilde)%*%t(hTilde) %*% invR %*% yT
      thetaHat <- as.matrix(hats[1:d,])
      betaHat <- as.matrix(hats[(d+1):(d+p+1),]) 
      sigInvMat <- solve(t(hTilde)%*%invR%*%hTilde)
      B <- matrix(0,nrow = d, ncol = m)
      Q <- invR %*% hTilde %*% solve(t(hTilde) %*% invR %*% hTilde) %*% t(hTilde) %*% invR - invR
      C <- solve(t(hTilde) %*% invR %*% hTilde)
      muLs <- matrix(c(0), nrow = numTest, ncol = g)
      sigma2Ls <- matrix(c(0), nrow = numTest, ncol = g)
      for (i in 1:g) {
        A <- invR%*%(yT - wT %*% thetaHat - HzT %*% betaHat)
        for(k in 1:numTest){
          print(k)  
          Sys.sleep(0.001)  
          flush.console()  
          hZ <- c(gLocsPred[k,(d+1):(d+p)],1)
          
          iVec <- iFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          B <- bFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          K <- cbind(t(B), iVec%*%t(hZ))
          J <- jFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          tempOmega <- sigVec[k,]^2
          if (length(tempOmega) == 1) {
            Omega <- as.matrix(tempOmega)
          } else { 
            Omega <- matrix(diag(tempOmega), ncol = length(tempOmega))
          }
          P <- blkdiag(Omega, matrix(0,nrow = length(hZ), ncol = length(hZ)))
          G <- cbind(t(muVec[k,]), t(hZ))
          G <- t(G)
          trQJ = 0
          for (bb in 1:ncol(Q)) {
            trQJ = trQJ + (Q[bb,] %*% J[,bb])
          }
          mu_L <- (t(muVec[k,]) %*% thetaHat) + (t(hZ) %*% betaHat) + (t(iVec) %*% A)
          sigma2_L <- t(A) %*% (J - iVec%*%t(iVec)) %*% A + 2*t(thetaHat)%*%(B - muVec[k,] %*% t(iVec)) %*% A + tr(thetaHat %*% t(thetaHat) %*% Omega) + gSigVec[i] *(1 + eta + trQJ + t(G) %*% C %*% G + tr(C %*% P - 2*C %*% t(hTilde) %*% invR %*% K))
          muLs[k, i] <- mu_L
          sigma2Ls[k, i] <- sigma2_L
          if(sigma2_L < 0) { 
            if(sigma2_L < -(10^-7)) {
              if(sigma2_L < -(10^-7)) {
                print(paste0("Negative variance at index ", k))
                print(paste0("Value of ", sigma2_L))
              }
            }
            sigma2Ls[k,i] = abs(sigma2_L)
          }
        }
      }
      samps = linkedSampling(fLocsPred, gLocsPred, muLs, sigma2Ls, fGamma, gGamma, numTest, gPred$mean)
    }
  }
  if (g > 1) {
    if (trend == "constant") {
      #call 2D constant trend
      muLs <- matrix(c(0), nrow = numTest, ncol = g)
      sigma2Ls <- matrix(c(0), nrow = numTest, ncol = g)
      sampsA <- array(0, dim = c(numTest, numTest, g))
      sampsSamples <- array(0, dim = c(numTest, 100, g))
      for (i in 1:g) {
        yT <- gOutput[,i]
        HzT <- matrix(1,m,1)
        hTilde <- matrix(1,m,1)
        hats <- solve(t(hTilde) %*% invR %*% hTilde)%*%t(hTilde) %*% invR %*% yT
        thetaHat <- rep(0, times = ncol(fLocs))
        betaHat <- hats
        
        sigInvMat <- solve(t(hTilde)%*%invR%*%hTilde) #saving computation time since this doesn't change each iteration
        B <- matrix(0,nrow = d, ncol = m)
        Q <- invR %*% hTilde %*% solve(t(hTilde) %*% invR %*% hTilde) %*% t(hTilde) %*% invR - invR
        C <- solve(t(hTilde) %*% invR %*% hTilde)
        A <- invR%*%(yT - wT %*% thetaHat - HzT %*% betaHat)
        
        for(k in 1:numTest){
          print(k)  
          Sys.sleep(0.001)  
          flush.console()  
          hZ <- 1
          
          iVec <- iFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          B <- bFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          K <- cbind(t(B), iVec%*%t(hZ))
          J <- jFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          tempOmega <- sigVec[k,]^2
          Omega <- matrix(diag(tempOmega), ncol = length(tempOmega))
          P <- blkdiag(Omega, matrix(0,nrow = length(hZ), ncol = length(hZ)))
          G <- cbind(t(muVec[k,]), t(hZ))
          G <- t(G)
          trQJ = 0
          for (bb in 1:ncol(Q)) {
            trQJ = trQJ + (Q[bb,] %*% J[,bb])
          }
          mu_L <- (t(muVec[k,]) %*% thetaHat) + (t(hZ) %*% betaHat) + (t(iVec) %*% A)
          sigma2_L <- t(A) %*% (J - iVec%*%t(iVec)) %*% A + gSigVec[i]*(1 + eta + sigInvMat + trQJ - 2*tr(C * t(hTilde) %*% invR %*% iVec))
          muLs[k, i] <- mu_L
          sigma2Ls[k, i] <- sigma2_L
          if(sigma2_L < 0) {
            print(paste0("Negative variance at index ", k))
            print(paste0("Value of ", sigma2_L))
            sigma2Ls[k,i] = abs(sigma2_L)
          }
        }
        tempSamps = linkedSampling(fLocsPred, gLocsPred, muLs[,i], sigma2Ls[,i], fGamma, gGamma, numTest, gPred$mean[,i])
        sampsA[,,i] = tempSamps$A
        sampsSamples[,,i] = tempSamps$samples
      }
      samps = list(A = sampsA, samples = sampsSamples)
    }
    if (trend == "linear") {
      #call 2D linear trend
      muLs <- matrix(c(0), nrow = numTest, ncol = g)
      sigma2Ls <- matrix(c(0), nrow = numTest, ncol = g)
      sampsA <- array(0, dim = c(numTest, numTest, g))
      sampsSamples <- array(0, dim = c(numTest, 100, g))
      for (i in 1:g) {
        yT <- gOutput[,i]
        HzT <- cbind(zT,matrix(1,m,1))
        hTilde <- cbind(wT, HzT)
        hats <- solve(t(hTilde) %*% invR %*% hTilde)%*%t(hTilde) %*% invR %*% yT
        thetaHat <- hats[1:d,]
        betaHat <- hats[(d+1):(d+p+1),]
        
        sigInvMat <- solve(t(hTilde)%*%invR%*%hTilde) #saving computation time since this doesn't change each iteration
        B <- matrix(0,nrow = d, ncol = m)
        Q <- invR %*% hTilde %*% solve(t(hTilde) %*% invR %*% hTilde) %*% t(hTilde) %*% invR - invR
        C <- solve(t(hTilde) %*% invR %*% hTilde)
        A <- invR%*%(yT - wT %*% thetaHat - HzT %*% betaHat)
        
        for(k in 1:numTest){
          print(k)
          Sys.sleep(0.001) 
          flush.console()
          hZ <- c(gLocsPred[k,(d+1):(d+p)],1)
          
          iVec <- iFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          B <- bFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          K <- cbind(t(B), iVec%*%t(hZ))
          J <- jFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
          tempOmega <- sigVec[k,]^2
          Omega <- matrix(diag(tempOmega), ncol = length(tempOmega))
          P <- blkdiag(Omega, matrix(0,nrow = length(hZ), ncol = length(hZ)))
          G <- cbind(t(muVec[k,]), t(hZ))
          G <- t(G)
          trQJ = 0
          for (bb in 1:ncol(Q)) {
            trQJ = trQJ + (Q[bb,] %*% J[,bb])
          }
          mu_L <- (t(muVec[k,]) %*% thetaHat) + (t(hZ) %*% betaHat) + (t(iVec) %*% A)
          sigma2_L <- sigma2_L <- t(A) %*% (J - iVec%*%t(iVec)) %*% A + 2*t(thetaHat)%*%(B - muVec[k,] %*% t(iVec)) %*% A + tr(thetaHat %*% t(thetaHat) %*% Omega) + gSigVec[i] *(1 + eta + trQJ + t(G) %*% C %*% G + tr(C %*% P - 2*C %*% t(hTilde) %*% invR %*% K))
          muLs[k, i] <- mu_L
          sigma2Ls[k, i] <- sigma2_L
          if(sigma2_L < 0) { #add some sort of flag here
            print(paste0("Negative variance at index ", k))
            print(paste0("Value of ", sigma2_L))
            sigma2Ls[k,i] = abs(sigma2_L)
          }
        }
        tempSamps = linkedSampling(fLocsPred, gLocsPred, muLs[,i], sigma2Ls[,i], fGamma, gGamma, numTest, gPred$mean[,i])
        sampsA[,,i] = tempSamps$A
        sampsSamples[,,i] = tempSamps$samples
      }
      samps = list(A = sampsA, samples = sampsSamples)
    }
  }
  #return muLs, sigma2Ls, samps
  outList <- list("means" = gPred$mean, "vars" = samps$A, "samples" = samps$samples, "gM" = gModel)
  return(outList)
}

f_1 <- function(x1, x2) {
  return (sin(x1) + x2^2)
}

f_2 <- function(x1, x2) {
  return (sin(3*x1*cos(pi*x2)))
}

g_1 <- function(w1, w2, z1, z2) {
  return (cos(3*z2)*sin(3*w1) + sin(3*z1)*cos(3*w2))
}

g_2 <- function(w1, w2, z1, z2) {
  return (cos(4*z2)*sin(2*w1) + sin(4*z1)*cos(2*w2))
}

#load other necessary packages (RobustGaSP, readxl, etc.)
library(lhs)
library(RobustGaSP)
library(readxl) #use if you'd like to import data through excel rather than generating it
library(R.matlab) #use if you have data from MATLAB

#generate data
set.seed(9)
N = 201
m = 300
numX = 2
x1=randomLHS(m,1)
x1 = x1*2
x2=randomLHS(m,1)
x2 = x2*2

w1=f_1(x1, x2)
w2=f_2(x1, x2)
fLocs <- cbind(x1, x2)
trendMat <- cbind(fLocs,matrix(1,m,1))

xx1 = seq(0,2, length.out = N)
xx2 <- rep(1, times = N)
xx2 = xx2 + (runif(N)*(10^-7))
fLocsPred <- cbind(xx1, xx2)
testTrendMat <- cbind(fLocsPred,matrix(1,N,1))

ygpVals <- randomLHS(N,2)
ygpVals <- ygpVals * 0
sdVals <- randomLHS(N,2)
sdVals <- sdVals * 0
fOutput <- cbind(w1,w2)

ndp=dim(fLocs)[1]
numX <- dim(fLocs)[2]


w <- fOutput
z1 <- randomLHS(m,1)
z1 <- z1 * 2
z2 <- randomLHS(m,1)
z2 <- z2 * 2
z <- cbind(z1, z2)
zLocs <- z
gOutput <- cbind(g_1(w[,1], w[,2], z[,1], z[,2]), g_2(w[,1], w[,2], z[,1], z[,2]))

zT1 <- rep(1, times = N)
zT1 = zT1 + (runif(N)*(10^-7))
zT2 <- rep(1, times = N)
zT2 = zT2 + (runif(N)*(10^-7))
zTest <- cbind(zT1, zT2)
zz <- zTest
zLocsPred <- zTest

ww1=f_1(xx1, xx2)
ww2=f_2(xx1, xx2)
ww <- cbind(ww1,ww2)
gOutputTest <- cbind(g_1(ww[,1], ww[,2], zz[,1], zz[,2]), g_2(ww[,1], ww[,2], zz[,1], zz[,2]))

#specify trend
trend = 'linear'
kernel = "mat52"
#call linked emulator
lE <- linkedEmulator(fLocs, fOutput, fLocsPred, zLocs, gOutput, zLocsPred, trend, kernel)

###Plotting###
sigma2s <- diag(lE$vars[,,1])
library(ggplot2)
plotXs <- matrix(rep(fLocsPred[,1], times=100), nrow = N)
plotRuns <- t(matrix(rep(1:100, times=N), ncol = N))
df <- data.frame(x = as.vector(plotXs), y = as.vector(lE$samples[,,1]), runId = as.vector(plotRuns))
g2 <- ggplot(df,aes(x = x, y = y)) + geom_line(aes(color=factor(runId)))
df2 <- data.frame(x = fLocsPred[,1], mu = lE$means[,1], lb = lE$means[,1] - 1.96*sqrt(sigma2s), ub = lE$means[,1] + 1.96*sqrt(sigma2s), trueVals = gOutputTest[,1])
g2 <- g2 + geom_line(data = df2, aes(x = x, y = mu), lwd=1)
g2 <- g2 + geom_line(data = df2, aes(x = x, y = lb), lwd=1, linetype = "dashed")
g2 <- g2 + geom_line(data = df2, aes(x = x, y = ub), lwd=1, linetype = "dashed")
g2 <- g2 + geom_line(data = df2, aes(x = x, y = trueVals), lwd=1, color = "blue")
g2 <- g2 + theme(legend.position = "none")
g2

sigma2s <- diag(lE$vars[,,2])
df <- data.frame(x = as.vector(plotXs), y = as.vector(lE$samples[,,2]), runId = as.vector(plotRuns))
g3 <- ggplot(df,aes(x = x, y = y)) + geom_line(aes(color=factor(runId)))
df2 <- data.frame(x = fLocsPred[,1], mu = lE$means[,2], lb = lE$means[,2] - 1.96*sqrt(sigma2s), ub = lE$means[,2] + 1.96*sqrt(sigma2s), trueVals = gOutputTest[,2])
g3 <- g3 + geom_line(data = df2, aes(x = x, y = mu), lwd=1)
g3 <- g3 + geom_line(data = df2, aes(x = x, y = lb), lwd=1, linetype = "dashed")
g3 <- g3 + geom_line(data = df2, aes(x = x, y = ub), lwd=1, linetype = "dashed")
g3 <- g3 + geom_line(data = df2, aes(x = x, y = trueVals), lwd=1, color = "blue")
g3 <- g3 + theme(legend.position = "none")
g3

