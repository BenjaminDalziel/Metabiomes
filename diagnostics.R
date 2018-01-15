# Diagnostic script for testing functions

rm(list=ls())
graphics.off()

setwd('~/Dropbox/Research/Metabiomes')
source('SampleCommunityMatrix.R')
source('DrawCommunityMatrix.R')

# Does SampleCommunityMatrix produce matrices that match input para? ------
n <- 30
S <- 50

Cseq <- seq(0, 1, length.out = n)
Pmseq <- seq(0, 1, length.out = n)
Pcseq <- seq(0, 1, length.out = n)
sigmaseq <- seq(0, 1, length.out = n)

Chat <- rep(NA, n)
Pmhat <- rep(NA, n)
Pchat <- rep(NA, n)
sigmahat <- rep(NA, n)


for (i in 1:n) {
  
  # Connectivity
  A <- SampleCommunityMatrix(S = S, sigma = 0.1, C = Cseq[i], Pm = 0.2, Pc = 0.3, s = -1)
  diag(A) <- NA
  Chat[i] <- sum(A != 0, na.rm=T) / (S^2 - S)
  
  # Proportion cooperative
  A <- SampleCommunityMatrix(S = S, sigma = 0.1, C = 0.2, Pm = Pmseq[i], Pc = 0, s = -1)
  diag(A) <- NA
  edges <- sum(A != 0, na.rm = T)
  count <- 0
  for (j in 2:S){
    for (k in 1:(j-1)){
      if (A[j,k] > 0 & A[k,j] > 0) {
        count <- count + 2
      }
    }
  }
  Pmhat[i] <- count/edges
  
  
  # Proportion competitive
  A <- SampleCommunityMatrix(S = S, sigma = 0.1, C = 0.2, Pm = 0, Pc = Pcseq[i], s = -1)
  diag(A) <- NA
  edges <- sum(A != 0, na.rm = T)
  count <- 0
  for (j in 2:S){
    for (k in 1:(j-1)){
      if (A[j,k] < 0 & A[k,j] < 0) {
        count <- count + 2
      }
    }
  }
  Pchat[i] <- count/edges
  
  
  # Sigma
  A <- SampleCommunityMatrix(S = S, sigma = sigmaseq[i], C = 0.2, Pm = 0.2, Pc = 0.3, s = -1)
  diag(A) <- 0
  sigmahat[i] <- sd(A[A != 0])
  
}

quartz(h = 6, w = 5)
par(mfrow = c(2,2))
plot(Cseq, Chat, xlab = "Specified", ylab = "Observed"); abline(0,1,col=2)
plot(Pmseq, Pmhat, xlab = "Specified", ylab = "Observed"); abline(0,1,col=2)
plot(Pcseq, Pchat, xlab = "Specified", ylab = "Observed"); abline(0,1,col=2)
plot(sigmaseq, sigmahat, xlab = "Specified", ylab = "Observed"); abline(0,1,col=2)




