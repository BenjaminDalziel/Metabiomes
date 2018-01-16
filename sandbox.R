# Sandbox script for developing functions




# Clean up and source project functions -----------------------------------

rm(list=ls())
graphics.off()

setwd('~/Dropbox/Research/Metabiomes')

source('SampleCommunityMatrix.R')
source('DrawCommunityMatrix.R')



# Try things --------------------------------------------------------------

# Testing DrawCommunityMatrix
quartz(h=6,w=6)
par(pin = c(5,5))
par(xpd = T)
plot(0, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n', xlim = c(-1,1), ylim = c(-1,1))

A <- SampleCommunityMatrix(S = 50, C = 0.1, sigma = 1, Pm = 0.25, Pc = 0.25, s = -1)
A <- SampleCommunityMatrix(S = 50, C = 0.05, sigma = 1, Pm = 0.25, Pc = 0.25, s = -1)
for(i in 1:10){
  DrawCommunityMatrix(A = A, x = runif(1, 0, 1), y = runif(1, 0, 1), cex = 0.05, lod = 1)
}

B <- SampleCommunityMatrix(S = 50, C = 0.1, sigma = 1, Pm = 0.25, Pc = 0.25, s = -1)
B <- SampleCommunityMatrix(S = 50, C = 0.05, sigma = 1, Pm = 0.25, Pc = 0.25, s = -1)
for(i in 1:10){
  DrawCommunityMatrix(A = B, x = runif(1, -1, 0), y = runif(1, -1, 0), cex = 0.05, lod = 1)
}



