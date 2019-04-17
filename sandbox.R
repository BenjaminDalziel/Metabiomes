# Sandbox script for developing functions




# Clean up and source project functions -----------------------------------

rm(list=ls())
graphics.off()

setwd('~/Dropbox/Research/Active/Metabiomes')

source('SampleCommunityMatrix.R')
source('DrawCommunityMatrix.R')



# Try things --------------------------------------------------------------

# Testing DrawCommunityMatrix
# quartz(h=6,w=6)
# par(pin = c(5,5))
# par(xpd = T)
# plot(0, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n', xlim = c(-1,1), ylim = c(-1,1))
# 
# A <- SampleCommunityMatrix(S = 50, C = 0.1, sigma = 1, Pm = 0.25, Pc = 0.25, s = -1)
# for(i in 1:10){
#   DrawCommunityMatrix(A = A, x = runif(1, 0, 1), y = runif(1, 0, 1), cex = 0.05, lod = 1)
# }
# 
# B <- SampleCommunityMatrix(S = 50, C = 0.1, sigma = 1, Pm = 0.25, Pc = 0.25, s = -1)
# for(i in 1:10){
#   DrawCommunityMatrix(A = B, x = runif(1, -1, 0), y = runif(1, -1, 0), cex = 0.05, lod = 1)
# }



# Testing an alternate version of DrawCommunityMatrix that uses an eigenvector of community matrix
# to deterine shape f a 'blob' representing the community in the metacommunity network
# so we can look and see same blob morphology = same ruleset (similar blob morph. = similar ruleset?)

# consider using the 'superformula' lol, and letting an eigensomething
# map to the parameters of the formula - number of axes of symmetry, etc.

A <- SampleCommunityMatrix(S = 50, C = 0.05, sigma = 1, Pm = 0.25, Pc = 0.25, s = 1)


x0 <- 0
y0 <- 0

r.base <- 0.005
r.gain <- 0.5
bw <- 5

eig <- eigen(A)
v <- abs(eig$vectors[,1])
n <- length(v)

vs <- as.numeric(filter(v, method = "conv", circular = TRUE, filter = rep(1/bw, bw)))

theta <- seq(0, 2*pi, length.out = n)
r <- r.base + r.gain * vs

x <- x0 + r * cos(theta)
y <- y0 + r * sin(theta)

x <- smooth.spline(x)$y
y <- smooth.spline(y)$y

x[n] <- x[1]
y[n] <- y[1]



quartz(h=6,w=6)
par(pin = c(5,5))
par(xpd = T)
plot(0, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n', xlim = c(-1,1), ylim = c(-1,1))
polygon(x,y)
