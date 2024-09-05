# Diagnostic script for testing functions

rm(list = ls())
graphics.off()

source("code/SampleCommunityMatrix.R")


# Choose which tests to run -----------------------------------------------

TestBasicMatrixProperties <- FALSE # Does SampleCommunityMatrix return matrices whose properties match input parameters?
TestMay1972 <- TRUE # Does SampleCommunityMatrix return matrices that reproduce main result of May 1972?
TestCoyte2015 <- FALSE # Does SampleCommunityMatrix return matrices that predict destabilizing effect of cooperation?



# BasicMatrixProperties ---------------------------------------------------
# Does SampleCommunityMatrix return matrices whose properties match input parameters?

if (TestBasicMatrixProperties) {
  n <- 50
  S <- 100

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
    matrix_para <- list(S = S, C = Cseq[i], sigma = 0.1, Pm = 0.2, Pc = 0.3, s = -1)
    A <- SampleCommunityMatrix(matrix_para)
    diag(A) <- NA
    Chat[i] <- sum(A != 0, na.rm = TRUE) / (S^2 - S)

    # Proportion cooperative
    matrix_para <- list(S = S, C = 0.2, sigma = 0.1, Pm = Pmseq[i], Pc = 0, s = -1)
    A <- SampleCommunityMatrix(matrix_para)
    diag(A) <- NA
    edges <- sum(A != 0, na.rm = TRUE)
    count <- 0
    for (j in 2:S) {
      for (k in 1:(j - 1)) {
        if (A[j, k] > 0 && A[k, j] > 0) {
          count <- count + 2
        }
      }
    }
    Pmhat[i] <- count / edges


    # Proportion competitive
    matrix_para <- list(S = S, C = 0.2, sigma = 0.1, Pm = 0, Pc = Pcseq[i], s = -1)
    A <- SampleCommunityMatrix(matrix_para)
    diag(A) <- NA
    edges <- sum(A != 0, na.rm = TRUE)
    count <- 0
    for (j in 2:S) {
      for (k in 1:(j - 1)) {
        if (A[j, k] < 0 && A[k, j] < 0) {
          count <- count + 2
        }
      }
    }
    Pchat[i] <- count / edges


    # Sigma
    matrix_para <- list(S = S, C = 0.2, sigma = sigmaseq[i], Pm = 0.2, Pc = 0.3, s = -1)
    A <- SampleCommunityMatrix(matrix_para)
    diag(A) <- 0
    sigmahat[i] <- sd(A[A != 0])
  }

  quartz(h = 6, w = 5)
  par(mfrow = c(2, 2))
  plot(Cseq, Chat, xlab = "Specified", ylab = "Observed", main = expression(paste(C)))
  abline(0, 1, col = 2)
  plot(Pmseq, Pmhat, xlab = "Specified", ylab = "Observed", main = expression(paste(P[m])))
  abline(0, 1, col = 2)
  plot(Pcseq, Pchat, xlab = "Specified", ylab = "Observed", main = expression(paste(P[c])))
  abline(0, 1, col = 2)
  plot(sigmaseq, sigmahat, xlab = "Specified", ylab = "Observed", main = expression(paste(sigma)))
  abline(0, 1, col = 2)
}









# May1972 -----------------------------------------------------------------
# Does SampleCommunityMatrix return matrices that reproduce main result of May 1972?


if (TestMay1972) {
  S <- 100
  Pm <- 0.25
  Pc <- 0.25
  s <- 1

  nrep <- 100
  lambda <- rep(NA, nrep)

  sigma <- runif(nrep, 0, 1)
  C <- runif(nrep, 0, 1)

  for (i in 1:nrep) {
    matrix_para <- list(S = S, C = C[i], sigma = sigma[i], Pm = Pm, Pc = Pc, s = s)
    A <- SampleCommunityMatrix(matrix_para)
    lambda[i] <- max(Re(eigen(A)$values)) # stable iff the eigenvalue with the largest real part is negative
    print(i)
  }

  n <- S
  alpha <- sigma^2
  quartz(h = 4, w = 4)
  par(pin = c(2, 2))
  plot(log((n * C)^-(1 / 2)), log(alpha), col = as.numeric(lambda > 1) + 1, ylab = expression(paste(log(alpha))), xlab = expression(paste(log((n * C)^-(1 / 2)))))
  abline(0, 1)
  legend("bottomright", col = c(2, 1), legend = c("unstable", "stable"), pch = 1, bty = "n")
}




# Coyte2015 ---------------------------------------------------------------
# Does SampleCommunityMatrix return matrices that predict destabilizing effect of cooperation?
# Reproducing Figure 2B right panel

if (TestCoyte2015) {
  S <- 100
  C <- 0.7
  s <- 1
  sigma <- 0.05

  n <- 50
  Pm <- seq(0, 1, length = n)
  Pc <- 1 - Pm

  lambda <- rep(NA, n)

  for (i in 1:n) {
    matrix_para <- list(S = S, C = C, sigma = sigma, Pm = Pm[i], Pc = Pc[i], s = s)
    A <- SampleCommunityMatrix(matrix_para)
    lambda[i] <- max(Re(eigen(A)$values)) # stable iff the eigenvalue with the largest real part is negative
  }

  quartz(h = 4, w = 4)
  plot(Pm, lambda, col = as.numeric(lambda > 1) + 1, xlab = expression(paste(P[m])), ylab = expression(paste(lambda)))
  legend("topleft", col = c(2, 1), legend = c("unstable", "stable"), pch = 1, bty = "n")
}
