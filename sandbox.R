# Sandbox script for developing functions

rm(list = ls())
graphics.off()

source("SampleCommunityMatrix.R")
source("SampleCommunityMatrices.R")


# Sketch of sampling community matrices
M <- 1
S <- 100
Pm <- 0.25
Pc <- 0.25
s <- 1
sigma <- 0.5
C <- 0.01

AA <- SampleCommunityMatrices(
    M = M, S = S, Pm = Pm, Pc = Pc,
    s = s, sigma = sigma, C = C
)



# Sketch of assembling metacommunity matrix
B <- matrix(NA, S * M, S * M)

subcommunity_mask <- function(a, b, S, M) {
    # Return a mask of the metacommunity matrix
    # corresponding to interactions between taxa
    # in subcommunities a and b

    rows <- ((a - 1) * S + 1):(a * S)
    cols <- ((b - 1) * S + 1):(b * S)

    mask <- matrix(FALSE, S * M, S * M)
    mask[rows, cols] <- TRUE

    return(mask)
}

# Place interaction matrices along the block diagonal
for (i in 1:M) {
    mask <- subcommunity_mask(i, i, S, M)
    B[mask] <- AA[, , i]
}


# TODO: Place migration terms along the off-block-diagonal



# Sketch of feasibility and stability search
FindEquilibrium <- function(A, r = rep(1, nrow(A))) {

    Nstar <- -solve(A) %*% r
    return(Nstar)

}

nrep <- 10000
for (i in 1:nrep) {

    print(i)
    A <- SampleCommunityMatrix(S, C, sigma, Pm, Pc, s)

    Nstar <- FindEquilibrium(A)

    D <- matrix(0, nrow(A), nrow(A))
    diag(D) <- Nstar

    lambda <- eigen(A)$values

    is_feasible <- all(Nstar > 0)
    is_stable <- all(Re(lambda) < 0)


    if (is_feasible) {
        print(paste(i, "is feasible!"))

        if (is_stable) {
            print("And it's stable!")
        } else {
            print("And it's unstable!")
        }

        break
    }
}








# SCRAPS

if (0) {
    # Try things --------------------------------------------------------------

    # Testing DrawCommunityMatrix
    source("DrawCommunityMatrix.R")


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
    v <- abs(eig$vectors[, 1])
    n <- length(v)

    vs <- as.numeric(filter(v, method = "conv", circular = TRUE, filter = rep(1 / bw, bw)))

    theta <- seq(0, 2 * pi, length.out = n)
    r <- r.base + r.gain * vs

    x <- x0 + r * cos(theta)
    y <- y0 + r * sin(theta)

    x <- smooth.spline(x)$y
    y <- smooth.spline(y)$y

    x[n] <- x[1]
    y[n] <- y[1]



    quartz(h = 6, w = 6)
    par(pin = c(5, 5))
    par(xpd = T)
    plot(0, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", xlim = c(-1, 1), ylim = c(-1, 1))
    polygon(x, y)
}
