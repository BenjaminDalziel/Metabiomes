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
sigma <- 0.05
C <- 0.7

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