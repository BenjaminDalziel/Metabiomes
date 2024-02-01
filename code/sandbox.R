# Sandbox script for developing functions
# Currently being used to sketch assembling the stability matrix

rm(list = ls())
graphics.off()

source("code/SampleCommunityMatrix.R")
source("code/SampleCommunityMatrices.R")



# Sample community matrices
L <- 3 # number of locations
S <- 100 # number of taxa
Pm <- 0.3  #proportion mutualist
Pc <- 0 # proportion competitive;  Pm + Pc must be less than or equal to 1
s <- 1
sigma <- 0.02
C <- 0.7
m <- 1 # migration rate; set to 0 if L = 1
r <- rep(1, L * S)

matrix_para <- list(
    S = S, Pm = Pm, Pc = Pc, s = s, sigma = sigma, C = C
)

AA <- SampleCommunityMatrices(L, matrix_para)



# Assemble the metacommunity interaction matrix A
A <- matrix(0, S * L, S * L)

subcommunity_mask <- function(a, b, S, L) {
    # Return a mask of the metacommunity matrix
    # corresponding to interactions between taxa
    # in subcommunities a and b

    rows <- ((a - 1) * S + 1):(a * S)
    cols <- ((b - 1) * S + 1):(b * S)

    mask <- matrix(FALSE, S * L, S * L)
    mask[rows, cols] <- TRUE

    return(mask)
}

# Place interaction matrices along the block diagonal
for (i in 1:L) {
    mask <- subcommunity_mask(i, i, S, L)
    A[mask] <- AA[, , i]
}


# Plot
mimage <- function(A, main = "") {
    image(A[, nrow(A):1], main = main)
}

mimage(A, main = "Metacommunity interaction matrix, A")




# Assemble metacommunity immigration matrix
M <- matrix(0, S * L, S * L)

row_taxa <- matrix(1:S, nrow = L * S, ncol = L * S, byrow = FALSE)
col_taxa <- t(row_taxa)

x <- as.numeric(gl(L, S, S * L))
row_subcom <- matrix(x, nrow = L * S, ncol = L * S, byrow = FALSE)
col_subcom <- t(row_subcom)


for (i in 1:nrow(M)) {
    for (j in 1:ncol(M)) {
        is_different_subcom <- row_subcom[i, j] != col_subcom[i, j]
        is_same_taxa <- row_taxa[i, j] == col_taxa[i, j]
        is_disperal <- is_different_subcom & is_same_taxa

        if (is_disperal) {
            M[i, j] <- m / (L - 1)
        }
    }
}

# Plot
mimage(M, main = "Metacommunity immigration matrix, M")



# Find feasibile equilibrium N^*
Ainv <- solve(A)
Nstar_nomigration <- -Ainv %*% r

g <- function(N) {

    -m + 1 / N * M %*% N

}

Nhatstar <- -Ainv %*% (g(Nstar_nomigration) + r)

par(cex = 2)
plot(Nstar_nomigration, Nhatstar)
abline(0, 1)


# TODO: pass Nhatstar to fsolve to try to get closer to Nstar
library(pracma)