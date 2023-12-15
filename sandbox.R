# Sandbox script for developing functions

rm(list = ls())
graphics.off()

source("SampleCommunityMatrix.R")
source("SampleCommunityMatrices.R")



# Sketch of sampling community matrices
M <- 1
S <- 100
Pm <- 0.6
Pc <- 0        #Pm + Pc must be less than or equal to 1
s <- 1
sigma <- 0.02
C <- 0.7

matrix_para <- list(
    S = S, Pm = Pm, Pc = Pc, s = s, sigma = sigma, C = C
)

AA <- SampleCommunityMatrices(M, matrix_para)





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


# Place migration terms along the off-block-diagonal
