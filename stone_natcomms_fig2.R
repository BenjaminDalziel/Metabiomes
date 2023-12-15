# Sandbox script for developing functions

rm(list = ls())
graphics.off()

source("SampleCommunityMatrix.R")
source("SampleCommunityMatrices.R")



# Parameters for community matrices
M <- 1
S <- 100
Pm <- 0.6
Pc <- 0        #Pm + Pc must be less than or equal to 1
s <- 1
sigma <- 0.02
C <- 0.7


# Write parameters to a list
matrix_para <- list(
    S = S, Pm = Pm, Pc = Pc, s = s, sigma = sigma, C = C
)




# Sketch of feasibility and stability analysis

# find equilibrium abundances
FindEquilibrium <- function(A, r = rep(1, nrow(A))) {
    Nstar <- -solve(A) %*% r
    return(Nstar)
}


# function to sample a community matrix and
# determine its feasibility and stability
# according to the method of Stone or May
AnalyzeMatrix <- function(matrix_para = matrix_para, method = "May") {

    A <- SampleCommunityMatrix(matrix_para)

    # feasibility
    Nstar <- FindEquilibrium(A)
    is_feasible <- all(Nstar > 0)

    # stability
    if (method == "May") {
        eig <- eigen(A)
        is_stable <- all(Re(eig$values) < 0)
    }

    if (method == "Stone") {
        D <- matrix(0, nrow(A), nrow(A))
        diag(D) <- Nstar

        Z <- D %*% A # "stability matrix"
        eig <- eigen(Z)
        is_stable <- all(Re(eig$values) < 0)
    }

    return(list(is_stable = is_stable, 
                is_feasible = is_feasible,
                eig = eig$values))
}

# function to run replicate simulations and tabulate results
TabulateSimulations <- function(matrix_para, method, nrep){

    tab <- matrix(0, 2, 2)

    crit <- rep(NA, nrep)
    is_feasible <- rep(NA, nrep)
    is_stable <- rep(NA, nrep)

    for (i in 1:nrep) {
        res <- AnalyzeMatrix(matrix_para, method)

        is_feasible[i] <- res$is_feasible
        is_stable[i] <- res$is_stable

        if (res$is_feasible) {
            if (res$is_stable) {
                tab[1, 1] <- tab[1, 1] + 1
            } else {
                tab[1, 2] <- tab[1, 2] + 1
            }
        } else {
            if (res$is_stable) {
                tab[2, 1] <- tab[2, 1] + 1
            } else {
                tab[2, 2] <- tab[2, 2] + 1
            }
        }

        crit[i] <- max(Re(res$eig))  # "critical stability eigenvalue"

        #print(nrep - i + 1)
    }

    rownames(tab) <- c("Feasible", "Infeasible")
    colnames(tab) <- c("Stable", "Unstable")

    return(list(tab = tab,
                is_feasible = is_feasible, 
                is_stable = is_stable, 
                crit = crit))
}



# Reproduce Stone Nature Comms fig 2

Pm_seq <- seq(0, 0.6, 0.01)
nlev <- length(Pm_seq)

crit <- rep(NA, nlev)

for(i in 1:nlev){

    matrix_para <- list(
        S = S, Pm = Pm_seq[i], Pc = Pc, s = s, sigma = sigma, C = C
    )

    out <- TabulateSimulations(matrix_para, "Stone", 50)
    crit[i] <- mean(out$crit)

    print(i)

}

par(pin = c(4, 4))
plot(Pm_seq, crit, cex = 3, type = "l", col = 2, ylim = c(-1, 0), lwd = 3,
    xlab = "Proportion mutualists", ylab = "Eigenvalue")
