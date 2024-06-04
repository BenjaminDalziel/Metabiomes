SampleCommunityMatrices <- function(M, matrix_para) {
    # Return an array of community interaction matrices
    # with a specified level of universality
    #
    # Arguments:
    # M       Number of communities
    #
    # Matrix parameters for each community:
    # (see SampleCommunityMatrix.R for details)
    #
    # S       Number of taxa
    # C       Community connectivity in prototype matrix
    # sigma   Standard deviation of interaction strength in prototype matrix
    # Pm      Proportion cooperative interactions (+/+) in prototype matrix
    # Pc      Proportion competitive interactions (-/-) in prototype matrix
    # s       Strength of intraspecific competition in prototype matrix
    # delta:  Universality

    #
    # Remarks:
    # The method follows Coyte et al. 2015.
    # SampleCommunityMatrix is the workhorse function.
    # If delta is NULL the community matrices are drawn independently.
    # If delta is not NULL the method of Bashan et al. is used.

    # References:
    # Coyte, K. Z., Schluter, J., & Foster, K. R. (2015).
    # Science, 350(6261), 663–666. http://doi.org/10.1126/science.aad2602

    # Bashan, A., Gibson, T.E., Friedman, J., Carey, V.J.,
    # Weiss, S.T., Hohmann, E.L. and Liu, Y.Y. (2016).
    # Universality of human microbial dynamics. Nature, 534(7606),
    # 259-262. http://doi.org/10.1038/nature18301


    S <- matrix_para$S
    C <- matrix_para$C
    sigma <- matrix_para$sigma
    Pm <- matrix_para$Pm
    Pc <- matrix_para$Pc
    s <- matrix_para$s
    
    A <- array(data = NA, dim = c(S, S, M))
    if (is.null(matrix_para$delta)) {
        for (i in 1:M) {
            A[, , i] <- SampleCommunityMatrix(matrix_para)
        }
    } else {
        
        delta <- matrix_para$delta
        
        B <- SampleCommunityMatrix(matrix_para)

        for (i in 1:M) {
            u <- runif(1, 1 - delta, 1 + delta)
            thisA <- u * B
            diag(thisA) <- -s
            A[, , i] <- thisA
        }
    }

    return(A)
}
