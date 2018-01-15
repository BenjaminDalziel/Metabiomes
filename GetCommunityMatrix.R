GetCommunityMatrix <- function(S, C, sigma, Pm, Pc, s) {
  
  # Return a community interaction matrix sampled per the arguments
  #
  # S       Number of taxa
  # C       Community connectivity i.e. proportion nonzero links
  # sigma   Standard deviation of interaction strength
  # Pm      Proportion cooperative interactions (+/+)
  # Pc      Proportion competitive interactions (-/-)
  # s       Strength of intraspecific competition
  #
  # The method follows Coyte et al. 2015, and can reproduce the criteria used in
  # May 1972 as a special case, by setting Pm = Pc = 0.25
  #
  # Coyte, K. Z., Schluter, J., & Foster, K. R. (2015). Science, 350(6261), 663–666.
  # http://doi.org/10.1126/science.aad2602
  #
  # May, R. M. (1972). Will a large complex system be stable? Nature, 238(5364), 413–414.
  
  A <- matrix(0, S, S)
  diag(A) <- -s
  
  for (i in 1:S) {
    for (j in 1:(i - 1)) {
      if (runif(1) <= C) {
        p <- runif(1)
        a <- sigma * abs(rnorm(1))
        b <- sigma * abs(rnorm(1))
        
        if (p > 1 - Pm) {
          A[i, j] <- a
          A[j, i] <- b
          
        } else {
          if (p <= Pc) {
            A[i, j] <- -a
            A[j, i] <- -b
            
          } else {
            if (runif(1) <= 0.5) {
              A[i, j] <-  a
              A[j, i] <- -b
            } else {
              A[i, j] <- -a
              A[j, i] <-  b
            }
            
          }
        }
      }
    }
  }
  
  return(A)
  
}
