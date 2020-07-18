#' Generate Sensitivity Parameters
#'
#' @description This function can generate a set of path coefficients from
#'     a phantom variable to variables in a structural equation model
#'     based on given distributions of the rank of optimization target
#'     (with probability of using a distribution based on its rank).
#'
#' @inheritParams sa.aco
#' @param dist.mean List of means - coordinates
#' @param dist.rank Rank of the archived values of optimization object
#' @param nl Neighborhood of the search area
#' @return
#'     Generated sensitivity parameter values (i.e., a matrix with n.of.ants
#'     rows and n.of.sens.pars columns)
#'
#' @export gen.sens.pars
#
#' @references
#'   Socha, K., & Dorigo, M. (2008). Ant colony optimization for
#'   continuous domains. European Journal of Operational Research,
#'   185(3), 1155-1173.
#'
#'   We thank Dr. Krzysztof Socha for providing us the
#'   original code (http://iridia.ulb.ac.be/supp/IridiaSupp2008-001/)
#'   for this function.
#'
#' @examples
#'
#' k <- 50 # size of archive
#' # Generate dist.mean and dist.rank
#' dist.mean <- cbind(rnorm(k), rnorm(k), rnorm(k), rnorm(k), rnorm(k))
#' y <- rowMeans(dist.mean)
#' dist.rank <- rank(-y, ties.method = "random")
#' # set up neighborhood
#' nl <- matrix(NA, k, k-1)
#' for (i in 1:k){
#'   nl[i,] <- (1:k)[1:k != i]
#'  }
#' my.sens.pars <- gen.sens.pars(dist.mean, dist.rank, n.of.ants = 10,
#'                               nl, q = 0.0001, k =50, xi = 0.50)
#' my.sens.pars
#'
gen.sens.pars <- function(dist.mean, dist.rank, n.of.ants, nl,
                          q = 0.0001, k = 50, xi = 0.50) {
  euc.dist <- function(d) { # Euclidean distance
    return(sqrt(sum(d^2)))
  }
  X <- array(dim = c(n.of.ants, dim(dist.mean)[2]))
  idx <- sample(dim(dist.mean)[1], size = n.of.ants,
    replace = TRUE, prob = stats::dnorm(dist.rank, 1, q*k))

  # iterate through the chosen distributions
  for (l in 1:length(idx)) {
    j <- idx[l]
    # rotate the coordinate system
    o.dist.mean <- t(t(dist.mean) - dist.mean[j,])  # transaltion of origin
    r.dist.mean <- o.dist.mean
    set <- nl[j,]  # set of available neighbours
    vec <- vector()
    for (m in 1:(dim(dist.mean)[2]-1)) {
      dis <- apply(matrix(r.dist.mean[set,m:dim(r.dist.mean)[2]],
        length(set),length(m:dim(r.dist.mean)[2])),1, euc.dist)
      if (sum(dis)==0.0)  return(NULL) # if the distribution have converged
      if (length(set)>1){ ## ADDED bracket
        choice <- sample(set,size=1,prob=dis^4)
      } else { ## ADDED bracket
        choice <- set
      }
      vec <- cbind(vec,o.dist.mean[choice,])
      R <- qr.Q(qr(vec), complete=TRUE) # rot. matrix after orthogonalization
      if (det(R)<0) {
        R[,1] <- -R[,1]
      }
      r.dist.mean <-  o.dist.mean %*% R # rotated coordinates
      set <- set[set!=choice]
    }

    dist.sd <- vector()
    for (i in 1:dim(dist.mean)[2]) {
      dist.sd <- c(dist.sd,sum(abs(r.dist.mean[nl[j,],i]-
        r.dist.mean[j,i]))/(k-1))
    }
    n.x <- stats::rnorm(dim(dist.mean)[2],r.dist.mean[j,],dist.sd*xi)
    n.x <- R %*% n.x
    n.x <- t(n.x + dist.mean[j,])
    X[l,] <- n.x
  }
  return(X)
}

