#' A Function that Generates a List of Neighbors in Tabu Search
#'
#' @description This function generates a list of neighbors for tabu
#'     search that can be used in tabu search for sensitivity analysis.
#'
#' @param current.param The center of the hypercubes.
#' @param maximum Logical. Maximize the objective function if TRUE,
#'     minimize the objective function if FALSE.
#' @param neigh.size Number of neighbors to search for in each
#'     iteration.
#' @param tabu.list The list of tabu.
#' @param max.len The length of the largest hypercube.
#' @param range The range for the parameter space in the tabu search.
#' @param r Radius of a tabu ball.
#' @param f The objective function to be optimized.
#' @param max.attempts The maximum number of attempts to find a neighbor that
#'    is not near the points in the tabu list, default is 10.
#' @return A list of information about the best neighbor, including
#'     best parameters, objective function values, and the model.
#' @export gen.neighbors.tabu

gen.neighbors.tabu <- function(current.param, maximum,
                               neigh.size, tabu.list,
                               max.len, range, r, f,
                               max.attempts = 10) {

  # initialization
  best.param <- NULL
  best.obj <- ifelse(maximum, -Inf, Inf)
  n.var <- length(current.param)
  neighbors <- matrix(nrow = neigh.size, ncol = n.var)

  # linear partitioning: len_{n+m-1} = (len_{n} (n-m+1)) / n for m = 1,...,n
  len <- max.len / neigh.size # length of hypercube

  for (i in 1:neigh.size) {
    n.attempts <- 1
    is.invalid <- TRUE
    while (is.invalid) {
      # generate samples from the uniform distribution over [0, len]
      # and use them as coordinates in a hypercube, then
      # translate the coordinates so that they are as if sampled from
      # the hypercube that is centered at current.param
      x <- stats::runif(n.var, 0, len) + current.param - len / 2
      if (i > 1) {
        # to ensure the sampled neighbor is completely outside the previous
        # square, we just need to make sure the sample range of one dimension
        # excludes the sample range of the previous square, the other
        # dimensions are free from this constraint

        # randomly pick a dimension to have the sample restriction
        j <- sample(1:n.var, 1)
        prev.len <- len / 2
        if (stats::runif(1) < 0.5) {
          x[j] <- stats::runif(1, current.param[j] - len / 2, current.param[j] - prev.len / 2)
        } else {
          x[j] <- stats::runif(1, current.param[j] + prev.len / 2, current.param[j] + len / 2)
        }
      }

      # force the neighbors to stay inside the specified range
      x <- pmin(x, range[2])
      x <- pmax(x, range[1])

      # check if the neighbor is in any of the tabu balls
      distance <- vapply(tabu.list, function(center) sum((x - center)^2), numeric(1))
      # the neighbor is invalid if the radius of the tabu ball is larger
      # than the distance between the current sample and the center of
      # each tabu ball
      is.invalid <- any(r >= distance)
      if (!is.invalid) {
        # try to evaluate the objective function using the neighbor
        obj.fun.res <- try(f(x), silent = TRUE)
        # the neighbor is invalid if the objective function raises an error
        is.invalid <- class(obj.fun.res)[1] == "try-error"
      }

      # increment the number of times that the neighbor generated
      # in this hypercube zone is inside a tabu ball
      n.attempts <- n.attempts + 1

      # halve the radius of the tabu balls if cannot find
      # a suitable neighbor in max.attempts iterations
      if (is.invalid & n.attempts > max.attempts) {
        n.attempts <- 1
        r <- r / 2
      }
    }

    obj <- obj.fun.res$y
    if ((maximum & obj > best.obj) || (!maximum & obj < best.obj)) {
      best.param <- x
      best.obj <- obj
      best.model <- obj.fun.res$model
    }

    # increase the length of the edges of the hypercube
    len <- len + 1 / neigh.size
  }
  return(list(best.param = best.param, best.obj = best.obj, best.model = best.model))
}
