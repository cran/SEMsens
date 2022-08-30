#' Tabu Search to Optimize Functions of Continuous Variables
#'
#' @description A helper function that implements the main logic of tabu search
#'     to optimize objective functions in continuous domains.
#' @inheritParams gen.neighbors.tabu
#' @param n.var  The dimension of search space.
#' @param max.len The length of the largest hypercube.
#' @param max.tabu.size The maximum size of the tabu list.
#' @param neigh.size The number of neighbors to search for
#'     in each iteration.
#' @param max.iter The maximum number of iterations.
#' @param max.iter.obj The maximum number of successive iterations
#'     without any improvement of the objective function value.
#' @param verbose Logical. Print the current best and overall best
#'     objective functions if TRUE, no printing if FALSE.
#' @return
#' A list with three components:
#' best.param (vector): the best set of parameters found;
#' best.obj (scalar): the value of obj.fun corresponding
#' to best.param; and
#' model.history: the histry of model results.
#' @export sa.tabu.helper
#' @references P., & Berthiau, G. (1997).
#'     Fitting of tabu search to optimize functions of continuous
#'     variables. International journal for numerical methods
#'     in engineering, 40(13), 2449-2457.
#'
sa.tabu.helper <- function(n.var, f, maximum = FALSE, max.len = 1, max.tabu.size = 5, neigh.size = NULL, max.iter = NULL, max.iter.obj = NULL, range = c(-1, 1), r = 1e-5, verbose = TRUE) {

  # Default values
  if (is.null(neigh.size)) {
    neigh.size <- min(n.var * 2, 10)
  }
  if (is.null(max.iter)) {
    max.iter <- n.var * 50
  }
  if (is.null(max.iter.obj)) {
    max.iter.obj <- n.var * 5;
  }

  # Turns warnings into errors
  options(warn = 2)

  tabu.list <- list()
  n.iter <- 1
  n.iter.obj <- 1
  model.history <- list()

  # Initial parameters to some random numbers
  max.attempts <- 50
  for (i in 1:max.attempts) {
    best.param <- current.param <- t(stats::runif(n.var, -1, 1))
    best.obj <- try(f(best.param), silent = TRUE)
    if (class(best.obj)[1] != "try-error") {
      break
    }
  }

  if (class(best.obj)[1] == "try-error") {
    stop("Can't find a valid set of initial parameters! Maybe try a different seed?")
  }
  best.obj <- best.obj$y

  # Add initial parameters to tabu list
  tabu.list[[1]] <- current.param

  if (verbose) {
    cat("  n   curr_obj   best_obj\n")
  }

  while ((n.iter <= max.iter) & (n.iter.obj <= max.iter.obj)) {
    # Find the best neighbor in the current iteration
    best.neighbor <- gen.neighbors.tabu(current.param, maximum, neigh.size, tabu.list, max.len, range, r, f)
    current.param <- best.neighbor$best.param
    current.obj <- best.neighbor$best.obj
    best.neighbor$best.model$evals <- n.iter
    model.history[[n.iter]] <- best.neighbor$best.model

    # Compare the best neighbor in the current iteration to the overall best
    if ((maximum & current.obj > best.obj) || (!maximum & current.obj < best.obj)) {
      best.obj <- current.obj
      best.param <- current.param
      # Reset the number of successive iterations without improvement
      n.iter.obj <- 1
    } else {
      n.iter.obj <- n.iter.obj + 1
    }

    # Push the best neighbor to the end of the tabu list
    tabu.list <- append(tabu.list, list(current.param))

    # Remove the neighbor at the top of the list if list size exceeds the maximum
    if (length(tabu.list) > max.tabu.size) {
      tabu.list <- tabu.list[-1]
    }

    # Print progress
    if (verbose) {
      cat(sprintf("%3d %10f %10f\n", n.iter, current.obj, best.obj))
    }

    n.iter <- n.iter + 1
  }

  return(list(
    best.param = best.param,
    best.obj = best.obj,
    model.history = do.call(rbind, model.history)
  ))
}
