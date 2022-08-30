#' Sensitivity Analysis for SEM using Tabu Search
#'
#' @description This function conducts sensitivity analysis
#'     for SEM using tabu search.
#' @inheritParams gen.neighbors.tabu
#' @inheritParams sa.tabu.helper
#' @inheritParams sa.aco
#' @return
#'     A list with five components:
#'     model: The old model;
#'     old.model.par: Parameters of the old model;
#'     model.results: Sensitivity analysis model results;
#'     best.param: Parameters that optimize the objective function;
#'     best.obj: The optimized objective function value;
#'     sens.par: NULL. Included for compatibility;
#'     outcome: NULL. Included for compatibility.
#' @export sa.tabu
#' @examples
#' library(lavaan)
#' # Generate data, this is optional as lavaan also takes variance covariance matrix
#' sim.model <- ' x =~ x1 + 0.8*x2 + 1.2*x3
#'                y =~ y1 + 0.5*y2 + 1.5*y3
#'                m ~ 0.5*x
#'                y ~ 0.5*x + 0.8*m'
#' set.seed(10)
#' data <- simulateData(sim.model, sample.nobs = 1000L)
#' # standardize dataset
#' data = data.frame(apply(data,2,scale))
#'
#' # Step 1: Set up the analytic model of interest
#' model <- 'x =~ x1 + x2 + x3
#'           y =~ y1 + y2 + y3
#'           m ~ x
#'           y ~ x + m'
#'
#' # Step 2: Set up the sensitivity analysis model.
#' #         The sensitivity parameters are phantom1, phantom2, and phantom3 in this example.
#' sens.model = 'x =~ x1 + x2 + x3
#'               y =~ y1 + y2 + y3
#'               m ~ x
#'               y ~ x + m
#'               x ~ phantom1*phantom
#'               m ~ phantom2*phantom
#'               y ~ phantom3*phantom
#'               phantom =~ 0 # added for mean of zero
#'               phantom ~~ 1*phantom' # added for unit variance
#' # Step 3: Set up the paths of interest to be evaluated in sensitivity analysis.
#' # Suppose we are interested in all direct and indirect paths.
#'   paths <- 'm ~ x
#'             y ~ x + m'
#' # Step 4: Perform sensitivity analysis
#' out <- sa.tabu(model = model,
#'                sens.model = sens.model,
#'                data = data,
#'                opt.fun = 1,
#'                max.iter = 2,
#'                max.iter.obj = 2)
#'  # Note, please specify larger numbers for
#'  # max.iter (e.g., 50) and max.iter.obj (e.g., 10)
#'
#' # Step 5: Summarize sensitivity analysis results.
#' # See sens.tables function for explanation of results.
#' tables <- sens.tables(out)
#'

sa.tabu <- function(model, sens.model,
                    data = NULL, sample.cov = NULL,
                    sample.nobs = NULL,
                    opt.fun = 1,
                    sig.level = 0.05, ...) {

  init.model <- model

  # Find non-phantom paths and store their names
  init.model.par.table <- lavaan::lavaanify(init.model, auto = T, model.type = "sem", fixed.x = TRUE)
  non.phan.path.ids <- which(init.model.par.table$op == "~")
  non.phan.path.names <- character(length(non.phan.path.ids))
  for (i in seq_along(non.phan.path.ids)) {
    j <- non.phan.path.ids[i]
    non.phan.path.names[i] <- paste(
      init.model.par.table$lhs[j],
      init.model.par.table$op[j],
      init.model.par.table$rhs[j]
    )
  }

  # Find phantom paths and store their names
  sens.model.par.table <- lavaan::lavaanify(sens.model, auto = T, model.type = "sem", fixed.x = TRUE)
  phan.path.ids <- which(sens.model.par.table$label != "")
  phan.path.names <- character(length(phan.path.ids))
  for (i in seq_along(phan.path.ids)) {
    j <- phan.path.ids[i]
    phan.path.names[i] <- paste(
      sens.model.par.table$lhs[j],
      sens.model.par.table$op[j],
      sens.model.par.table$rhs[j]
    )
  }

  # Fit initial model
  init.model.sem <- lavaan::sem(model = init.model.par.table, data = data, sample.cov = sample.cov, sample.nobs = sample.nobs)
  init.model.params <- lavaan::standardizedSolution(init.model.sem, type = "std.all")

  # Run sensitivity analysis
  sens.model.template <- sens.model
  # Define objective function
  f <- function(phantom.coef) {
    # Replace the placeholders with phantom.coef
    for (j in 1:length(phantom.coef)) {
      sens.model.template <- gsub(
        paste0("phantom", j),
        paste(phantom.coef[j]),
        sens.model.template
      )
    }
    # Fit an SEM model with the given phantom.coef
    sens.model.template.par.table <- lavaan::lavaanify(sens.model.template, auto = T, model.type = "sem", fixed.x = TRUE)
    sens.model.sem <- try(lavaan::sem(model = sens.model.template.par.table, data = data, sample.cov = sample.cov, sample.nobs = sample.nobs), silent = TRUE)
    #if (class(sens.model.sem)[1] == "try-error") {
    #  stop("Error!")
    #}
    sens.model.params <- lavaan::standardizedSolution(sens.model.sem, type = "std.all")

    if (opt.fun == 1) {
      # Ratio between the mean absolute coefficients of initial model and
      # new model
      y <- mean(abs(
        sens.model.params$est[non.phan.path.ids] -
          init.model.params$est[non.phan.path.ids]), na.rm = TRUE) /
        mean(abs(init.model.params$est[non.phan.path.ids]), na.rm = TRUE)
    } else if (opt.fun == 2) {
      # Standard deviation of the difference between the coefficients in
      # the initial model and the new model, standardized by the mean
      # absolute values of the coefficients in the initial model
      y <- stats::sd(sens.model.params$est[non.phan.path.ids] -
                init.model.params$est[non.phan.path.ids], na.rm = TRUE) /
        mean(abs(init.model.params$est[non.phan.path.ids]), na.rm = TRUE)
    } else if (opt.fun == 3) {
      # Average change in p-values
      y <- mean(abs(sens.model.params$pvalue[non.phan.path.ids] -
                      init.model.params$pvalue[non.phan.path.ids]), na.rm = TRUE)
    } else if (opt.fun == 4) {
      # Average distance of p-values from the significance level
      y <- mean(abs(sens.model.params$pvalue[non.phan.path.ids] -
                      rep(sig.level, length(non.phan.path.ids))), na.rm = TRUE)
    } else if (opt.fun == 5) {
      # Change in RMSEA
      y <- abs(unname(lavaan::fitmeasures(sens.model.sem)["rmsea"]) -
                 unname(lavaan::fitmeasures(init.model.sem)["rmsea"]))
    } else if (opt.fun == 6) {
      # Distance of RMSEA in new model from 0.05
      y <- abs(unname(lavaan::fitmeasures(sens.model.sem)["rmsea"]) - 0.05)
    }
    return(list(y = y, model = sens.model.params))
  }
  res <- sa.tabu.helper(length(phan.path.ids), f, maximum = TRUE, ...)
  colnames(res$best.param) <- phan.path.names

  out <- list(
    model = model,
    old.model.par = init.model.params,
    model.results = res$model.history,
    best.param = res$best.param[1, ],
    best.obj = res$best.obj,
    sens.par = NULL,
    outcome = NULL
  )
  return(out)
}
