#' Sensitivity Analysis for Structural Equation Modeling Using
#' colony optimization (ACO)
#'
#' @description This function can perform sensitivity analysis for
#'     structural equation modeling using ant colony optimization (ACO).
#'
#' @param data The data set used for analysis.
#' @param model The analytic model of interest.
#' @param sens.model Sensitivity analysis model template for
#'     structural equation modeling
#'     with a phantom variable. This is the model of interest
#'     with a phantom variable and sensitivity parameters added.
#'     See examples provided.
#' @param opt.fun Customized or preset optimization function.
#'     The argument can be customized as a function, e.g., opt.fun =
#'     quote(new.par$pvalue[paths]-old.par$pvalue[paths]), where new.par and old.par
#'     are the parameter estimates from the sensitivity analysis and analytic models,
#'     respectively.
#'     When opt.fun is
#'     1, the optimization function is the average departure of new estimate
#'     from the old estimate divided by the old estimate
#'     y <-  mean(abs(new.par$est[paths] -
#'     old.par$est[paths]))/mean(abs(old.par$est[paths])); When opt.fun is
#'     2, the optimization function is the standard deviation of deviance
#'     divided by the old estimate
#'     y <-  stats::sd(new.par$est[paths] - old.par$est[paths])/
#'     mean(abs(old.par$est[paths]));
#'     When opt.fun is 3, the optimization function is the average
#'     p value changed or
#'     y <-  mean(abs(new.par$pvalue[paths] - old.par$pvalue[paths]));
#'     When opt.fun is 4, the optimization function is the average distance
#'     from significance level or y <-  mean(abs(new.par$pvalue[paths] -
#'     rep(sig.level,length(paths))));
#'     When opt.fun is 5, we assess the change of RMSEA or
#'     y <-  abs(unname(lavaan::fitmeasures(new.out)["rmsea"]) -
#'     unname(lavaan::fitmeasures(old.out)["rmsea"]));
#'     When opt.fun is 6, we optimize how close RMSEA is to 0.05 or
#'     y <-  1/abs(unname(lavaan::fitmeasures(new.out)["rmsea"]) - 0.05).
#'
#' @param paths Paths in the model to be evaluated in a sensitivity analysis. If not
#'     specified, all paths will be evaluated. It can be specified in a
#'     numeric format or in a model format. For example, if we evaluate the changes (in p value
#'     or parameter estimation) for paths in an analytic model, we may specify
#'     paths in a model format, e.g.,
#'     paths = 'm ~ x
#'         y ~ x + m'.
#'       Or, alternatively, as specify paths = c(1:3) if these paths present in line 1 to 3 in the
#'       sensitivity analysis model results.
#' @param sig.level Significance level, default value is 0.05.
#' @param d Domains for initial sampling, default is c(-1 ,1) for all
#'     sensitivity analysis parameters. It can
#'     be specified as a list of ranges. For example,
#'     d = list(-0.8, 0.8, -0.9, 0.9) for two
#'     sampling domains with the first from -0.8 to 0.8 and
#'     the second from -0.9 to 0.9.
#' @param sample.cov covariance matrix for SEM analysis
#'     when data are not available.
#' @param sample.nobs Number of observations for covariance matrix.
#' @param n.of.ants Number of ants used in each iteration after the initialization
#'     of k converged sensitivity analysis models, default value is 10.
#' @param e Maximum error value used when solution quality used as
#'     the stopping criterion, default is 1e-10.
#' @param max.value Maximal value of optimization when used as
#'     the stopping criterion. Default is infinite.
#' @param max.iter Maximal number of function evaluations when used as
#'     the stopping criterion.
#' @param k  Size of the solution archive, default is 100.
#' @param q Locality of the search (0,1), default is 0.0001.
#' @param xi  Convergence pressure (0, Inf), suggested: (0,1), default is 0.5.
#' @param seed Random seed if specified, default is NULL.
#' @param verbose Print out evaluation process if TRUE, default is TRUE.
#' @param rate.of.conv The convergence rate threshold for sensitivity analysis models,
#'     default is .10.
#' @param  measurement Logical. If TRUE, the argument paths will
#'     include measurement paths in the lavaanify format. Default is FALSE.
#' @param ... Additional arguments from the lavaan package.
#'
#' @return
#'     Sensitivity analysis results, including the number of evaluations (n.eval),
#'     number of iterations (n.iter), the maximum value of the objective function (max.y) and
#'     associated sensitivity parameters values (phantom.coef), analytic model (old.model),
#'     its results (old.model.par) and fit measures (old.model.fit),
#'     sensitivity analysis model (sens.model), its fit measures (sens.fit),
#'     outcome of the objective function (outcome), sensitivity parameters across all
#'     converged evaluations (sens.pars),
#'     sensitivity analysis model results (model.results),
#'     analytic model results (old.out),
#'     and the first converged sensitivity analysis model results (sens.out).
#'
#' @export sa.aco
#'
#' @references
#'   Leite, W., & Shen, Z., Marcoulides, K., Fish, C., & Harring, J. (accepted).
#'   Using ant colony optimization for sensitivity analysis in structural equation modeling.
#'   Structural Equation Modeling: A Multidisciplinary Journal.
#'
#'   Socha, K., & Dorigo, M. (2008). Ant colony optimization for
#'   continuous domains. \emph{European Journal of Operational Research,
#'   185}(3), 1155-1173. <doi:10.1016/j.ejor.2006.06.046>
#'
#'   Harring, J. R., McNeish, D. M., & Hancock, G. R. (2017).
#'   Using phantom variables in structural equation modeling to
#'   assess model sensitivity to external misspecification.
#'   \emph{Psychological Methods, 22}(4), 616-631.
#'   <doi:10.1080/10705511.2018.1506925>
#'
#'   We thank Dr. Krzysztof Socha for providing us the
#'   ACO code for continuous domains (http://iridia.ulb.ac.be/supp/IridiaSupp2008-001/)
#'   that the current function is based on.
#'
#'
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
#'
#' # Step 3: Set up the paths of interest to be evaluated in sensitivity analysis.
#' # Suppose we are interested in all direct and indirect paths.
#'   paths <- 'm ~ x
#'             y ~ x + m'
#'
#' # Step 4: Perform sensitivity analysis
#' my.sa <- sa.aco(data, model = model, sens.model = sens.model,
#'                 opt.fun = 3, k = 5, #p-value
#'                 paths = paths,
#'                 max.iter = 30)
#' #Note, please specify larger numbers for k (e.g., 100) and max.iter (e.g., 1000)
#'
#' # Step 5: Summarize sensitivity analysis results.
#' # See sens.tables function for explanation of results.
#' tables <- sens.tables(my.sa)
#'

#'
sa.aco <- function(data = NULL, sample.cov, sample.nobs, model, sens.model,
                   opt.fun, d = NULL, paths = NULL, verbose = TRUE,
                   max.value = Inf, max.iter = 1000,  e = 1e-10,
                   n.of.ants = 10, k = 100, q = 0.0001, sig.level = 0.05,
                   rate.of.conv = .1, measurement = FALSE,
                   xi = 0.5,
                   seed = NULL, ...) {

  # Calculate the number of sensitivity parameters
  for.n.of.sens.pars <- lavaan::lavaanify(sens.model, fixed.x = TRUE)
  n.of.sens.pars <-length(for.n.of.sens.pars[which(
    for.n.of.sens.pars$lhs!="phantom" &
      for.n.of.sens.pars$rhs=="phantom"), ]$lhs)
  if (n.of.sens.pars < 2)
    stop ("Sensitivity model must have at least two sensitivity parameters or phantom coefficients.")

  # Run analytic model and pull results
   if (is.null(data)){
    old.out = lavaan::sem(model = model, sample.cov = sample.cov,
                              sample.nobs = sample.nobs, ...)
  } else {
    old.out = lavaan::sem(model = model, data = data, ...)
  }
  old.par = lavaan::standardizedSolution(old.out, type = "std.all")
  old.fit <- lavaan::fitMeasures(old.out)
  if (!is.null(seed)) {set.seed(seed)} # set up seed if any
  if (is.null(paths)) {paths <- old.par} # if paths are not set, all paths are included
  if (is.character(paths)) {paths <- lavaan::lavaanify(paths, fixed.x = TRUE)} # if paths are in model format

  e.abs <- e # absolute error
  e.rel <- e # relative error

  # initiate parameters
  eval <- 0
  iter <- 0
  last.impr <- max.iter
  nl <- matrix(NA, k, k-1)
  sens.pars <- data.frame()
  outcome <- vector()
  model.results <- data.frame()
  max.X <- rep(NA, n.of.sens.pars)
  max.y <- -Inf
  p.X <- vector()
  sens.fit <- vector()
  p <- data.frame(v = numeric(), sd = numeric(), gr = numeric());

  # initiate a number of k sensitivity analysis models with random sensitivity parameters
  #  sampled from domains
  if (is.null(d)) {d <- list(rep(c(-1, 1), n.of.sens.pars))} else {
    if(!is.list(d)) stop("d (domain) must be in a list format; e.g.,
    d = list(-1, 1,
             -1, 1,
             -1, 1,
             -1, 1)")
  }

  if (rate.of.conv <= 0 | rate.of.conv > 1)
  stop ("Convergence rate (rate.of.conv) must be in (0, 1]")
  for (i in 1:(round(1/rate.of.conv*k, 0))) {
    X <- vector()
    for (j in 1:n.of.sens.pars) {  # sample sensitivity parameters from domains
      X <- c(X, stats::runif(1, d[[1]][2*j-1], d[[1]][2*j]))
    }
    X <- t(X)
    new.model = sens.model
    for (l in 1:n.of.sens.pars) {
      new.model = gsub(paste("phantom", l, sep = ""), paste(X[l]), new.model,
                       ignore.case = FALSE, perl = FALSE,
                       fixed = FALSE, useBytes = FALSE)
    }

    iter <- iter + 1
    if (verbose) {cat('Number of tried evaluations is ', iter, ".\n", sep = "")}
    warnings <- options(warn = 2)
    if (is.null(data)){
      new.out = try(lavaan::sem(model = new.model, sample.cov = sample.cov,
                                sample.nobs = sample.nobs, ...), silent = TRUE)
    } else {
      new.out = try(lavaan::sem(model = new.model, data = data, ...), silent = TRUE)
    }
    if(isTRUE(class(new.out)=="try-error")) { next }
    on.exit(options(warnings))

    new.par = lavaan::standardizedSolution(new.out, type="std.all")
    eval <- eval + 1
    if (verbose) {cat('Number of converged evaluations is ', eval, ".\n", sep = "")}
    new.par$lines <- 1:length(new.par[, 1])
    new.par$evals <- eval
    model.results <- rbind(model.results, new.par)
    if (eval == 1) {
      sens.out <- new.out
      model.1 <- model.results
      model.1$path <- paste(model.1$lhs, model.1$op, model.1$rhs, sep="")
      phan.names <- model.1[which(model.1$evals == 1 &
                                    model.1$op=="~" & model.1$rhs=="phantom"),]$path
      if(is.data.frame(paths)){
        if (measurement){
          paths <- which(model.1$lhs %in% paths$lhs & model.1$rhs %in% paths$rhs)
        } else {
          paths <- which(model.1$lhs %in% paths$lhs & model.1$op =="~" & model.1$rhs %in% paths$rhs)
        }
      }
    }
    sens.par <- c(X, eval = eval)
    sens.pars <- rbind(sens.pars, sens.par)
    fit <- c(lavaan::fitMeasures(new.out), eval = eval)
    sens.fit <- rbind(sens.fit, fit)
    if (!is.numeric(opt.fun)){
      y <-  eval(opt.fun)
    } else if (opt.fun == 1) {
        # if opt.fun==1, we assess the average departure of estimator in the
        #   sensitivity analysis model from the analytic model divided by
        #   the estimator in the analytic model
        y <-  mean(abs(old.par$est[paths]), na.rm = TRUE)/
          mean(abs(new.par$est[paths]), na.rm = TRUE)
      } else if (opt.fun == 2){
        # if opt.fun==2, we assess the standard deviation of estimate in
        #    the sensitivity analysis model from the analytic model divided by
        #    the estimate in the analytic model
        y <-  stats::sd(new.par$est[paths] - old.par$est[paths], na.rm = TRUE)/
          mean(abs(old.par$est[paths]), na.rm = TRUE)
      } else if (opt.fun == 3) {
        # if opt.fun==3, we assess the average p-value changed
        y <-  mean(abs(new.par$pvalue[paths] - old.par$pvalue[paths]), na.rm = TRUE)
      } else if (opt.fun == 4){
        # if opt.fun==4, we assess the average distance of p-value from the significance level
        y <-  1 / mean(abs(new.par$pvalue[paths] - rep(sig.level,length(paths))), na.rm = TRUE)
      } else if (opt.fun == 5){
        # if opt.fun==5, we assess the change of RMSEA
        y <-  abs(unname(lavaan::fitmeasures(new.out)["rmsea"]) -
                    unname(lavaan::fitmeasures(old.out)["rmsea"]))
      } else if (opt.fun == 6){
        # if opt.fun==6, we optimize how close RMSEA is to 0.05
        y <-  1/abs(unname(lavaan::fitmeasures(new.out)["rmsea"]) - 0.05)
      }
    outcome <- c(outcome, y)
    p.X <- rbind(p.X, X)
    p <- rbind(p, data.frame(v = y, sd = 0, gr = 0))
    if (eval == k){break} # break out if we have a number of k converged models
  }

  # The algorithm will stop if sensitivity analysis models converge below
  # the convergence rate
  if (length(p.X)==0 | length(p$v) < k) stop(
  "\n
   Sensitivity analysis models do not reach the specified convergence rate.
   Please set a lower convergence rate threshhold (i.e., rate.of.conv) or reduce model complexicity")

  p$gr <- rank(-p$v, ties.method = "random")
  for (i in 1:k){
    nl[i,] <- (1:k)[1:k!=i]
  }

  while (TRUE) { # the algorithm will stop if one of the stop criteria is met
    dist.mean <- p.X
    # the algorithm will stop if it converges (no more change in sensitivity parameter values)
    if (sum(apply(dist.mean, 2, stats::sd)) == 0) {
      colnames(sens.pars) <- c(phan.names, "eval")
      return(list(n.eval = eval, n.iter = iter, max.y = max.y,
                             phantom.coef = max.X, old.model.par = old.par, old.model.fit = old.fit,
                  model = model, sens.model = sens.model, sens.fit = sens.fit,
                  outcome = outcome, sens.pars = sens.pars,
                  model.results = model.results, old.out = old.out, sens.out = sens.out))
    }
    dist.rank <- p$gr
    dim(dist.mean) <- c(length(p$v), n.of.sens.pars)
    o.X <- vector()
    o.X <- gen.sens.pars(dist.mean, dist.rank, n.of.ants, nl,  q, k, xi)
    # the algorithm will stop if it converges (no more available random samples)
    if (length(o.X) == 0) {
      colnames(sens.pars) <- c(phan.names, "eval")
      return(list(n.eval = eval, n.iter = iter, max.y = max.y,
                             phantom.coef = max.X, old.model.par = old.par, old.model.fit = old.fit,
                  model = model, sens.model = sens.model, sens.fit = sens.fit,
                  outcome = outcome, sens.pars = sens.pars,
                  model.results = model.results, old.out = old.out, sens.out = sens.out))
    }
    X <- o.X
    dim(X) <- c(length(X)/n.of.sens.pars, n.of.sens.pars)
    for (j in 1:dim(X)[1]) { # refit the models for n.of.ants times
      X.sens <- X[j, ]
      X.model <- as.vector(X.sens)
      new.model = sens.model
      for (i in 1:dim(X)[2]) { # refit the model with sensitivity analysis parameters
        new.model = gsub(paste("phantom",i, sep=""), paste(X.model[i]),
                         new.model, ignore.case = FALSE, perl = FALSE,
                         fixed = FALSE, useBytes = FALSE)
      }

      iter <- iter + 1
      if (verbose) {cat('Number of tried evaluations is ', iter, ".\n", sep = "")}
      warnings <- options(warn = 2)
      on.exit(options(warnings))
      if (is.null(data)){
        new.out = try(lavaan::sem(model = new.model, sample.cov = sample.cov,
                                  sample.nobs = sample.nobs, ...), TRUE)
      } else {
        new.out =  try(lavaan::sem(model = new.model, data = data, ...), TRUE)
      }
      if(isTRUE(class(new.out)!="try-error")) {

      new.par <- lavaan::standardizedSolution(new.out, type="std.all")
        eval <- eval + 1
        if (verbose) {cat('Number of converged evaluations is ', eval, ".\n", sep = "")}
        p.X <- rbind(p.X, X.sens)
        new.par$lines <- 1:length(new.par[,1])
        new.par$evals <- eval
        model.results <- rbind(model.results, new.par)
        fit <- c(lavaan::fitMeasures(new.out), eval = eval)
        sens.fit <- rbind(sens.fit, fit)
        sens.par <- c(X.sens, eval = eval)
        sens.pars <- rbind(sens.pars, sens.par)
        if (!is.numeric(opt.fun)){
          y <-  eval(opt.fun)
        } else if (opt.fun == 1) {
          # if opt.fun ==1, we assess the average departure of estimator in
          #   the sensitivity analysis model from the analytic model divided by
          #   the estimator in the analytic model
          y <-  mean(abs(old.par$est[paths]), na.rm = TRUE)/
            mean(abs(new.par$est[paths]), na.rm = TRUE)
        } else if (opt.fun == 2){
          # if opt.fun==2, we assess the standard deviation of estimator in
          #   the sensitivity analysis model departs from the analytic model divided by the estimator in
          #   the analytic model
          y <-  stats::sd(new.par$est[paths] - old.par$est[paths], na.rm = TRUE)/
            mean(abs(old.par$est[paths]), na.rm = TRUE)
        } else if (opt.fun == 3) {
          # if opt.fun==3, we assess the average p-value changed
          y <-  mean(abs(new.par$pvalue[paths] - old.par$pvalue[paths]), na.rm = TRUE)
        } else if (opt.fun == 4){
          # if opt.fun==4, we assess the average distance from 0.05
          y <-  1 / mean(abs(new.par$pvalue[paths] - rep(sig.level, length(paths))), na.rm = TRUE)
        } else if (opt.fun == 5){
          # if opt.fun==5, we assess the change of RMSEA
          #        y[j] <-  abs(new.out@test$standard$pvalue - old.out@test$standard$pvalue)
          y <-  abs(unname(lavaan::fitmeasures(new.out)["rmsea"]) -
                      unname(lavaan::fitmeasures(old.out)["rmsea"]))
        } else if (opt.fun == 6){
          # if opt.fun==6, we assess how far RMSEA is to the significance level
           y <-  abs(unname(lavaan::fitmeasures(new.out)["rmsea"]) - 0.05)
        }
        outcome <- c(outcome, y)
        p <- rbind(p, data.frame(v = y, sd = 0, gr = 0))

      p$gr <- rank(-p$v, ties.method = "random") # calculate the rank of the solutions
      idx.final <- p$gr <= k
      p <- p[idx.final,]
      p.X <- p.X[idx.final,]
      dim(p.X) <- c(length(p.X)/n.of.sens.pars, n.of.sens.pars)
    }}

    # recalculate the ranking
    p$gr <- rank(-p$v, ties.method="random")
    for (i in 1:k) {nl[i,] <- (1:k)[1:k!=i]}

    # check if the required accuracy have been obtained
    if (max(outcome, na.rm = TRUE) > max.y) {
      max.y <- max(outcome, na.rm = TRUE)
      max.X <- sens.pars[which.max(outcome), ]
      colnames(max.X) <- c(phan.names, "eval")
      last.impr <- eval}

    if ((abs(max.y - max.value) < abs(e.rel * max.value + e.abs)) |
          (max.y > max.value)) {
        colnames(sens.pars) <- c(phan.names, "eval")
        return(list(n.eval = eval, n.iter = iter, max.y = max.y,
                               phantom.coef = max.X, old.model.par = old.par, old.model.fit = old.fit,
                    model = model, sens.model = sens.model, sens.fit = sens.fit,
                    outcome = outcome, sens.pars = sens.pars,
                    model.results = model.results, old.out = old.out, sens.out = sens.out))
    }

    # check if the maximum allowed number of objective function
    # evaluations has not been exceeded
    if (max.iter > 0 & iter >= max.iter) {
      colnames(sens.pars) <- c(phan.names, "eval")
      return(list(n.eval = eval, n.iter = iter, max.y = max.y,
                  phantom.coef = max.X, old.model.par = old.par, old.model.fit = old.fit,
                  model = model, sens.model = sens.model, sens.fit = sens.fit,
                  outcome = outcome, sens.pars = sens.pars,
                  model.results = model.results, old.out = old.out, sens.out = sens.out))
    }
  }
}



