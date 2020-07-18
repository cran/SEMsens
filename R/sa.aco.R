#' Sensitivity Analysis for Structural Equation Modeling Using ACO
#'
#' @description This function can perform sensitivity analysis for
#'     structural equation modeling using ant colony optimization (ACO).
#'
#' @param data The data set used for analysis.
#' @param model The analytic model of interest set up as a lavaan format.
#' @param sens.model Sensitivity analysis model template for
#'     structural equation modeling
#'     with a phantom variable. This is the model of interest
#'     with phantom variable and sensitivity parameters added.
#'     See examples provided.
#' @param n.of.sens.pars number of sensitivity parameters
#'    added in the sens.model.
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
#'     y <-  mean(abs(new.par$pvalue[paths] - old.par$pvalue[paths]))
#'
#'     When opt.fun is 4, the optimization function is the average distance
#'     from significance level or y <-  mean(abs(new.par$pvalue[paths] -
#'     rep(sig.level,length(paths))))#'
#'
#' @param paths Paths in the model to be evaluated in a sensitivity analysis.
#' @param sig.level Significance level, default value is 0.05.
#' @param d Domains for initial sampling, default is c(-1 ,1) for all. It can
#'     be specified as a list of ranges (e.g., d = list(-1, 1, -1, 1) for two
#'     sampling domains).
#' @param sample.cov covariance matrix
#' @param sample.nobs Number of observations for covariance matrix
#' @param n.of.ants Number of ants used in each iteration after the initialization
#'     of k length, default value is 10.
#' @param e Maximum error value used when solution quality used as
#'     stopping criterion, default is 1e-10.
#' @param max.value Maximal value of optimization when used as
#'     the stopping criterion
#' @param max.iter Maximal number of function evaluations when used as
#'     the stopping criterion
#' @param k  Size of the solution archive, default is 50.
#' @param q Locality of the search (0,1), default is 0.0001
#' @param xi  Convergence pressure (0,Inf), suggested: (0,1), default is 0.5
#' @param seed Random seed if specified, default is NULL.
#' @param verbose Print out evaluation process if true, default is TRUE.
#'
#' @return
#'     Sensitivity analysis results
#'
#' @export sa.aco
#'
#' @references
#'   Socha, K., & Dorigo, M. (2008). Ant colony optimization for
#'   continuous domains. \emph{European Journal of Operational Research,
#'   185}(3), 1155-1173.
#'
#'   Harring, J. R., McNeish, D. M., & Hancock, G. R. (2017).
#'   Using phantom variables in structural equation modeling to
#'   assess model sensitivity to external misspecification.
#'   \emph{Psychological Methods, 22}(4), 616.
#'
#'   We thank Dr. Krzysztof Socha for providing us the
#'   original code (http://iridia.ulb.ac.be/supp/IridiaSupp2008-001/)
#'   that the current function is based on.
#'
#'
#' @examples
#' library(lavaan)
#' # generate data
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
#' # Step 2: Set up sensitivity analysis model
#' #         the sensitivity parameters are phantom1, phantom2 and phantom3
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
#' # Step 3: check the analytic model results and decide parameter of interests
#' #         for sensitivity analysis
#' old.model = lavaan::lavaanify(model = model, auto = TRUE, model.type="sem")
#' old.out = lavaan::lavaan(model = old.model, data = data)
#' old.par = lavaan::standardizedSolution(old.out, type="std.all")
#' old.par # we are interested in lines 7, 8 and 9 for the indirect and direct effects
#'
#' # Step 4: perform sensitivity analysis
#' my.sa <- sa.aco(data, model = model, sens.model = sens.model,
#'                 n.of.sens.pars = 3, k = 5,
#'                 opt.fun = quote(1/abs(new.par$pvalue[9]-0.05)), #p-value
#'                 paths = 9,
#'                 max.iter = 10)
#'    # Note: We run with k = 5 and max.iter = 10 for illustration purpose in 5 seconds,
#'   #  please specify them as larger numbers (e.g., default value of k = 50 and mat.iter = 1000)
#'
#'
#' # Step 5: summarize sensitivity analysis results
#' tables <- sens.tables(my.sa)
#' tables
#'
sa.aco <- function(data = NULL, sample.cov, sample.nobs, model, sens.model, n.of.sens.pars = NULL,
                   opt.fun, d = NULL, paths = NULL, verbose = TRUE,
                   max.value = Inf, max.iter = 1000,  e = 1e-10,
                   n.of.ants = 10, k = 50, q = 0.0001, sig.level = 0.05,
                   xi = 0.5,
                   seed = NULL) {
  # run analytic model and pull results
  model = lavaan::lavaanify(model = model, auto = TRUE, model.type = "sem")
  if (is.null(data)){
    old.out = lavaan::sem(model = model, sample.cov = sample.cov,
                              sample.nobs = sample.nobs)
  } else {
    old.out = lavaan::lavaan(model = model, data = data)
  }
  old.par = lavaan::standardizedSolution(old.out, type = "std.all")
  old.fit <- lavaan::fitMeasures(old.out)
  if (!is.null(seed)) {set.seed(seed)} # set up seed if any
  if (is.null(paths)) {paths <- 1: length(old.par[, 1])} # if paths are not set, all paths are included
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

  for (i in 1:(100*k)) {
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
    new.model = lavaan::lavaanify(new.model, auto = T, model.type="sem")
    iter <- iter + 1
    if (verbose) {cat('Number of tried evaluations is ', iter, ".\n", sep = "")}
#    options(warn = 2)
    warnings <- options(warn = 2)
    on.exit(options(warnings))
    if (is.null(data)){
      new.out = try(lavaan::sem(model = new.model, sample.cov = sample.cov,
                                sample.nobs = sample.nobs), silent = TRUE)
      if(isTRUE(class(new.out)=="try-error")) { next }
    } else {
      new.out = try(lavaan::lavaan(model = new.model, data = data), TRUE)
      if(isTRUE(class(new.out)=="try-error")) { next }
    }
#    options(warn = 0)
    new.par = lavaan::standardizedSolution(new.out, type="std.all")
    eval <- eval + 1
    if (verbose) {cat('Number of converged evaluations is ', eval, ".\n", sep = "")}
    new.par$lines <- 1:length(new.par[, 1])
    new.par$evals <- eval
    model.results <- rbind(model.results, new.par)
    if (eval == 1) {
      model.1 <- model.results
      model.1$path <- paste(model.1$lhs, model.1$op, model.1$rhs, sep="")
      phan.names <- model.1[which(model.1$evals == 1 &
                                    model.1$op=="~" & model.1$rhs=="phantom"),]$path
#      colnames(X) <- phan.names
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
        y <-  mean(abs(old.par$est[paths]))/mean(abs(new.par$est[paths]))
      } else if (opt.fun == 2){
        # if opt.fun==2, we assess the standard deviation of estimate in
        #    the sensitivity analysis model from the analytic model divided by
        #    the estimate in the analytic model
        y <-  stats::sd(new.par$est[paths] - old.par$est[paths])/
          mean(abs(old.par$est[paths]))
      } else if (opt.fun == 3) {
        # if opt.fun==3, we assess the average p-value changed
        y <-  mean(abs(new.par$pvalue[paths] - old.par$pvalue[paths]))
      } else if (opt.fun == 4){
        # if opt.fun==4, we assess the average distance of p-value from the significance level
        y <-  1 / mean(abs(new.par$pvalue[paths] - rep(sig.level,length(paths))))
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

  p$gr <- rank(-p$v, ties.method = "random")
  for (i in 1:k){
    nl[i,] <- (1:k)[1:k!=i]
  }

  while (TRUE) { # the algorithm will stop if one of the stop criteria is met
    dist.mean <- p.X
    # the algorithm will stop if it converges (no more change in sensitivity parameter value)
    if (sum(apply(dist.mean, 2, stats::sd)) == 0) {
      colnames(sens.pars) <- c(phan.names, "eval")
#      colnames(max.X) <- c(phan.names, "eval")
      return(list(n.eval = eval, n.iter = iter, max.y = max.y,
                             phantom.coef = max.X, old.model.par = old.par, old.model.fit = old.fit,
                  model = model, sens.model = sens.model, sens.fit = sens.fit,
                  outcome = outcome, sens.pars = sens.pars,
                  model.results = model.results))
    }
    dist.rank <- p$gr
    dim(dist.mean) <- c(length(p$v), n.of.sens.pars)
#    count <- 0
    o.X <- vector()
    o.X <- gen.sens.pars(dist.mean, dist.rank, n.of.ants, nl,  q, k, xi)
    # the algorithm will stop if it converges (no more available random samples)
    if (length(o.X) == 0) {
      colnames(sens.pars) <- c(phan.names, "eval")
#      colnames(max.X) <- c(phan.names, "eval")
      return(list(n.eval = eval, n.iter = iter, max.y = max.y,
                             phantom.coef = max.X, old.model.par = old.par, old.model.fit = old.fit,
                  model = model, sens.model = sens.model, sens.fit = sens.fit,
                  outcome = outcome, sens.pars = sens.pars,
                  model.results = model.results))
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
      new.model = lavaan::lavaanify(new.model, auto = T, model.type = "sem")
      iter <- iter + 1
      if (verbose) {cat('Number of tried evaluations is ', iter, ".\n", sep = "")}
      warnings <- options(warn = 2)
      on.exit(options(warnings))
      if (is.null(data)){
        new.out = try(lavaan::sem(model = new.model, sample.cov = sample.cov,
                                  sample.nobs = sample.nobs), TRUE)
        if(isTRUE(class(new.out)=="try-error")) { next }
      } else {
        new.out =  try(lavaan::lavaan(model = new.model, data = data), TRUE)
        if(isTRUE(class(new.out)=="try-error")) { next }
      }
#      options(warn = 0)
      new.par <- lavaan::standardizedSolution(new.out, type="std.all")
      #if (new.out@optim$converged==TRUE) {}
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
          y <-  mean(abs(old.par$est[paths]))/mean(abs(new.par$est[paths]))
        } else if (opt.fun == 2){
          # if opt.fun==2, we assess the standard deviation of estimator in
          #   the sensitivity analysis model departs from the analytic model divided by the estimator in
          #   the analytic model
          y <-  stats::sd(new.par$est[paths] - old.par$est[paths])/mean(abs(old.par$est[paths]))
        } else if (opt.fun == 3) {
          # if opt.fun==3, we assess the average p-value changed
          y <-  mean(abs(new.par$pvalue[paths] - old.par$pvalue[paths]))
        } else if (opt.fun == 4){
          # if opt.fun==4, we assess the average distance from 0.05
          y <-  1 / mean(abs(new.par$pvalue[paths] - rep(sig.level, length(paths))))
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
    }

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
                    model.results = model.results))
      }
    # check if the maximum allowed number of objective function
    # evaluations has not been exceeded
    if (max.iter > 0 & iter >= max.iter) {
      colnames(sens.pars) <- c(phan.names, "eval")
      return(list(n.eval = eval, n.iter = iter, max.y = max.y,
                  phantom.coef = max.X, old.model.par = old.par, old.model.fit = old.fit,
                  model = model, sens.model = sens.model, sens.fit = sens.fit,
                  outcome = outcome, sens.pars = sens.pars,
                  model.results = model.results))
    }
  }
}



