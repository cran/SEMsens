#'Sensitivity Analysis for Structural Equation Modeling Using Simulated
#' Annealing (SA)
#'
#' @description This function can perform sensitivity analysis for
#'     structural equation modeling using simulated annealing (SA)
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
#'     y <-  mean(abs(new.par$est.std[paths] -
#'     old.par$est.std[paths]))/mean(abs(old.par$est.std[paths])); When opt.fun is
#'     2, the optimization function is the standard deviation of deviance
#'     divided by the old estimate
#'     y <-  stats::sd(new.par$est.std[paths] - old.par$est.std[paths])/
#'     mean(abs(old.par$est.std[paths]));
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
#'     sensitivity analysis parameters.
#' @param sample.cov covariance matrix for SEM analysis
#'     when data are not available.
#' @param sample.nobs Number of observations for covariance matrix.
#' @param e Maximum error value used when solution quality used as
#'     the stopping criterion, default is 1e-10.
#' @param n.iter Maximal number of function evaluations within each temperature.
#' @param k  Size of the solution archive, default is 100.
#' @param verbose Print out evaluation process if TRUE, default is TRUE.
#' @param Ntemps Number of temperatures that the algorithm visits. Default value is 10.
#' @param C.criteria Convergence criterion. Default value is 1.
#' @param steepness Steepness of cooling schedule. Default value is 6.
#' @param  measurement Logical. If TRUE, the argument paths will
#'     include measurement paths in the lavaanify format. Default is FALSE.
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
#' @export sa.sa
#'
#' @references
#'   Fisk, C., Harring, J., Shen, Z., Leite, W., Suen, K., & Marcoulides, K.
#'    (2022). Using simulated annealing to investigate sensitivity of
#'    SEM to external model misspecification.
#'    Educational and Psychological Measurement.
#'   <doi:10.1177/00131644211073121>
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
#' mysa <- sa.sa(data = data, model = model,
#' sens.model = sens.model, paths = paths,
#'  n.iter = 3,  Ntemps = 2)
#' # We set Ntemps = 2 and n.iter = 3 to reduce the running time.
#' # You may leave them as default values or specify larger numbers.


sa.sa <- function(data = NULL, sample.cov, sample.nobs, model, sens.model,
                  opt.fun = 1, d = NULL, paths = NULL, verbose = TRUE,
                  n.iter = 10,  e = 1e-10,
                  k = 10, sig.level = 0.05,
                  Ntemps = 10, C.criteria = 1, steepness = 6,
                  measurement = FALSE) {

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
                          sample.nobs = sample.nobs)
  } else {
    old.out = lavaan::sem(model = model, data = data)
  }
  old.par = lavaan::standardizedSolution(old.out, type = "std.all")
  old.fit <- lavaan::fitMeasures(old.out)
  if (is.null(paths)) {paths <- old.par} # if paths are not set, all paths are included
  if (is.character(paths)) {paths <- lavaan::lavaanify(paths, fixed.x = TRUE)} # if paths are in model format

  # initiate a number of k sensitivity analysis models with random sensitivity parameters
  #  sampled from domains
  if (is.null(d)) {
    d <- rep(list(seq(-1, 1, by=.01)), n.of.sens.pars)
  } else {
    if(!is.list(d)) stop("d (domain) must be in a list format; e.g.,
    d =     list(x1 = seq(-1, 1, by=.01),
         x2 = seq(-1, 1, by=.01),
         x3 = seq(-1, 1, by=.01),
         x4 = seq(-1, 1, by=.01),
         x5 = seq(-1, 1, by=.01))")
  }

  Neighbour <- function(x) {
    change.what <- sample.int(length(d), 1)
    x[change.what] <- sample(d[[change.what]], 1)
    return(x)
  }

  s0<- list() # sensitivity parameters
  F0<- list() # objective function values

  #Sens.Param.Vector <-s_prop
  OBJECTIVEfunction <- function(Sens.Param.Vector,
                                opt.fun = opt.fun,
                                paths = paths){
    Sens.Param.Vector <- as.vector(Sens.Param.Vector)
    new.model = sens.model
    for (i in 1:n.of.sens.pars) {
      new.model = gsub(paste("phantom",i, sep=""), paste(Sens.Param.Vector[i]),
                       new.model, ignore.case = FALSE, perl = FALSE,
                       fixed = FALSE, useBytes = FALSE)
    }

      new.model = lavaan::lavaanify(new.model, auto = T, model.type = "sem")

      warnings <- options(warn = 2)
      if (is.null(data)){
        new.out = try(lavaan::sem(model = new.model, sample.cov = sample.cov,
                                  sample.nobs = sample.nobs), silent = TRUE)
      } else {
        new.out = try(lavaan::sem(model = new.model, data = data), silent = TRUE)
      }
      on.exit(options(warnings))
      if(isTRUE(class(new.out)=="try-error")) {return(NA)}


      new.par <- lavaan::standardizedSolution(new.out, type="std.all")


        sens.out <- new.out
        model.1 <- new.par
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

       #if (!is.numeric(opt.fun)){
        #y <-  eval(opt.fun)
      #} else
        if (opt.fun == 1) {
        # if opt.fun==1, we assess the average departure of estimator in the
        #   sensitivity analysis model from the analytic model divided by
        #   the estimator in the analytic model
        y <-  mean(abs(old.par$est.std[paths]), na.rm = TRUE)/
          mean(abs(new.par$est.std[paths]), na.rm = TRUE)
      } else if (opt.fun == 2){
        # if opt.fun==2, we assess the standard deviation of estimate in
        #    the sensitivity analysis model from the analytic model divided by
        #    the estimate in the analytic model
        y <-  stats::sd(new.par$est.std[paths] - old.par$est.std[paths], na.rm = TRUE)/
          mean(abs(old.par$est.std[paths]), na.rm = TRUE)
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
    return(y)
    ### Need to check to see if model converged to avoid the algorithm crashing.
   # if (is.null(y) == FALSE){
    #  return(y)
    #}
    #else{return(NA)}
  }

  # Finds K starting points that do not return NA and assigns them to
  # s0 (parameters) and F0 (objective function values)
  for (i in 1:k) {
    sens.par <- NULL
    for (j in 1:n.of.sens.pars){
      sens.par  <- c(sens.par, sample(d[[j]], 1))
    }
    s_prop <- sens.par
    #opt.fun=1
    f_prop <- OBJECTIVEfunction(Sens.Param.Vector = s_prop, opt.fun = opt.fun,
                                paths = paths)
    while (is.na(f_prop)== T ) {
      sens.par <- NULL
      for (j in 1:n.of.sens.pars){

        sens.par  <- c(sens.par, sample(d[[j]], 1))
      }
      s_prop <- sens.par
      s_prop
      f_prop <-OBJECTIVEfunction(Sens.Param.Vector = s_prop, opt.fun = opt.fun,
                                 paths = paths)
    }
    s0[[i]] <-   s_prop
    F0[[i]] <-   f_prop
  }

  for.sens.paths <- lavaan::lavaanify(sens.model, fixed.x = TRUE)
  phan.names <- for.sens.paths[which(
    for.sens.paths$lhs!="phantom" &
      for.sens.paths$rhs=="phantom"), ]$label

  simulated_annealing_SPEAR <- function(Ntemps, C.criteria, steepness) {

    outTable <- matrix(data=0, nrow = 1, ncol = n.of.sens.pars + 1)
    colnames(outTable) <- c("obj", phan.names)

    s_c <-  list() ###to keep track of the current states
    f_c <-  list()

    s_c <- s0
    f_c<- F0
    #M <- 10
    M<-k
    s_b <-  s_n<- s_c[[M]]
    f_b <-  f_n <- f_c[[M]]



    ## un-comment if you would like to see each step run.
    message("It\tBest\tCurrent\tNeigh\tTemp")
    message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.20f", 0L, f_b, f_c[[M]], f_n, 1))
    #n.iter=20
    #J<-1
    for (J in 1:Ntemps) {


      Temp <- 50/log(1 + J)^steepness - 50/log(1+Ntemps)^steepness

      count<-1

      ## for convergence criteria.
      N1<- (f_c[[M]])+(f_c[[M-1]])+(f_c[[M-2]])+(f_c[[M-3]])
      N2<-(f_c[[M-4]])+(f_c[[M-5]])+(f_c[[M-6]])+(f_c[[M-7]])
      MinVal<-min((f_c[[M]]),(f_c[[M-1]]),(f_c[[M-2]]),(f_c[[M-3]]),
                  (f_c[[M-4]]),(f_c[[M-5]]),(f_c[[M-6]]),(f_c[[M-7]]))



      ### Convergence criteria set as the difference between reps 10:7 and 6:3 < C.criteria*temp
      ### OR if its has been at the same temperature for 300 reps
      while(((f_c[[M]]== MinVal && abs(N1-N2) < C.criteria*Temp) && count>20)==F && (count < n.iter) ){


        s_n <- Neighbour(s_c[[M]])
        f_n <- OBJECTIVEfunction(Sens.Param.Vector = s_n, opt.fun = opt.fun,
                                 paths = paths)


        if( is.na(f_n)==F ){
          # update current state
          if (f_n <= f_c[[M]] ||
              stats::runif(1, 0, 1) < exp(-(f_n - f_c[[M]]) / Temp)){
            s_c[[M+1]] <- s_n
            f_c[[M+1]] <- f_n
            M <- M+1

            ### Update convergence criteria counters
            N1<- (f_c[[M]])+(f_c[[M-1]])+(f_c[[M-2]])+(f_c[[M-3]])
            N2<-(f_c[[M-4]])+(f_c[[M-5]])+(f_c[[M-6]])+(f_c[[M-7]])
            MinVal<-min((f_c[[M]]),(f_c[[M-1]]),(f_c[[M-2]]),(f_c[[M-3]]),
                        (f_c[[M-4]]),(f_c[[M-5]]),(f_c[[M-6]]),(f_c[[M-7]]))


          }

          if ( f_n < f_b ) {
            s_b <- s_n
            f_b <- f_n
          }


        }

        message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.12f", J, f_b, f_c[[M]], f_n, Temp))
        count<- count + 1
      }

      ### At the end of each temperature while loop,
      ### reinitialize the loop using the last 8 elements of the s_c list from the previous loop


      last<- length(s_c) - 7
      s_old<- s_c
      f_old<- f_c
      s_c<-  list()
      f_c<-  list()

      for (i in 1:8) {

        s_c[[i]] <- s_old[[last+(i-1)]]   #### New s_c list for next temperature while loop
        f_c[[i]] <- f_old[[last+ (i-1)]]   ### New f_c list for next temperature while loop
      }

      M <- 8   ### Reset while loop counter
    }

    outTable[1,] <- c(f_b,s_b)


    return(outTable)
  }
  # run the function
  output.SA <- simulated_annealing_SPEAR(Ntemps = Ntemps, C.criteria = C.criteria, steepness = steepness)


  #### Take the solution from simulated annealing and run it through Lavaan to output the final model

  #finalVector<- outTable[1,1:n.of.sens.pars+1]
  finalVector<- output.SA[1,1:n.of.sens.pars+1]
  final.model = sens.model
  for (i in 1:n.of.sens.pars) {
    final.model = gsub(paste("phantom",i, sep=""), paste(finalVector[i]),
                       final.model, ignore.case = FALSE, perl = FALSE,
                     fixed = FALSE, useBytes = FALSE)
  }

  final.model = lavaan::lavaanify(final.model, auto = T, model.type = "sem")



  if (is.null(data)){
    FinalModel = try(lavaan::sem(model = final.model, sample.cov = sample.cov,
                              sample.nobs = sample.nobs), silent = TRUE)
  } else {
    FinalModel = try(lavaan::sem(model = final.model, data = data), silent = TRUE)
  }

  return(list(output.SA, FinalModel))
}
