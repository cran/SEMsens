#'  Summary of sensitivity analysis results
#'
#' @description  This function can summarize the sensitivity analysis results from
#'    \code{\link{sa.aco}} function.
#'
#' @param expr Returned object of \code{\link{sa.aco}} function.
#' @param sig.level Significance level, default value is 0.05.
#' @param path Logical, if TRUE, the function only present results for structural paths.
#'             If FALSE, the function will present results for all paths (including
#'             structural paths and measurement paths). Default value is TRUE.
#' @param sort Logical, if TRUE, the function will present sorted results.
#'             If FALSE, the function will present unsorted results. Default value is TRUE.
#' @return
#'     Lists of 5 summary tables. The first table (sens.summary) provides
#'     analytic model results (model path coefficient/model.est, p value/pvalue),
#'     mean, minimum, and maximum values of estimated path coefficients across
#'     all sensitivity analysis models (mean.est.sens, min.est.sens, and max.est.sens).
#'     The second table (phan.paths) provides the summary of sensitivity parameters,
#'     including the mean, minimum, and maximum values of each sensitivity parameters (
#'     mean.phan, min.phan, max.phan). The third table (phan.min) provides the sensitivity
#'     parameters that lead to the minimum path coefficient estimation in a sensitivity
#'     analysis model. The fourth table (phan.max) provides the sensitivity
#'     parameters that lead to the maximum path coefficient estimation in a sensitivity
#'     analysis model. The fifth table (p.paths) provides the sensitivity parameters, if
#'     any, that lead to the change of p value across the significance level.
#'
#'
#' @export sens.tables
#'
#' @references
#'   Leite, W., & Shen, Z., Marcoulides, K., Fish, C., & Harring, J. (in press).
#'   Using ant colony optimization for sensitivity analysis in structural equation modeling.
#'   Structural Equation Modeling: A Multidisciplinary Journal.
#'
#' @examples
#' # see examples in the \code{\link{sa.aco}} function
#'
#'
sens.tables <- function(expr = NULL, sig.level = 0.05, path = TRUE, sort = TRUE){
  # pull results from the sa.aco function
    old.model.par <- expr$old.model.par
    sens.par <- expr$sens.pars
    outcome <- expr$outcome
    par.est <- expr$model.results
    old.model <- expr$model

  n.of.sens <- length(par.est[which(par.est$evals==1 &
                                      par.est$op=="~" & par.est$rhs=="phantom"),][,1])
  par.est$paths <- paste(par.est$lhs, par.est$op, par.est$rhs, sep="")
  phan.names <- par.est[which(par.est$evals==1 &
                               par.est$op=="~" & par.est$rhs=="phantom"),]$paths
  old.model.par$paths <-paste(old.model.par$lhs, old.model.par$op,
                                  old.model.par$rhs, sep="")
  # omitting the report of measurement results if path is TRUE
  if(path){
    paths <- old.model.par[which(old.model.par$op=="~"),]$paths
  } else {
    paths <- old.model.par$paths
  }

  # to produce the first summary table
  sens.summary <-NULL
  for (i in paths) {
    par.est.vector <- par.est[which(par.est$paths==i),]$est.std
    if(all(is.na(par.est.vector))){
      par <- c(old.model.par[which(old.model.par$paths==i),]$est.std,
               old.model.par[which(old.model.par$paths==i),]$pvalue,
               NA, NA, NA)
    }else{
      par <- c(old.model.par[which(old.model.par$paths==i),]$est.std,
               old.model.par[which(old.model.par$paths==i),]$pvalue,
               mean(par.est.vector, na.rm = T),
               min(par.est.vector, na.rm = T),
               max(par.est.vector, na.rm = T))
    }
    sens.summary <- rbind(sens.summary, par)
  }
  colnames(sens.summary) <- c("model.est", "model.pvalue",
                              "mean.est.sens", "min.est.sens", "max.est.sens")
  rownames(sens.summary) <- paths

  # sort the results if TRUE
  sens.summary <- as.data.frame(sens.summary)
  if(sort){
    sens.summary <- sens.summary[order(sens.summary$mean.est.sens),]
  }

  # to produce the second summary table (sensitivity parameters)
  phan.paths <-NULL
  for (i in phan.names) {
    par.est.vector <- par.est[which(par.est$paths==i),]$est.std
    if(all(is.na(par.est.vector))){
      par <- c(NA, NA, NA)
    }else{
      par <- c(mean(par.est.vector, na.rm = T),
               min(par.est.vector, na.rm = T),
               max(par.est.vector, na.rm = T))
    }
   phan.paths <- rbind(phan.paths, par)
  }
  rownames(phan.paths) <- phan.names
  colnames(phan.paths) <- c("mean.phan", "min.phan", "max.phan")

  # sort the results if TRUE
  phan.paths <- as.data.frame(phan.paths)
  if(sort){
    phan.paths <- phan.paths[order(phan.paths$mean.phan),]
  }

  # to produce the third summary table (minimum sensitivity parameters)
  phan.min <-NULL
  for (i in paths){
    est.std.data <- par.est[which(par.est$paths==i),]
    if (all(is.na(est.std.data$est.std))){
      min <- NA
    }else{
      min <- min(est.std.data$est.std, na.rm = T)
    }
    if (!is.na(min)){
    evals <- par.est[which(par.est$est.std==min),]$evals
    if(length(evals)==1){
      eval <- evals
      } else {eval <- sample(evals,1)}
    model.min <- par.est[which(par.est$evals==eval),]
    X.min <- model.min[which(model.min$paths %in% phan.names),]$est.std
    # X.min <-  model.min[which(model.min[,2]=="~" & model.min[,3]=="phantom"),]
    # X.min <- X.min$est.std
    } else {X.min <- rep(NA, n.of.sens)}
    phan.min <- rbind(phan.min, X.min)
  }
   rownames(phan.min) <- paths
   colnames(phan.min) <- phan.names

   # to produce the fourth summary table (maximum sensitivity parameters)
   phan.max <-NULL
   for (i in paths) {
     est.std.data <- par.est[which(par.est$paths==i),]
     if(all(is.na(est.std.data$est.std))){
       max <- NA
     }else{
       max <- max(est.std.data$est.std, na.rm = T)
     }
     if (!is.na(max)){
     evals <- par.est[which(par.est$est.std==max),]$evals
     if(length(evals)==1){
       eval = evals
     } else {eval <- sample(evals,1)}
     model.max <-par.est[which(par.est$evals==eval),]
     X.max <- model.max[which(model.max$paths %in% phan.names),]$est.std
     # X.max <-  model.max[which(model.max[,2]=="~" & model.max[,3]=="phantom"),]
     # X.max <- X.max$est.std
     } else {X.max <- rep(NA, n.of.sens)}
     phan.max <- rbind(phan.max, X.max)
   }
   rownames(phan.max) <- paths
   colnames(phan.max) <- phan.names

   # to produce the fifth summary table (p-values)
   p.paths <- NULL
   old.pvalue <- NULL
   for (i in paths) { # loop through each path
     est.std.data <- par.est[which(par.est$paths==i),] # select the same path over all evaluations
     pvalue <- old.model.par[which(old.model.par$paths==i),]$pvalue # get the original p value
     old.pvalue <- c(old.pvalue, pvalue) # store the old p value
     if (is.na(pvalue)) { # if original p value is NULL, set phan.sig as NAs
       phan.sig <- rep(NA, n.of.sens + 1)
       p.paths <- rbind(p.paths, phan.sig)
     } else { # If original p value exists, check whether p value changes across significance level
       if (pvalue <= sig.level) { # If original p value is less than or equal to significance level
         if(all(is.na(est.std.data$pvalue))){ # if sensitivity model p values are NULL, put NA for sens pars
           phan.sig <- rep(NA, n.of.sens + 1)
           p.paths <- rbind(p.paths, phan.sig)
         }else{ # if there is at least one sensitivity model with a p value,
           if(sum(est.std.data$pvalue[!is.na(est.std.data$pvalue)]>sig.level)>0){
             # If there is at least one non significance p value across all sensitivity model results,
             # we will find the smallest non significant p value and associated sensitivity parameters
             non.sig <- est.std.data[which(est.std.data$pvalue>sig.level),]$pvalue
             non.sig.p <- min(non.sig)
             evals <- est.std.data[which(est.std.data$pvalue==non.sig.p),]$evals
             if(length(evals)==1){
               eval = evals
             } else {eval <- sample(evals,1)}
             model.non.sig <- par.est[which(par.est$evals==eval),]
             non.sig <-  model.non.sig[which(model.non.sig[,2]=="~" & model.non.sig[,3]=="phantom"),]
             non.sig  <- non.sig$est.std
             p.paths <- rbind(p.paths, c(non.sig.p, non.sig))
           } else {
             # If these is no change in significance level, put NAs in sensitivity parameters
             # to represents no sensitivity parameters lead to p value change across significance level
             phan.sig <-rep(NA, n.of.sens + 1)
             p.paths <- rbind(p.paths, phan.sig)
         }
         }
         } else { # if orignial p value is greater than significance level
           if(all(is.na(est.std.data$pvalue))) {
             phan.sig <-rep(NA, n.of.sens + 1)
             p.paths <- rbind(p.paths, phan.sig)
           }else{
             if(sum(est.std.data$pvalue[!is.na(est.std.data$pvalue)]<sig.level)>0){
               sig <- est.std.data[which(est.std.data$pvalue<sig.level),]$pvalue
               sig.p <- max(sig)
               evals <- est.std.data[which(est.std.data$pvalue==sig.p),]$evals
               if(length(evals)==1){
                 eval = evals
               } else {eval <- sample(evals,1)}
               model.sig <- par.est[which(par.est$evals==eval),]
               sig <-  model.sig[which(model.sig[,2]=="~" & model.sig[,3]=="phantom"),]
               sig  <- sig$est.std
               p.paths <- rbind(p.paths, c(sig.p, sig))
             } else {
               phan.sig <-rep(NA, n.of.sens + 1)
               p.paths <- rbind(p.paths, phan.sig)
             }
           }
         }
       }
     }
   p.paths <- cbind(old.pvalue, p.paths)
   rownames(p.paths) <- paths
   colnames(p.paths) <- c("p.value", "p.changed", phan.names)

   # sort the results if TRUE
   p.paths <- as.data.frame(p.paths)
   if(sort){
     p.paths <- p.paths[order(p.paths$p.changed),]
   }


   return(list(sens.summary = sens.summary, phan.paths = phan.paths,
               phan.min = phan.min, phan.max = phan.max, p.paths = p.paths))
}

