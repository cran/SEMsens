#'  Summary of sensitivity analysis results
#'
#' @description  This function can summarize the sensitivity analysis results from
#'    \code{\link{sa.aco}} function.
#'
#' @param expr Returned object of \code{\link{sa.aco}} function.
#' @param sig.level Significance level, default value is 0.05.
#' @param choice Set up the length of summary; default is all.
#' @return
#'     List of summary tables
#'
#' @export sens.tables
#' @examples
#' # see examples in the \code{\link{sa.aco}} function
#'
sens.tables <- function(expr = NULL, sig.level = 0.05, choice = NULL){
  # pull results from the sa.aco function
    old.model.par <- expr$old.model.par
    sens.par <- expr$sens.pars
    outcome <- expr$outcome
    par.est <- expr$model.results
    old.model <- expr$model
    old.model.par <- expr$old.model.par

  n.of.sens <- length(par.est[which(par.est$evals==1 &
                                      par.est$op=="~" & par.est$rhs=="phantom"),][,1])
  par.est$paths <- paste(par.est$lhs, par.est$op, par.est$rhs, sep="")
  phan.names <- par.est[which(par.est$evals==1 &
                               par.est$op=="~" & par.est$rhs=="phantom"),]$paths
  old.model.par$paths <-paste(old.model.par$lhs, old.model.par$op,
                                  old.model.par$rhs, sep="")
  paths <- old.model.par$paths
  if(!is.null(choice)){paths <- paths[choice]}
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
  rownames(sens.summary) <- paths
  colnames(sens.summary) <- c("model.est", "model.pvalue",
                              "mean.est.sens", "min.est.sens", "max.est.sens")
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
   for (i in paths) {
     est.std.data <- par.est[which(par.est$paths==i),]
     pvalue <- old.model.par[which(old.model.par$paths==i),]$pvalue
     old.pvalue <- c(old.pvalue, pvalue)
     if (is.na(pvalue)) {
       phan.sig <- rep(NA, n.of.sens + 1)
       p.paths <- rbind(p.paths, phan.sig)
     } else {
       if (pvalue <= sig.level) {
         if(all(is.na(est.std.data$pvalue))){
           phan.sig <- rep(NA, n.of.sens + 1)
           p.paths <- rbind(p.paths, phan.sig)
         }else{
           if(sum(est.std.data$pvalue[!is.na(est.std.data$pvalue)]>sig.level)>0){
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
             phan.sig <-rep(NA, n.of.sens + 1)
             p.paths <- rbind(p.paths, phan.sig)
         }
         }
         } else {
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
               sig <-  model.sig[which(model.non.sig[,2]=="~" & model.non.sig[,3]=="phantom"),]
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

   return(list(sens.summary = sens.summary, phan.paths = phan.paths,
               phan.min = phan.min, phan.max = phan.max, p.paths = p.paths))
}

