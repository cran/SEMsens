## ---- message=FALSE, warning=FALSE--------------------------------------------
# Load lavaan and SEMsens packages
require(lavaan)
require(SEMsens)

## -----------------------------------------------------------------------------

smith10.use <- data(smith19.use, package = "SEMsens")

#Indicate the ordinal variables and then factorize these variables
factornames <- c("P2BEPARN", "P2CHDOES","P2HRDWRM","P2FLTRAP","P2FEELAN","P2CHHARD",
                 "P2MOREWK", "P2HITCHO","P2HITAPO","P2HITPRV","P2HITWAR","P2ATTENB",
                 "P2ATTENP", "P2PARADV","P2PARGRP","P2ATTENS","P2VOLUNT","P2FUNDRS",
                 "P2TVRULE", "P2TVRUL3","P2TVRUL2","P1TELLST","P1SINGSO","P1HELPAR",
                 "P1CHORES", "P1GAMES","P1NATURE","P1BUILD","P1SPORT","P2LIBRAR",
                 "P2CONCRT", "P2MUSEUM","P2ZOO")


## -----------------------------------------------------------------------------
#A model mimicked from the original model in Smith-Adcock et al. (2019). 
model1 <- "   
           # regressions
          C2RSCALE~T2LEARN+T2EXTERN+T2INTERN+T2CONTRO+T2INTERP+f1+f2+f3+f4+f5+f6
          T2LEARN~f1 +f2 +f3 +f4 +f5 +f6
          T2EXTERN~f1 +f2 +f3 +f4 +f5 +f6
          T2INTERN~f1 +f2 +f3 +f4 +f5 +f6
          T2CONTRO~f1 +f2 +f3 +f4 +f5 +f6
          T2INTERP~f1 +f2 +f3 +f4 +f5 +f6

           # latent variables
          f1=~P2BEPARN+P2CHDOES+P2HRDWRM+P2FLTRAP+P2FEELAN+P2CHHARD+P2MOREWK
          f2=~P2HITCHO+P2HITAPO+P2HITPRV+P2HITWAR
          f3=~P2ATTENB+P2ATTENP+P2PARADV+P2PARGRP+P2ATTENS+P2VOLUNT+P2FUNDRS
          f4=~P2TVRULE+P2TVRUL3+P2TVRUL2
      	  f5=~P1TELLST+P1SINGSO+P1HELPAR+P1CHORES+P1GAMES+P1NATURE+P1BUILD+P1SPORT
      	  f6=~P2LIBRAR+P2CONCRT+P2MUSEUM+P2ZOO
            
      	    # residual correlations
      	  T2LEARN ~~ T2EXTERN+T2INTERN+T2CONTRO+T2INTERP
          T2EXTERN ~~ T2INTERN+T2CONTRO+T2INTERP
          T2INTERN ~~ T2CONTRO+T2INTERP
          T2CONTRO ~~ T2INTERP

          f1 ~~ f2+f3+f4+f5+f6
          f2 ~~ f3+f4+f5+f6
          f3 ~~ f4+f5+f6
          f4 ~~ f5+f6
          f5 ~~ f6
         "

#The sensitivity model 
#This sensitivity model used 6 phantom variables for a quicker result
sens.model1 <- '
            # regressions
          C2RSCALE~T2LEARN+T2EXTERN+T2INTERN+T2CONTRO+T2INTERP+f1+f2+f3+f4+f5+f6
          T2LEARN~f1 +f2 +f3 +f4 +f5 +f6
          T2EXTERN~f1 +f2 +f3 +f4 +f5 +f6
          T2INTERN~f1 +f2 +f3 +f4 +f5 +f6
          T2CONTRO~f1 +f2 +f3 +f4 +f5 +f6
          T2INTERP~f1 +f2 +f3 +f4 +f5 +f6

            # latent variable definitions
          f1=~P2BEPARN+P2CHDOES+P2HRDWRM+P2FLTRAP+P2FEELAN+P2CHHARD+P2MOREWK
          f2=~P2HITCHO+P2HITAPO+P2HITPRV+P2HITWAR
          f3=~P2ATTENB+P2ATTENP+P2PARADV+P2PARGRP+P2ATTENS+P2VOLUNT+P2FUNDRS
          f4=~P2TVRULE+P2TVRUL3+P2TVRUL2
          f5=~P1TELLST+P1SINGSO+P1HELPAR+P1CHORES+P1GAMES+P1NATURE+P1BUILD+P1SPORT
          f6=~P2LIBRAR+P2CONCRT+P2MUSEUM+P2ZOO

            # residual correlations
          T2LEARN ~~ T2EXTERN+T2INTERN+T2CONTRO+T2INTERP
          T2EXTERN ~~ T2INTERN+T2CONTRO+T2INTERP
          T2INTERN ~~ T2CONTRO+T2INTERP
          T2CONTRO ~~ T2INTERP

          f1 ~~ f2+f3+f4+f5+f6
          f2 ~~ f3+f4+f5+f6
          f3 ~~ f4+f5+f6
          f4 ~~ f5+f6
          f5 ~~ f6
          
            # phantom variables 
          T2LEARN ~ phantom1*phantom
          T2EXTERN ~ phantom2*phantom
          C2RSCALE ~ phantom3*phantom
          
          f3 ~ phantom4*phantom
          f4 ~ phantom5*phantom
          f5 ~ phantom6*phantom

          phantom =~ 0   # added for mean of zero
          phantom ~~ 1*phantom   # added for unit variance
          '


## -----------------------------------------------------------------------------
old.model = model1
old.out  = sem(model = model1, data = smith19.use, estimator="WLSMV",ordered = factornames)
summary(old.out, standardized=TRUE)
old.model.par <- standardizedSolution(old.out,type = "std.all")
old.model.par
#The lines that we are interested are line 1 to line 11

## ----eval=FALSE, message=FALSE, warning=FALSE---------------------------------
#  my.sa1 <- sa.aco(data= smith19.use,
#                   model = model1,
#                   sens.model = sens.model1,
#                   k = 10,
#                   opt.fun = quote(1/abs(new.par$pvalue[9]-0.05)),
#                   rate.of.conv = .05,
#                   paths = c(1:11),
#                   seed=119,
#                   max.iter = 100,
#                   estimator="WLSMV",
#                   ordered = factornames)
#  #15 out of 175 evaluations converged.
#  

## ---- eval=FALSE, echo=T------------------------------------------------------
#  #error message
#  Error in sa.aco(data= smith19.use,
#                   model = model1,
#                   sens.model = sens.model1,
#                   k = 10,rate.of.conv = .1, :
#  
#                  Sensitivity analysis models do not reach the specified convergence rate.
#                  Please set a lower convergence rate threshold (i.e., rate. of. conv) or
#                  reduce model complexicity)
#  

## ----eval=FALSE---------------------------------------------------------------
#  my.sa.table = sens.tables(my.sa1) # get results

## ----echo=FALSE---------------------------------------------------------------
load("smith19_6_my.sa.table.rdata")

## -----------------------------------------------------------------------------
my.sa.table$sens.summary
my.sa.table$phan.paths
my.sa.table$phan.min
my.sa.table$phan.max
my.sa.table$p.paths

