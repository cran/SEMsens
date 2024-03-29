% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sa.aco.R
\name{sa.aco}
\alias{sa.aco}
\title{Sensitivity Analysis for Structural Equation Modeling Using Ant
Colony Optimization (ACO)}
\usage{
sa.aco(
  data = NULL,
  sample.cov,
  sample.nobs,
  model,
  sens.model,
  opt.fun,
  d = NULL,
  paths = NULL,
  verbose = TRUE,
  max.value = Inf,
  max.iter = 1000,
  e = 1e-10,
  n.of.ants = 10,
  k = 100,
  q = 1e-04,
  sig.level = 0.05,
  rate.of.conv = 0.1,
  measurement = FALSE,
  xi = 0.5,
  seed = NULL,
  ...
)
}
\arguments{
\item{data}{The data set used for analysis.}

\item{sample.cov}{covariance matrix for SEM analysis
when data are not available.}

\item{sample.nobs}{Number of observations for covariance matrix.}

\item{model}{The analytic model of interest.}

\item{sens.model}{Sensitivity analysis model template for
structural equation modeling
with a phantom variable. This is the model of interest
with a phantom variable and sensitivity parameters added.
See examples provided.}

\item{opt.fun}{Customized or preset optimization function.
The argument can be customized as a function, e.g., opt.fun =
quote(new.par$pvalue[paths]-old.par$pvalue[paths]), where new.par and old.par
are the parameter estimates from the sensitivity analysis and analytic models,
respectively.
When opt.fun is
1, the optimization function is the average departure of new estimate
from the old estimate divided by the old estimate
y <-  mean(abs(new.par$est[paths] -
old.par$est[paths]))/mean(abs(old.par$est[paths])); When opt.fun is
2, the optimization function is the standard deviation of deviance
divided by the old estimate
y <-  stats::sd(new.par$est[paths] - old.par$est[paths])/
mean(abs(old.par$est[paths]));
When opt.fun is 3, the optimization function is the average
p value changed or
y <-  mean(abs(new.par$pvalue[paths] - old.par$pvalue[paths]));
When opt.fun is 4, the optimization function is the average distance
from significance level or y <-  mean(abs(new.par$pvalue[paths] -
rep(sig.level,length(paths))));
When opt.fun is 5, we assess the change of RMSEA or
y <-  abs(unname(lavaan::fitmeasures(new.out)["rmsea"]) -
unname(lavaan::fitmeasures(old.out)["rmsea"]));
When opt.fun is 6, we optimize how close RMSEA is to 0.05 or
y <-  1/abs(unname(lavaan::fitmeasures(new.out)["rmsea"]) - 0.05).}

\item{d}{Domains for initial sampling, default is c(-1 ,1) for all
sensitivity analysis parameters. It can
be specified as a list of ranges. For example,
d = list(-0.8, 0.8, -0.9, 0.9) for two
sampling domains with the first from -0.8 to 0.8 and
the second from -0.9 to 0.9.}

\item{paths}{Paths in the model to be evaluated in a sensitivity analysis. If not
specified, all paths will be evaluated. It can be specified in a
numeric format or in a model format. For example, if we evaluate the changes (in p value
or parameter estimation) for paths in an analytic model, we may specify
paths in a model format, e.g.,
paths = 'm ~ x
    y ~ x + m'.
  Or, alternatively, as specify paths = c(1:3) if these paths present in line 1 to 3 in the
  sensitivity analysis model results.}

\item{verbose}{Print out evaluation process if TRUE, default is TRUE.}

\item{max.value}{Maximal value of optimization when used as
the stopping criterion. Default is infinite.}

\item{max.iter}{Maximal number of function evaluations when used as
the stopping criterion.}

\item{e}{Maximum error value used when solution quality used as
the stopping criterion, default is 1e-10.}

\item{n.of.ants}{Number of ants used in each iteration after the initialization
of k converged sensitivity analysis models, default value is 10.}

\item{k}{Size of the solution archive, default is 100.}

\item{q}{Locality of the search (0,1), default is 0.0001.}

\item{sig.level}{Significance level, default value is 0.05.}

\item{rate.of.conv}{The convergence rate threshold for sensitivity analysis models,
default is .10.}

\item{measurement}{Logical. If TRUE, the argument paths will
include measurement paths in the lavaanify format. Default is FALSE.}

\item{xi}{Convergence pressure (0, Inf), suggested: (0,1), default is 0.5.}

\item{seed}{Random seed if specified, default is NULL.}

\item{...}{Additional arguments from the lavaan package.}
}
\value{
Sensitivity analysis results, including the number of evaluations (n.eval),
    number of iterations (n.iter), the maximum value of the objective function (max.y) and
    associated sensitivity parameters values (phantom.coef), analytic model (old.model),
    its results (old.model.par) and fit measures (old.model.fit),
    sensitivity analysis model (sens.model), its fit measures (sens.fit),
    outcome of the objective function (outcome), sensitivity parameters across all
    converged evaluations (sens.pars),
    sensitivity analysis model results (model.results),
    analytic model results (old.out),
    and the first converged sensitivity analysis model results (sens.out).
}
\description{
This function can perform sensitivity analysis for
    structural equation modeling using ant colony optimization (ACO).
}
\examples{
library(lavaan)
# Generate data, this is optional as lavaan also takes variance covariance matrix
sim.model <- ' x =~ x1 + 0.8*x2 + 1.2*x3
               y =~ y1 + 0.5*y2 + 1.5*y3
               m ~ 0.5*x
               y ~ 0.5*x + 0.8*m'
set.seed(10)
data <- simulateData(sim.model, sample.nobs = 1000L)
# standardize dataset
data = data.frame(apply(data,2,scale))

# Step 1: Set up the analytic model of interest
model <- 'x =~ x1 + x2 + x3
          y =~ y1 + y2 + y3
          m ~ x
          y ~ x + m'

# Step 2: Set up the sensitivity analysis model.
#         The sensitivity parameters are phantom1, phantom2, and phantom3 in this example.
sens.model = 'x =~ x1 + x2 + x3
              y =~ y1 + y2 + y3
              m ~ x
              y ~ x + m
              x ~ phantom1*phantom
              m ~ phantom2*phantom
              y ~ phantom3*phantom
              phantom =~ 0 # added for mean of zero
              phantom ~~ 1*phantom' # added for unit variance

# Step 3: Set up the paths of interest to be evaluated in sensitivity analysis.
# Suppose we are interested in all direct and indirect paths.
  paths <- 'm ~ x
            y ~ x + m'

# Step 4: Perform sensitivity analysis
my.sa <- sa.aco(data, model = model, sens.model = sens.model,
                opt.fun = 3, k = 5, #p-value
                paths = paths,
                max.iter = 10)
#Note, please specify larger numbers for k (e.g., 100) and max.iter (e.g., 1000)

# Step 5: Summarize sensitivity analysis results.
# See sens.tables function for explanation of results.
tables <- sens.tables(my.sa)


}
\references{
Leite, W., & Shen, Z., Marcoulides, K., Fish, C., & Harring, J. (accepted).
  Using ant colony optimization for sensitivity analysis in structural equation modeling.
  Structural Equation Modeling: A Multidisciplinary Journal.

  Socha, K., & Dorigo, M. (2008). Ant colony optimization for
  continuous domains. \emph{European Journal of Operational Research,
  185}(3), 1155-1173. <doi:10.1016/j.ejor.2006.06.046>

  Harring, J. R., McNeish, D. M., & Hancock, G. R. (2017).
  Using phantom variables in structural equation modeling to
  assess model sensitivity to external misspecification.
  \emph{Psychological Methods, 22}(4), 616-631.
  <doi:10.1080/10705511.2018.1506925>

  We thank Dr. Krzysztof Socha for providing us the
  ACO code for continuous domains (http://iridia.ulb.ac.be/supp/IridiaSupp2008-001/)
  that the current function is based on.
}
