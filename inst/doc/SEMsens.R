## -----------------------------------------------------------------------------
library(lavaan)
library(SEMsens)

## -----------------------------------------------------------------------------
# Lower diagonal correlation matrix in the study by Kim & Schatschneider (2017)
lower = '
1.00
.40 1.00
.40 .64 1.00
.41 .66 .61 1.00
.42 .52 .53 .61 1.00
.34 .50 .46 .53 .48 1.00
.42 .47 .41 .43 .47 .55 1.00
.39 .46 .39 .30 .21 .30 .37 1.00
.24 .31 .30 .31 .26 .32 .27 .56 1.00
.33 .35 .35 .40 .31 .25 .35 .51 .42 1.00
.30 .42 .36 .32 .24 .37 .43 .44 .37 .49 1.00'

# Convert to full covariance matrix
sample.cov = getCov(lower, sds = c(5.64,14.68,6.57,6.07,3.39,10.16,6.11,4.91,15.59,0.96,0.99),
              names = c("Working_memory",
                        "Vocabulary",
                        "Grammar",
                        "Inference",
                        "ToM",
                        "TNL",
                        "Expository",
                        "Spelling",
                        "Sentence_copying",
                        "One_day",
                        "Castle"))

  # The original analytic model
model <-'Vocabulary~Working_memory
  Grammar~Working_memory
  Inference~Vocabulary+Grammar+Working_memory
  ToM~Vocabulary+Grammar+Working_memory
  Spelling~Working_memory
  Sentence_copying~Working_memory
  Discourse~Inference+ToM+Vocabulary+Grammar+Working_memory
  Writing~Spelling+Sentence_copying+Discourse

  Discourse=~TNL+Expository
  Writing=~One_day+Castle

  Vocabulary~~Grammar
  Grammar~~Sentence_copying
  Vocabulary~~Sentence_copying
  Grammar~~Spelling
  Vocabulary~~Spelling
  Inference~~ToM
  Discourse~~Sentence_copying
  Discourse~~Spelling
  Spelling~~Sentence_copying'

  # A sensitivity analysis model template, which additionally includes paths from a phantom
  #    variable to a set of variables (= number of sensitivity parameters) in the analytic model.
sens.model <- 'Vocabulary~Working_memory
  Grammar~Working_memory
  Inference~Vocabulary+Grammar+Working_memory
  ToM~Vocabulary+Grammar+Working_memory
  Spelling~Working_memory
  Sentence_copying~Working_memory
  Discourse~Inference+ToM+Vocabulary+Grammar+Working_memory
  Writing~Spelling+Sentence_copying+Discourse

  Discourse=~TNL+Expository
  Writing=~One_day+Castle

  Vocabulary~~Grammar
  Grammar~~Sentence_copying
  Vocabulary~~Sentence_copying
  Grammar~~Spelling
  Vocabulary~~Spelling
  Inference~~ToM
  Discourse~~Sentence_copying
  Discourse~~Spelling
  Spelling~~Sentence_copying

  Working_memory ~ phantom1*phantom
  Grammar ~ phantom2*phantom
  Vocabulary ~ phantom3*phantom
  ToM ~ phantom4*phantom
  Inference ~ phantom5*phantom
  Spelling ~ phantom6*phantom
  Sentence_copying  ~ phantom7*phantom
  Discourse ~ phantom8*phantom
  Writing ~ phantom9*phantom
  phantom =~ 0 # added for mean of zero
  phantom ~~ 1*phantom'

  # Perform sensitivity analysis
my.sa <-sa.aco(model = model, sens.model = sens.model, sample.cov = sample.cov,
               sample.nobs = 193, n.of.sens.pars = 9, k = 5, max.value= 2000, max.iter = 10,
               opt.fun = 4, ## to just significant
               paths = c(1:18), seed = 1, verbose = FALSE)
# Note: We run with k = 5 and max.iter = 10 for illustration purpose in 5 seconds, 
# please specify them as larger numbers (e.g., default value of k = 50 and mat.iter = 1000)



## -----------------------------------------------------------------------------
# Summary of the sensitivity analysis for each path
sens.tables(my.sa)[[1]]

# Summary of the sensitivity parameters
sens.tables(my.sa)[[2]]

# The sensitivity parameters lead to the minimum coefficient for each
sens.tables(my.sa)[[3]]

# The sensitivity parameters lead to the maximum coefficient for each path
sens.tables(my.sa)[[4]]

# The sensitivity parameters lead to change in significance for each path
sens.tables(my.sa)[[5]]


