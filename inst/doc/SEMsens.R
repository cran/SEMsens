## ---- message=FALSE, warning=FALSE--------------------------------------------
# Load lavaan and SEMsens packages
require(lavaan)
require(SEMsens)
set.seed(1)

## -----------------------------------------------------------------------------
# STEP 1: Prepare data:
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
 
# STEP 2: Set up analytic model and sensitivity anaysis model
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

  # A sensitivity analysis model template, which additionally includes paths
  #  from a phantom variable to a set of variables (= number of sensitivity parameters)
  #   in the analytic model.
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

# STEP 3: Set up the paths of interest to be evaluated in sensitivity analysis.
  paths <- 'Vocabulary~Working_memory
  Grammar~Working_memory
  Inference~Vocabulary+Grammar+Working_memory
  ToM~Vocabulary+Grammar+Working_memory
  Spelling~Working_memory
  Sentence_copying~Working_memory
  Discourse~Inference+ToM+Vocabulary+Grammar+Working_memory
  Writing~Spelling+Sentence_copying+Discourse'

  # STEP 4: Perform sensitivity analysis.
my.sa <-sa.aco(model = model, sens.model = sens.model, sample.cov = sample.cov,
               sample.nobs = 193, k = 50, max.value= 2000, max.iter = 100,
               opt.fun = 4, ## from significant to just significant
               paths = paths, seed = 1, verbose = FALSE)
# We set up a max iteration of 100 and solution archive length of 50 for 
# illustration purpose. Please specify a larger number of iteration (e.g., 1000),
#  and a larger k (e.g., 100).


## -----------------------------------------------------------------------------
my.table <- sens.tables(my.sa)
# Table 1: Summary of the sensitivity analysis for each path
my.table[[1]]

# Table 2:  Summary of the sensitivity parameters
my.table[[2]]

# Table 3: The sensitivity parameters lead to the minimum coefficient for each
my.table[[3]]

# Table 4: The sensitivity parameters lead to the maximum coefficient for each path
my.table[[4]]

# Table 5: The sensitivity parameters lead to change in significance for each path
my.table[[5]]


