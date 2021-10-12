## ---- message=FALSE, warning=FALSE--------------------------------------------
# Load lavaan and SEMsens packages
require(lavaan)
require(SEMsens)
set.seed(1)

## -----------------------------------------------------------------------------
# Correlation matrix from published data
cormat <- diag(1,29) # start with the diagonal
cormat[lower.tri(cormat)] <- c( # add lower triangle data
  .46,.22,.30,-.20,-.30,-.20,-.22,-.33,-.19,-.25,-.63,-.60,-.48,-.17,-.42,-.09,-.36,.09,
  .04,-.12,0,-.20,-.20,-.43,-.33,-.49,-.18,-.03,.21,.30,-.15,-.25,-.10,-.20,-.31,-.15,
  -.24,-.40,-.37,-.32,-.08,-.24,-.02,-.27,.06,.07,-.08,0,-.10,-.16,-.25,-.18,-.35,-.06,
  -.05,.04,-.05,-.09,-.12,-.10,-.01,.03,.02,-.18,-.14,-.09,-.02,-.12,-.05,-.08,.03,.07,
  -.08,-.06,-.10,-.07,-.10,-.03,-.11,.04,.03,-.07,-.24,0,-.12,-.50,-.32,-.45,-.38,-.29,
  -.31,-.16,-.24,-.16,-.37,.12,.05,-.28,-.13,-.06,-.03,-.16,-.12,-.23,-.08,.06,.29,.31,
  .10,.08,-.02,.07,.21,.24,.25,.06,.09,.04,.14,-.02,-.01,.10,.04,.04,.13,.19,.03,.18,
  .14,-.04,.38,.41,.23,.13,.22,.36,.35,.35,.09,.14,.11,.14,-.09,-.04,.11,.04,.11,.15,
  .30,.17,.33,.12,-.02,.31,.05,.03,.15,.22,.20,.21,-.05,.10,-.05,.03,-.04,-.06,.03,.05,
  .03,.07,.24,.19,.23,.13,-.08,.17,.15,.19,.28,.23,.23,.04,.18,.08,.10,-.09,-.06,.07,
  .01,.11,.08,.21,.20,.28,.06,-.03,.43,.56,.35,.32,.26,.16,.27,.13,.45,-.04,-.04,.16,
  .07,.14,.15,.23,.12,.29,.09,-.07,.67,.24,.26,.27,.15,.21,.13,.25,-.10,-.02,.12,0,.12,
  .10,.18,.21,.27,.17,-.06,.36,.32,.31,.11,.29,.16,.40,-.10,-.05,.22,.10,.17,.15,.27,
  .23,.33,.17,-.08,.78,.77,.11,.49,.09,.40,-.23,-.14,.19,-.07,.29,.27,.57,.40,.68,.27,0,
  .82,.10,.57,.08,.42,-.29,-.10,.13,-.12,.33,.31,.66,.46,.75,.32,.02,.10,.49,.08,.38,
  -.31,-.11,.18,-.06,.33,.34,.54,.38,.67,.31,-.08,.14,.52,.20,.22,.31,.42,.39,.01,.01,
  .07,.06,.13,-.03,.06,.14,.42,-.29,-.12,.24,-.06,.25,.25,.38,.28,.48,.18,-.01,.16,.16,
  .28,.68,.48,.06,.10,.01,.12,.15,-.05,.04,.02,.03,.24,.11,.14,.16,.31,.22,.32,.09,-.01,
  .37,.12,.18,-.20,-.15,-.22,-.08,-.23,-.11,.07,.13,.18,-.07,0,-.10,.04,-.06,-.03,.15,
  .54,.11,.12,.07,.10,.16,-.01,-.02,-.06,-.01,-.09,-.05,-.10,-.14,.06,.57,.28,.17,.36,
  .27,-.05,.32,.21,.36,.18,.05,.42,.68,.26,-.05,.55,.23,.04,.28,-.03,-.07
)
cormat[upper.tri(cormat)] <- t(cormat)[upper.tri(cormat)] # mirror lower to upper

# Standard deviations from published data
sds <- diag(c(1,.38,.73,.65,.58,.78,.75,.55,.59,.54,.45,.74,.85,.80,.57,
              .74,.69,.42,.57,.63,.62,.70,1.01,.94,.65,.37,.73,.40,1))

# Create covariance matrix covmat with row and column names
covmat <- sds %*% cormat %*% sds
rownames(covmat) <- colnames(covmat) <- c(
  "FACES","Duty","Team","Hours","Absent","Effort","Sick","Quality","Family",
  "Community","Personal","Thoughts","SI","QI","Stress","Jobs","Costs",
  "Benefits","Joblessness","Moving","Impact","Interference","V1","V2",
  "Prepare","Look","Intensity","Resign","Unemployed"
)

## -----------------------------------------------------------------------------
# Original model from Hom and Kinicki (2001)
model <- "JSat =~ FACES + Duty + Team
  IC =~ Family + Personal + Community
  JAvoid =~ Quality + Absent + Effort + Sick
  WC =~ QI + Thoughts + SI
  WEU =~ Stress + Benefits + Impact + Jobs
  JSearch =~ Prepare + Look + Intensity
  CA =~ V1 + V2
  Turnover =~ Resign
  Resign ~~ 0*Resign
  UR =~ Unemployed
  Unemployed ~~ 0*Unemployed
  
  JSat ~ IC
  JAvoid ~ JSat
  WC ~ JAvoid + JSat + IC
  WEU ~ WC + JSat + UR
  JSearch ~ WEU
  CA ~ JSearch
  Turnover ~ CA + WC + UR"

## -----------------------------------------------------------------------------
# Sensitivity model
sens.model <- "JSat =~ FACES + Duty + Team
  IC =~ Family + Personal + Community
  JAvoid =~ Quality + Absent + Effort + Sick
  WC =~ QI + Thoughts + SI
  WEU =~ Stress + Benefits + Impact + Jobs
  JSearch =~ Prepare + Look + Intensity
  CA =~ V1 + V2
  Turnover =~ Resign
  Resign ~~ 0*Resign
  UR =~ Unemployed
  Unemployed ~~ 0*Unemployed
  
  JSat ~ IC
  JAvoid ~ JSat
  WC ~ JAvoid + JSat + IC
  WEU ~ WC + JSat + UR
  JSearch ~ WEU
  CA ~ JSearch
  Turnover ~ CA + WC + UR

  IC ~ phantom1*phantom
  UR ~ phantom2*phantom
  JSat ~ phantom3*phantom
  JAvoid ~ phantom4*phantom
  WC ~ phantom5*phantom
  WEU ~ phantom6*phantom
  JSearch ~ phantom7*phantom
  CA ~ phantom8*phantom
  Turnover ~ phantom9*phantom

  phantom =~ 0
  phantom ~~ 1*phantom"

## -----------------------------------------------------------------------------
ptab = lavaanify(model = model, auto = TRUE, model.type="sem")
ptab[26:40,1:4]

## ---- message=FALSE, results="hide"-------------------------------------------
my.sa <- sa.aco(model = model, sens.model = sens.model, sample.cov = covmat,
                sample.nobs = 410, opt.fun = 3, paths = 27:39, 
                seed = 1, k = 5, max.iter = 20)

## -----------------------------------------------------------------------------
my.sa.results = sens.tables(my.sa, path=TRUE) # get results
my.sa.results = lapply(my.sa.results, round, digits = 3) # round results (optional)

## -----------------------------------------------------------------------------
my.sa.results$sens.summary

## -----------------------------------------------------------------------------
my.sa.results$phan.paths

## -----------------------------------------------------------------------------
my.sa.results$phan.min

## -----------------------------------------------------------------------------
my.sa.results$phan.max

## -----------------------------------------------------------------------------
my.sa.results$p.paths

