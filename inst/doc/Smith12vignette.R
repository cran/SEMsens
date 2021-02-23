## ----message=F----------------------------------------------------------------
#Load the packages
require(SEMsens)
require(lavaan)

## -----------------------------------------------------------------------------
#Set a correlation matrix
lower = '
1.00
0.68 1.00
0.54 0.55 1.00
0.65 0.63 0.67 1.00
0.33 0.37 0.68 0.54 1.00
-0.01 0.00 0.03 -0.04 0.07 1.00'

#convert to full covariance matrix, using function from lavaan
full = getCov(lower, sds= c(4.61,5.37,7.25,3.44,8.91,8.80),
              names = c("PU", "PEU", "SP", "CI","SOC","Gains"))


## -----------------------------------------------------------------------------
# Original model
lav_model <-  'SP~SOC
Gains~SP
PU~SP+PEU
PEU~SP
CI~SP+PU+PEU+SOC
Gains ~~ 0*CI 
'
# Fit the original model with sem function
modelFit <-  sem(lav_model, sample.nobs=517, sample.cov=full, fixed.x=TRUE, std.lv=TRUE)
summary(modelFit, standardized = TRUE) #look at Std.all
fitMeasures(modelFit)

## -----------------------------------------------------------------------------
smith_original <- lavaan::lavaanify(model = lav_model, auto = TRUE, model.type = "sem", fixed.x = TRUE)
smith_original <- lavaan::lavaan(model = smith_original, sample.cov = full, sample.nobs = 517)
smith_original_par <- lavaan::standardizedSolution(smith_original, type = "std.all")
smith_original_par #4th row and 7th column of table : smith_original_par[1:4,1:7]

## -----------------------------------------------------------------------------

# Sensitivity model, with sensitivity parameters for all variables
sens_model <-  'SP~SOC
    Gains ~ SP
    PU ~ SP+PEU
    PEU ~ SP
    CI ~ SP+PU+PEU+SOC
    Gains ~~ 0*CI
    SP ~ phantom1*phantom
    Gains ~ phantom2*phantom
    PU ~ phantom3*phantom
    PEU ~ phantom4*phantom
    CI ~ phantom5*phantom
    SOC ~ phantom6*phantom
    phantom =~ 0  #mean of zero
    phantom ~~ 1*phantom  # variance of one'

## -----------------------------------------------------------------------------
smith_example <- sa.aco(
  sample.cov = full,
  sample.nobs = 517,
  model = lav_model,
  sens.model = sens_model, 
  opt.fun = 1,
  paths = c(1:9), 
  max.iter = 20, 
  k = 5)

## -----------------------------------------------------------------------------
smith_tables <- sens.tables(smith_example)
smith_tables$sens.summary

## -----------------------------------------------------------------------------
smith_tables$phan.paths

## -----------------------------------------------------------------------------
smith_tables$phan.min

## -----------------------------------------------------------------------------
smith_tables$phan.max

## -----------------------------------------------------------------------------
smith_tables$p.paths

