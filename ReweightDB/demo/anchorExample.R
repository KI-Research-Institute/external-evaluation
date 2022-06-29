rm(list=ls())
library(glmnet)
library(pROC)
library(WeightedROC)
library(glue)

library(ReweightDB)

# Optimization definitions
divergence <- 'entropy'
lambda <- 1e-6
minW <- 1e-6
optimizationMethod <- 'dual'

# Load datasets
load(file=system.file('data/internalTrain.RData', package = "ReweightDB"))
load(file=system.file('data/internalTest.RData', package = "ReweightDB"))
load(file=system.file('data/external.RData', package = "ReweightDB"))


# Train and predict in internal and external sets
xFeatures <- colnames(dExt)[1:(ncol(dExt)-1)]
m <- cv.glmnet(sapply(dIntTrain[xFeatures], as.numeric), dIntTrain[['Y']],
               family = "binomial", type.measure = "auc", alpha = 0.5)
p1 <- predict(m, sapply(dIntTest[xFeatures], as.numeric), type = "response", s = "lambda.1se")[,1]
intAuc <- auc(roc(dIntTest[['Y']], p1, quiet = TRUE, direction='<'))
p2 <- predict(m, sapply(dExt[xFeatures], as.numeric), type = "response", s = "lambda.1se")[,1]
extAuc <- auc(roc(dExt[['Y']], p2, quiet = TRUE, direction='<'))

# Compute interactions between features and outcomes and squared features
dBalanceInt <- computeTable1LikeTransformation(dIntTest, outcomeBalance=TRUE)
dBalanceExt <- computeTable1LikeTransformation(dExt, outcomeBalance=TRUE)
muExt <- colMeans(dBalanceExt)

# Reweight
w <- reweightByMeans(
  dBalanceInt, muExt, divergence = divergence, lambda = lambda, minW = minW, optimizationMethod=optimizationMethod,
  verbose = T)

# Print results
cat(glue('\nInternal AUC = {format(intAuc, digits=3)}'), '\n')
cat(glue('\nExternal AUC = {format(extAuc, digits=3)}'), '\n')
if (!any(is.na(w))) {
  p <- w/length(w)
  cat(glue('Weight entropy = {format(- t(p) %*% log(p), digits=3)}'),'\n')
  wauc <- WeightedAUC(WeightedROC(p1, dIntTest[['Y']], w))
  cat(glue('Estimated AUC = {format(wauc, digits=3)}'),'\n')
} else
  cat('Did not find a feasible solution\n')
