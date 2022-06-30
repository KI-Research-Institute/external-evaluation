test_that("reweightByMeans works", {
  library(pROC)
  library(WeightedROC)
  library(glmnet)

  # Optimization definitions
  divergence <- 'entropy'
  lambda <- 1e-2
  minW <- 1e-6

  # Load data and model objects
  load(file=system.file('data/internalTrain.RData', package = "ReweightDB"))
  load(file=system.file('data/internalTest.RData', package = "ReweightDB"))
  load(file=system.file('data/external.RData', package = "ReweightDB"))  # External data
  load(file=file.path(system.file('data', package = "ReweightDB"), 'model1.RData')) # model


  # Train and predict in internal and external sets
  xFeatures <- colnames(dExt)[1:(ncol(dExt)-1)]
  p1 <- predict(model1, sapply(dIntTest[xFeatures], as.numeric), type = "response", s = "lambda.1se")[,1]
  intAuc <- auc(roc(dIntTest[['Y']], p1, quiet = TRUE, direction='<'))
  p2 <- predict(model1, sapply(dExt[xFeatures], as.numeric), type = "response", s = "lambda.1se")[,1]
  extAuc <- auc(roc(dExt[['Y']], p2, quiet = TRUE, direction='<'))

  # Compute interactions between features and outcomes and squared features
  dBalanceInt <- computeTable1LikeTransformation(dIntTest, outcomeBalance=TRUE)
  dBalanceExt <- computeTable1LikeTransformation(dExt, outcomeBalance=TRUE)
  muExt <- colMeans(dBalanceExt)

  # Reweight
  w <- reweightByMeans(
    dBalanceInt, muExt, divergence = divergence, lambda = lambda, minW = minW, optimizationMethod='dual',
    verbose = T)

  # Print results
  cat(glue('\nInternal AUC = {format(intAuc, digits=3)}'), '\n')
  cat(glue('\nExternal AUC = {format(extAuc, digits=3)}'), '\n')
  if (!any(is.na(w))) {
    p <- w/length(w)
    cat(glue('Weight entropy = {format(- t(p) %*% log(p), digits=3)}'),'\n')
    dwauc <- WeightedAUC(WeightedROC(p1, dIntTest[['Y']], w))
    cat(glue('Estimated AUC = {format(dwauc, digits=3)}'),'\n')
  } else
    cat('Did not find a feasible solution\n')


  cat('Testing primal optimization\n')
  w <- reweightByMeans(
    dBalanceInt, muExt, divergence = divergence, lambda = lambda, minW = minW, optimizationMethod='primal',
    verbose = T)
  p <- w/length(w)
  cat(glue('Weight entropy = {format(- t(p) %*% log(p), digits=3)}'),'\n')
  pwauc <- WeightedAUC(WeightedROC(p1, dIntTest[['Y']], w))
  cat(glue('Estimated AUC = {format(pwauc, digits=3)}'),'\n')
  expect_equal(format(dwauc, digits=3), '0.716')
  expect_equal(format(pwauc, digits=3), '0.716')
})
