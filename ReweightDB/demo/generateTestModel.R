# Generate a test model for unit tests and vignettes

rm(list=ls())
library(glmnet)
library(glue)

library(ReweightDB)

# Load datasets
load(file=system.file('data/internalTrain.RData', package = "ReweightDB"))
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

# Train and predict in internal and external sets
xFeatures <- colnames(dIntTrain)[1:(ncol(dIntTrain)-1)]
model1 <- cv.glmnet(sapply(dIntTrain[xFeatures], as.numeric), dIntTrain[['Y']],
               family = "binomial", type.measure = "auc", alpha = 0.5)

save(model1, file=file.path(script_dir, '../data/model1.RData'))
