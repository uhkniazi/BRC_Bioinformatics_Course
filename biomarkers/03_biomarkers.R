# Name: 03_biomarkers.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 26/04/2018
# Desc: perform variable selection for biomarker discovery

setwd('biomarkers/')

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

## load the training and test data sets
load('results/lData.train.rds')
str(lData.train)

load('results/lData.test.rds')
str(lData.test)

## some gene names are not common so remove those from training data
i = match(rownames(lData.test$data), rownames(lData.train$data))
# sanity check
identical(rownames(lData.train$data)[i], rownames(lData.test$data))
lData.train$data = lData.train$data[i,]

## merge the HC and LTB groups
fGroups = rep('ATB', times=length(lData.train$grouping))
fGroups[lData.train$grouping != 'ATB'] = 'Other'
fGroups = factor(fGroups, levels = c('Other', 'ATB'))
table(fGroups, lData.train$grouping)

## perform Random Forest based variable selection step
dfData = data.frame(t(lData.train$data))
str(dfData)
dfData[1:5, 1:5]
## adjust boot.num as desired
set.seed(123) # for replication purposes
oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100, big.warn = F)
save(oVar.r, file='results/oVar.r.rds')

# plot the top 20 genes based on importance scort with 95% confidence interval for standard error
plot.var.selection(oVar.r)

# get the variables
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
str(dfRF)
head(dfRF)

## how does the conditional distribution look like, 
## grouping the scores for the genes based on the coefficient of variation
coplot(ivMean ~ ivSD | group.lab, data=dfRF, pch=19, xlab=c('Standard Error of Mean Score', 'Given Coef of Variation'),
       ylab='Mean Score', cex=0.5)

# select the top 30 variables
cvTopGenes = rownames(dfRF)[1:30]

# use the top 30 genes to find top combinations of genes
dfData = data.frame(t(lData.train$data[cvTopGenes, ]))

set.seed(123)
oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

## guesstimate of model size before we start to overfit
log(nrow(dfData))

# print variable combinations
for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
}

#### check the performance using 10 fold cross validation
## 10 fold nested cross validation with various variable combinations
par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  dfData.train = data.frame(t(lData.train$data))
  dfData.train = data.frame(dfData.train[,cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  
  oCV = CCrossValidation.LDA(test.dat = dfData.train, train.dat = dfData.train, test.groups = fGroups,
                             train.groups = fGroups, level.predict = 'ATB', boot.num = 50)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}

## choose the 1 variable model
dfData = data.frame(t(lData.train$data))
dfData = data.frame(dfData.train[,CVariableSelection.ReduceModel.getMinModel(oVar.sub, 1)])
colnames(dfData) = CVariableSelection.ReduceModel.getMinModel(oVar.sub, 1)

head(dfData)
dfData$fGroups = fGroups
### fit a model to this data
fit.1 = lda(fGroups ~ ., data=dfData)

## check the prediction on this data
ivPredict = predict(fit.1, newdata = dfData)$posterior[, 'ATB']

## how good is this predictor
library(lattice)
densityplot(~ ivPredict, groups=fGroups, type='n', xlab='Predicted Score', main='Model predicted score', auto.key = list(columns=2))
xyplot(ivPredict ~ lData.train$grouping, xlab='Actual Group', ylab='Predicted Probability of Being ATB (1)',
       main='Predicted scores vs Actual groups')

### plot the score performance
cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, 1)
dfData.train = data.frame(t(lData.train$data))
dfData.train = data.frame(dfData.train[,cvTopGenes.sub])
colnames(dfData.train) = cvTopGenes.sub

oCV = CCrossValidation.LDA(test.dat = dfData.train, train.dat = dfData.train, test.groups = fGroups,
                           train.groups = fGroups, level.predict = 'ATB', boot.num = 50)

plot.cv.performance(oCV)

### check the cutoff tables
df = getCutoffTprFprCrossValidation(oCV)

fPredict = rep('reject', times=length(ivPredict))
fPredict[ivPredict >= 0.6] = 'ATB'
table(fPredict, fGroups)


## format the test data
fGroups.test = rep('ATB', times=length(lData.test$grouping))
fGroups.test[lData.test$grouping != 'ATB'] = 'Other'
fGroups.test = factor(fGroups.test, levels = c('Other', 'ATB'))
table(fGroups.test, lData.test$grouping)

## format the dataset
dfData.test = data.frame(t(lData.test$data))
dfData.test = data.frame(dfData.test[,cvTopGenes.sub])
colnames(dfData.test) = cvTopGenes.sub
head(dfData.test)

## does the batch show any difference here
xyplot(dfData.test$GBP6 ~ fGroups.test | lData.test$batch, xlab='Groups', ylab='Expression Level')

dfData.test$fGroups.test = fGroups.test
head(dfData.test)
### fit a model to this data
fit.2 = lda(fGroups.test ~ ., data=dfData.test)

## check the prediction on this data
ivPredict = predict(fit.2, newdata = dfData.test)$posterior[, 'ATB']

xyplot(ivPredict ~ lData.test$grouping, xlab='Actual Group', ylab='Predicted Probability of Being ATB (1)',
       main='Predicted scores vs Actual groups in Test data')

## predict in the test data
fPredict = rep('reject', times=length(ivPredict))
fPredict[ivPredict >= 0.6] = 'ATB'
table(fPredict, fGroups.test)

