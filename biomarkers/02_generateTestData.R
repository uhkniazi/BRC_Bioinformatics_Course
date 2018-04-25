# Name: 02_generateTestData.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 25/04/2018
# Desc: generate test data

### libraries to load
library(GEOquery)
library(Biobase)
library(lumi)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)
library(annotate)
library(downloader)

setwd('biomarkers/')
dir.create('dataExternal')

## internal functions 
# Function: f_lGetPCAClusterCount
# Desc: Takes first two components of the PCA and counts possible clusters. The function does this by 
#       binning the vector of data into X bins, assigning a class label to each bin, counting how many
#       observations in each bin, any total number of bins with at least one observations. This is calculated
#       for both the components of the pca matrix, and the max number of bins with at least one observation, in
#       first or second dimension is reported along with a data.frame with cluster labels.
# Args: pr.out = principal component object returned by prcomp function
# Rets: returns list with 2 elements: 
#       1 - cluster.count = possible number of clusters in the data
#       2 - cluster.label = data.frame with cluster labels
f_lGetPCAClusterCount = function(pr.out){
  # how many clusters in data, using first 2 components
  x1 = pr.out$x[,1]
  x2 = pr.out$x[,2]
  # bin the data from the 2 components
  h1 = hist(x1, plot=F)
  # give a class label to each bin
  c1 = cut(x1, h1$breaks, labels = 1:(length(h1$mids)))
  h2 = hist(x2, plot=F)
  c2 = cut(x2, h2$breaks, labels = 1:(length(h2$mids)))
  # labels for vectors and the class labels
  dfClust = data.frame(lab=names(x1), c1, c2)
  # get contingency table
  mClust = as.matrix(table(c1 = dfClust$c1, c2 = dfClust$c2))
  # count the max of row and col sums that are not zero
  ir = length(which(rowSums(mClust) != 0))
  ic = length(which(colSums(mClust) != 0))
  iClust.count = ifelse(ir > ic, ir, ic)
  lRet = list(cluster.count=iClust.count, cluster.label=dfClust)
  return(lRet)
}

## download the data set
## go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19491

# url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19491/matrix/GSE19491_series_matrix.txt.gz'
# download(url, destfile='dataExternal/GSE19491_series_matrix.txt.gz')
oExp = getGEO(filename = 'dataExternal/GSE19491_series_matrix.txt.gz')
## ignore the warning as it is a bug
# add lumi nuIDs - converting probe ids to identify genes
oExp = addNuID2lumi(oExp, lib.mapping = 'lumiHumanIDMapping' )

# remove any NA data
i = which(is.na(rowSums(exprs(oExp))))
oExp = oExp[-i,]

## check the metadata
df = pData(oExp)
str(df)
head(df)

# print samples
as.data.frame(table(oExp$source_name_ch1))
as.data.frame(table(oExp$characteristics_ch1.3))

# get the whole blood data
i = grep('^Whole blood from healthy control$|^Whole Blood from patient with Active TB$|^Whole Blood from patient with Latent TB$', 
         x = oExp$source_name_ch1, ignore.case = T, perl = T)
oExp = oExp[,i]

# sanity check
as.data.frame(table(oExp$source_name_ch1))
as.data.frame(table(oExp$characteristics_ch1.3))

## create short names for these factors
df = pData(oExp)
fSamples = rep(NA, times=nrow(df))
f = as.character(oExp$source_name_ch1)
# create factors
i = grep('healthy control', f)
fSamples[i] = 'HC'

i = grep('Active TB', f)
fSamples[i] = 'ATB'

i = grep('Latent TB', f)
fSamples[i] = 'LTB'

fSamples = factor(fSamples, levels = c('HC', 'LTB', 'ATB'))
table(fSamples)
levels(fSamples)

df$fSamples = fSamples
df = droplevels.data.frame(df)
## update the table
pData(oExp) = df

#### examine the data matrix
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

mData = exprs(oExp)
range(mData)
dim(mData)
i = sample(1:nrow(mData), 2000, replace = F)
mData = mData[i,]
oDiag.1 = CDiagnosticPlots(mData, 'Raw Matrix')
fBatch = df$fSamples
## check data quality
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.PCA(oDiag.1, fBatch, cex.main=1, csLabels = '')
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.5, cex.main=0.7)

#################################### data normalization
# remove negative values first and set minimum value to 1
exprs(oExp) = exprs(oExp) + abs(min(exprs(oExp))) + 1
range(exprs(oExp))

mData = log(exprs(oExp))
range(mData)
dim(mData)
i = sample(1:nrow(mData), 2000, replace = F)
mData = mData[i,]
## data has negative values, make positive before further analysis
oDiag.1 = CDiagnosticPlots(mData, 'Log Raw Matrix')
fBatch = df$fSamples

## check normalisation
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.PCA(oDiag.1, fBatch, cex.main=1, csLabels = '')
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.5, cex.main=0.7)

# normalize and log2 transform the data using lumi
oExp.lumi = lumiT(oExp, 'log2')
## data normalization
oExp = lumiN(oExp.lumi, method='rsn')
rm(oExp.lumi)
gc()

## check data quality after normalisation
mData = exprs(oExp)
range(mData)
dim(mData)
set.seed(123)
i = sample(1:nrow(mData), 2000, replace = F)
mData = mData[i,]
## data has negative values, make positive before further analysis
oDiag.2 = CDiagnosticPlots(mData, 'Normalised Matrix')
fBatch = df$fSamples

## check normalisation
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)
plot.PCA(oDiag.2, fBatch, cex.main=1, csLabels = '')
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.5, cex.main=0.7)
## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.2)
l$PCA.jitter = F
l$HC.jitter = F
## this should give an error if scaling can't be done
## if all the vector 0 for PCA
oDiag.2.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)
plot.PCA(oDiag.2.2, fBatch, legend.pos = 'topright', csLabels = '')

l = f_lGetPCAClusterCount(oDiag.2.2@lData$PCA)
l$cluster.count
table(z2 = l$cluster.label$c2, z1 = l$cluster.label$c1)
fBatch = factor(ifelse(l$cluster.label$c2 %in% as.character(1:5), 'Batch1', 'Batch2'))

# sanity check
plot.PCA(oDiag.2.2, fBatch, legend.pos = 'top', csLabels = '')

oExp$fBatch = fBatch
xtabs(~ (oExp$characteristics_ch1.4) + oExp$fBatch + oExp$fSamples)

plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

# which part of experimental design is causing this?
df = pData(oExp)
f = grepl('PerfectPure', df$extract_protocol_ch1)
plot.PCA(oDiag.2.2, factor(as.numeric(f)), legend.pos = 'right', csLabels = '')

f = as.character(oExp$extract_protocol_ch1)
f[grepl('PerfectPure', df$extract_protocol_ch1)] = 'PerfectPure'
f[!grepl('PerfectPure', df$extract_protocol_ch1)] = 'MagMAX'

table(oExp$fBatch, f)

table(oExp$fBatch, oExp$fSamples)
rm(mData)
## assign symbols to genes
mDat = exprs(oExp)
df = select(lumiHumanAll.db, keys = rownames(mDat), columns=c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
## drop columns with NA
df = na.omit(df)
## drop columns with duplicate genes
i = duplicated(df$SYMBOL)
table(i)
df = df[!i, ]
head(df)
## put both tables in same order
i = match(df$PROBEID, rownames(mDat))
mDat = mDat[i,]
# sanity check
identical(rownames(mDat), df$PROBEID)
rownames(mDat) = df$SYMBOL

## load the training data and select the matching genes
load('results/lData.train.rds')

# how many genes match
table(rownames(lData.train$data) %in% rownames(mDat))
mDat = mDat[rownames(mDat) %in% rownames(lData.train$data), ]
dim(mDat)

lData.test = list(data=mDat, grouping=oExp$fSamples, batch=oExp$fBatch)
save(lData.test, file='results/lData.test.rds')
