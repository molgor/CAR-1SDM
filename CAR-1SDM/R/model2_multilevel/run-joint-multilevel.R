# This is the script that implements the Model 2 (common GMRF)
## Cross validation
library(pROC)
library('caret')

## Set Working directory
## Import code:
#setwd('/app/external_plugins/biospytial_rwrapper/CAR-1SDM')
## remove any previous object
#rm(list=ls())

## Init data source with variables
#print('Load data source and preprocess')

formula_sample =  sample ~ Disttoroadm + Populationm
formula_presence = species ~ Elevationm + Precipitationm



## Build dataframes, S<- Sample, P<- Presence
S <- model.frame(formula_sample, DataFrame,na.action='na.pass')
P <- model.frame(formula_presence, DataFrame,na.action='na.pass')

## Split lhs and rhs from the design matrix
SX <- select(S, -c(1))
PX <- select(P, -c(1))
Sy <- select(S, c(1))
Py <- select(P, c(1))

names(Sy) <- 'response'
names(Py) <- names(Sy)
## Stack both processes into same dataframe
## First the responses (Y) will be concatenated by row.
Y = rbind(Sy,Py)

###
T1 <- matrix(rep(0,4), ncol = 2)
T2 <- matrix(rep(0,4), ncol = 2)
T1[1,1] <- 1
T2[2,2] <- 1

## Perform Kronnecker with different covariates (Block diagonal)
X <- data.frame((T1 %x% as.matrix(SX)) + (T2 %x% as.matrix(PX)))
names(X) <- c(names(SX),names(PX))

DD <- cbind(Y,X)

nK <- dim(M_bis)[1]
## make sequence vector for id.area
ida <- data.frame(seq(nK))
idarea <- unlist(rbind(ida,ida))
## make sequence vector for correlation 
corx <- rep(x = 1,times = nK)
cory <- rep(x = 2,times = nK)

indre <- c(corx,cory)
## A general formula for all the covariates
#formula <- response ~ Disttoroadm + Populationm + Elevationm + MeanTempm

formula <- response ~ Disttoroadm + Populationm + Elevationm + Precipitationm
###### Runnning the model
## now, assuming that the order in M_bis is the same as in cellids (OOOORDEEER, not value)
## Run the model
trials = rep(1,2 * nK)
burnin = 50000
n.sample = 100000
thin = 50




#trains = createFolds(y = DataFrame$species, k=7, returnTrain = TRUE)
#validate = createFolds(y = DataFrame$species, k=7, returnTrain = FALSE)

#DataFrame$presences <- DataFrame$species
#model2 <- S.CARmultilevel(formula,family = 'binomial',
#                          trials=trials, 
#                          W=M_bis, 
#                          ind.area = idarea,
#                         ind.re=factor(idarea),
#                          rho = 1,
#                          burnin = burnin,
#                          n.sample = n.sample,
#                          data = DD
#                         )
#
#

