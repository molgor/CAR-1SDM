x <- 1:5
square <- function(x) {
  x^2
}
square(x)

X <- 40:30
square(X)

source("../R/preprocess_data.R")
# load the building function
source("../R/joint.binomial.bymCARModel1.R")

ntot = length(DataFrame$sample)
npresence_1 = length(na.omit(DataFrame$species[DataFrame$species == 1])) / ntot*100
nsample_1 = length(na.omit(DataFrame$sample[DataFrame$sample == 1]))/ntot*100
npresence_0 = length(na.omit(DataFrame$species[DataFrame$species == 0]))/ntot*100
nsample_0 = length(na.omit(DataFrame$sample[DataFrame$sample == 0]))/ntot*100
n_miss_presence = sum(is.na(DataFrame$species))/ntot*100
n_miss_sample = sum(is.na(DataFrame$sample))/ntot*100
cbind(npresence_0,npresence_1,nsample_0,nsample_1,n_miss_presence,n_miss_sample)

## MCMC parameters burnin and sample to be consistent with the other models
n.sample = 1000
thin = 1
verbose = TRUE

## Make ROC curve
library(pROC)
library('caret')
#trains = createFolds(y = DataFrame$species, k=7, returnTrain = TRUE)

nonas = which(! is.na(DataFrame$species) )
Y_withoutNA = DataFrame$species[nonas]
validate = createFolds(y = Y_withoutNA, k=7, returnTrain = FALSE)


DataFrame$presences <- DataFrame$species

l <- list()
i = 1
DataFrame$predicted_values_CV <- NA
DataFrame$predicted_valuesBernoulli <- NA

formula_sample =  sample ~ Disttoroadm + Populationm
formula_presence = species ~ Elevationm + Precipitationm

for (fold in validate) {

    observed.presences <- DataFrame$species[fold]
    ## Substitue by NA
    DataFrame$species[fold] <- NA
    results  <- joint.binomial.bymCARModel1(formula_S = formula_sample,
                                        formula_P = formula_presence,
                                        n.sample=n.sample,
                                        data = DataFrame,
                                        burnin=burnin,
                                        postburnin=postburnin,
                                        thin=thin,
                                        verbose=TRUE)

    DataFrame$species <- DataFrame$presences
    ## Return original values
    DataFrame$predicted_values_CV[fold] <- results$fitted.values[fold]
    predicted.probability = results$fitted.values[fold]

    ## Generate Bernoulli sample [ Only for the fold data set]
    print("Generating bernoulli sampling...")

    post.joint = data.frame(results$samples$fitted.joint[fold])

    ptot <- post.joint %>% mutate_all(function(p) rbernoulli(1,p))

    sumpt <- colSums(ptot)

    nsamples = dim(ptot)[1]


    ProbPS <- sumpt / nsamples

    DataFrame$predicted_valuesBernoulli[fold] <- ProbPS


    pROC_obj <- roc(observed.presences,predicted.probability)

    l[[i]] <- pROC_obj
    i = i + 1
}

library(lattice)
xyplot(1:10 ~ 1:10)
