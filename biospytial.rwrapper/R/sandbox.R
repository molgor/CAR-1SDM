# This is a script for loading all the necessary data for processing 
# i.e. it is a preprocessing script.
# Change it appropriately

my.binomial.bymCAR <- function(formula, data=NULL, trials, W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, MALA=TRUE, verbose=TRUE)


my.binomial.bymCAR(formula=formula_sample,W=M_bis,trials = trials,data=TDF,burnin=10000,n.sample=15000,verbose = TRUE)



    model.sample <-my_CARB(formula=formula_sample,family="binomial",W=M_bis,trials = trials,data=TDF,burnin=10000,n.sample=15000,verbose = TRUE,my_func=my.binomial.bymCAR)

    model.sample <-S.CARbym(formula=formula_sample,family="binomial",W=M_bis,trials = trials,data=TDF,burnin=10000,n.sample=15000,verbose = TRUE)
