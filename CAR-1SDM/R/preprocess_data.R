## [[file:preprocess_data.org::*Necessary libraries][Necessary libraries:1]]
    setwd("R")
    #library("CAR-1SDM")
## Necessary libraries:1 ends here

## [[file:preprocess_data.org::*Necessary libraries][Necessary libraries:2]]
    library(CARBayes)
## Necessary libraries:2 ends here

## [[file:preprocess_data.org::*Necessary libraries][Necessary libraries:3]]
    library(dplyr)
    library(purrr)
    library(reticulate)
## Necessary libraries:3 ends here

## [[file:preprocess_data.org::*Loading the /choosing principle/ functions][Loading the /choosing principle/ functions:1]]
source("ChoosingPrinciples.R")
## Loading the /choosing principle/ functions:1 ends here

## [[file:preprocess_data.org::*Read the adjancency matrix of the spatial lattice][Read the adjancency matrix of the spatial lattice:1]]
#file = '../data/presence_only_models/predictors/dataset100x100-puebla-p9/0-pred.csv'
#PDF = read.csv(file)
## REad adjancency matrix
# Import adjancency matrix generated from region
mat_filename = "../data/training_data_sample_puebla_p9_abies_pinophyta_adjmat.npy"
# Use numpy functions
np <- import("numpy")
M <- np$load(mat_filename)
## Read the adjancency matrix of the spatial lattice:1 ends here

## [[file:preprocess_data.org::*Read /training/ dataframe][Read /training/ dataframe:1]]
TDF = read.csv("../data/training_data_sample_puebla_p9_abies_pinophyta.csv")
#TDF = read.csv("../data/training_data_sample_puebla_p9_tyrannidae_birds.csv")
## Order it according to the id of the cell
TDF = TDF[order(TDF$cell_ids),]
# Convert the columns to numeric

## Beware older implementations of as.numeric function turns strings to (apparently) random numeric values
TDF = mutate_at(TDF,vars(Dist.to.road_m,Elevation_m,
                         MaxTemp_m,MeanTemp_m,
                         MinTemp_m,Population_m,
                         Precipitation_m,
                         SolarRadiation_m,
                         VaporPres_m,
                         WindSp_m),na_if,"N.A.")



TDF = mutate_at(TDF,vars(Dist.to.road_m,Elevation_m,
                         MaxTemp_m,MeanTemp_m,
                         MinTemp_m,Population_m,
                         Precipitation_m,
                         SolarRadiation_m,
                         VaporPres_m,
                         WindSp_m),as.numeric)

# Remove unnecessary symbols in variable names
names(TDF) = lapply(names(TDF),function(x) gsub("_","",x))
names(TDF) = lapply(names(TDF),function(x) gsub("\\.","",x))
## Read /training/ dataframe:1 ends here

## [[file:preprocess_data.org::*Generating relative absences with the /choosing principle/][Generating relative absences with  the /choosing principle/:1]]
# Change the name of a column that for some reason is called the same
names(TDF)[23] <- 'code_id2'

## Treatment I, missing values in X and Y, comment this if using treatment II
DataFrame = TDF %>% rowwise() %>%
              mutate(sample=pseudo_absence_naive(Plantae,LUCA),
                     species=pseudo_absence_naive(Pinophyta,Plantae))

## Uncomment this for treatment II (i.e. missing values only in X)
#DataFrame = TDF %>% rowwise() %>%
#            mutate(sample=pseudo_absence_naive(Plantae,LUCA),
#                   species=pseudo_absence_trivial(Pinophyta,Plantae))

## Uncomment this if you want to assume that all missing data are absences
## i.e. remove NAs
#
## Remove entries in the adjacency matrix that correspond to missing data
#rr <- DataFrame %>%
#    filter(!is.na(species) & !is.na(sample))
#
#sam_idx_nan <- which(is.na(DataFrame$sample))
#
#M = M[-c(sam_idx_nan),-c(sam_idx_nan)]

## Remove missig data in DataFrame
#DataFrame = TDF %>% rowwise() %>%
#  mutate(sample=pseudo_absence_trivial(Plantae,LUCA),
#         species=pseudo_absence_trivial(Pinophyta,Plantae))
## Generating relative absences with  the /choosing principle/:1 ends here

## [[file:preprocess_data.org::*Preprocess adjancency matrix $M$][Preprocess adjancency matrix $M$:1]]
## Remove entries with zero neighbours (adjacency matrix)
### Calculates number of neighbours in D (sum)
D = apply(M,MARGIN = 1,sum)
### get index with 0 neighbours
idx = which(D == 0)

### select cells with no neighbours
cell_with_no_neighbour = TDF$cellids[idx]
## Remove island for TDF
TDF <- TDF[-c(idx),]
## Erase idx for M and for TDF
M_bis = M[-c(idx),-c(idx)]

## remove rows that have no neighbours (islands)
DataFrame <- DataFrame[-c(idx),]

n <- dim(M_bis)[1]
trials <- rep(1,n)
## Preprocess adjancency matrix $M$:1 ends here

## [[file:preprocess_data.org::*Preprocess adjancency matrix $M$][Preprocess adjancency matrix $M$:2]]
## Replace missing values with mean, Of course we could do this using other more fancy method
covs2work = c("Disttoroadm","Populationm","Elevationm","Precipitationm","MeanTempm")
for(i in covs2work){
    DataFrame[,i][is.na(DataFrame[,i])] <- mean(DataFrame[,i][!is.na(DataFrame[,i])],na.rm=TRUE)
}
## Preprocess adjancency matrix $M$:2 ends here
