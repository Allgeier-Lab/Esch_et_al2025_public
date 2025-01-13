##--------------------------------------------##
##    Katrina Sky Munsterman                  ##
##    kmunster@umich.edu                      ##
##    www.github.com/kmunsterman              ##
##    fish production                         ##
##--------------------------------------------##

## code adapted from Hamilton et al, 2022. 
## "Climate impacts alter fisheries productivity and turnover on coral reefs"

## load packages 
library(devtools)
devtools::install_github("renatoamorais/rfishprod")
library(rfishprod)
library(plyr)
library(tidyverse)
library(Rlab)


###import survey data with associated parameters from Fishbase or other sources###
###Example data###
comp<-read.csv("/Users/lanie/Documents/UM/AR_Prod-BioDiv/Example.prod.data.csv")


# case format
sum(comp$Count) 

fish <- comp[rep(seq_len(nrow(comp)), comp$Count), 1:12]
head(fish)

# check if any Lmeas are larger or equal to Linf
fish[fish$Lmeas>=fish$Linf,]
fish[fish$Lmeas2==fish$Lmax,]
fish[fish$Lmeas>=fish$Lmax,]

# change these individuals to Lmax - 0.1cm:
fish$Lmeas2 <- ifelse(fish$Lmeas==fish$Lmax | fish$Lmeas>fish$Lmax, fish$Lmax-0.1, fish$Lmeas)
#fish$Lmeas2 <- ifelse(fish$Lmeas==fish$Lmax, fish$Lmax-0.1, fish$Lmeas)
#fish$Lmeas2<-ifelse(fish$Lmeas>fish$Lmax, fish$Lmax-0.1, fish$Lmeas)

#### Calculate KMax, Operational Age, Length ####

# calculate KMax (Morais equation 6)
# FIX: check trait-based approach to calculate KMax Morais and Bellwood 2020
fish$index <- log10(fish$k) - (-2.31) * log10(fish$Linf)
fish$Kmax <- 10 ^ (fish$index + (-2.31) * log10(fish$Lmax))

# calculate age estimates
fish$EstAge <- (1/fish$Kmax)*log((fish$Lmax)/((1-fish$Lmeas2/fish$Lmax)*fish$Lmax))

# calculate growth with VBGF, using ages

# age to add for each day of the year:
age <- (1:365)/365

# create table for new length: 1 column per day, 1 row per individual fish
VB_lngth <-  matrix(ncol=length(age), nrow=nrow(fish), dimnames=list(NULL, paste("Day", 1:365, sep="_")))

# for each individual, calculate new length for each day using VBGF formula

for(i in 1:nrow(fish)) {
  VB_lngth[i, ] <- fish$Lmax[i]*(1-exp(-fish$Kmax[i]*(fish$EstAge[i] + age)))
}

VB_lngth

# matrix of lengths for each day of the year per fish

# convert lengths to weights using a & b coefficents:
VB_wt <-  apply(VB_lngth, 2, function(x) as.numeric(fish$a.final)*(x^as.numeric(fish$b.final)))   

# units = grams (per fish)

#### Predict Instantaneous Mortality ####

# predicting M: the instantaneous mortality rate

fish$Md <- with (fish,
                 rfishprod::predM (Lmeas = Lmeas2,
                                   t = 1,   
                                   Lmax = Lmax,
                                   Kmax = Kmax,
                                   method = 'Lorenzen'))  
# Md changes with individual length

#### Calculate Survival Prob ####

fish$SurvivalDay <- exp(-fish$Md)

# make a table of daily survival probabilities for the whole year:
prob.surv_tab <- matrix(ncol=365, nrow=nrow(fish), dimnames=list(NULL, paste("Day", 1:365)))

# FIX: length-based mortality instead of age-based mortality
for(i in 1:nrow(fish)) {
  prob.surv_tab [i, ] <- exp(-fish$Md[i] * ((fish$EstAge[i]+age)/(fish$EstAge[i]+1)))
}
#View(prob.surv_tab)

# proof of concept - cumulative survival 
prob.surv_cum.tab <- prob.surv_tab

for (i in 2:ncol(prob.surv_cum.tab))
  for (e in 1:nrow(prob.surv_cum.tab)) {
    prob.surv_cum.tab[e,i] <- prob.surv_cum.tab[e,i]*prob.surv_cum.tab[e,(i-1)]
  }

# probability of survival DECREASES through time, because this is now CUMULATIVE survival probability for each individual

##### Simulate Stochastic Mortality of each fish per day ####
# Bernoulli trials: run multiple times and take the mean

#### 1) Single Iteration of Survival ####
surv.outcome_day <- matrix(ncol=365, nrow=nrow(prob.surv_tab))

for (i in 1:nrow(prob.surv_tab))
  for (e in 1:ncol(prob.surv_tab)) {
    surv.outcome_day[i,e] <- Rlab::rbern(1,prob=prob.surv_tab[i,e])
  }

# 1 = survives to next day, 0 = fish dies   
# This is one iteration only

# create function: once a 0 appears in a row, all subsequent values are 0 (i.e. fish stays dead!)
zero_row_replace <- function(x) {
  x[which(cumany(x==0))] <- 0
  return(x)
}

# apply function
surv.outcome_day2 <- plyr::adply(surv.outcome_day, 1, zero_row_replace)
#View(surv.outcome_day2)

# remove first column, X1: (MAKE SURE TO ONLY RUN THIS ONCE!)
surv.outcome_day2 <- surv.outcome_day2[,-1]

# calculate fish mass from length using formula: W = a * L^b
head(fish)

biomass.g <- (as.numeric(fish$a.final) * fish$Lmeas2^as.numeric(fish$b.final))

# add mass as first column (Day 0) in matrix:
VB_wt2 <- cbind(biomass.g, VB_wt)
colnames(VB_wt2)[1] <- "Day_0"

# convert fish masses to mass gained per day (productivity) by subtracting previous day:

# create new object:
VB_prod_wt <- VB_wt2

for(i in 2:ncol(VB_prod_wt))
  for(e in 1:nrow(VB_prod_wt)) {
    VB_prod_wt[e,i] <- VB_wt2[e,i] - VB_wt2[e,(i-1)]
  }

# delete first column:
VB_prod_wt2 <- VB_prod_wt[,2:ncol(VB_prod_wt)]

# table of daily production in g/individual over a year (NOT CUMULATIVE)

# multiply survival table by production table to get production estimates:

# check matrices are the same size
dim(VB_prod_wt2)       
dim(surv.outcome_day2)

est_prod_day <- surv.outcome_day2*VB_prod_wt2
# all days where fish did not survive have 0 productivity

# sum daily production over the course of a year:
prod_yr <- rowSums(est_prod_day)
# 1 value per row (grams produced over 1 year)

prod_day <- est_prod_day$V1    #value on day one
# 1 value per row (i.e. grams produced after 1 day)

# combine daily and annual productivity  estimates
prod.per.fish <- cbind(fish,prod_day,prod_yr)

##### 2) Multiple Iterations of Survival ####

prob.surv_tab
y <-matrix(ncol=1, nrow= nrow(prob.surv_tab))

for (x in 1:100){   # FIX: check error to determine number of iterations to run
  
  surv_outcome_stoc <- matrix(ncol=365, nrow=nrow(prob.surv_tab))
  
  for (i in 1:nrow(prob.surv_tab))
    for (e in 1:ncol(prob.surv_tab)) {
      surv_outcome_stoc [i,e] <- Rlab::rbern(1,prob= prob.surv_tab[i,e])
    }
  
  surv_outcome_stoc2<-adply(surv_outcome_stoc, 1, zero_row_replace)
  surv_outcome_stoc2<-dplyr::select(surv_outcome_stoc2, -c(X1))
  
  est_daily_prod_tab<-surv_outcome_stoc2*VB_prod_wt2
  
  prod_year<-rowSums(est_daily_prod_tab)
  
  y   <- cbind(y,prod_year)
}

head(y)
# annual productivity per individual (each column = separate 1 year simulation)

# remove first column of NAs
y <- y[,-1]
head(y)

# get mean per fish across iterations:
mean.prod2 <- data.frame(prod_yr=rowSums(y)/100)
head(mean.prod)

#### Combine DFs, Save ####

# combine "fish" dataset with individual productivity estimates:
final.prod <- cbind(prod.per.fish, mean.prod2)

prod.sum<- final.prod %>% unite("reef.date", c(Reef,Date),  sep = "-") %>% group_by(reef.date) %>% summarise_at(vars(prod_yr), list(sum)) %>%
  filter(!reef.date %in% c("H5-6/5/18", "H6-6/5/18", "LH2-5/24/17", "LH3-5/21/18", "ES1-5/21/17", "MOW1-5/21/17", "MOW2-5/21/17", "MOW3-5/21/17", "MOW4-5/20/18", "MOW4-5/21/17") )

write.csv(prod.sum, "Comp_Prod.csv", row.names=F)

