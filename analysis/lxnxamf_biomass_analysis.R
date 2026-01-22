# lxnxamf_biomass_analysis.R
## script to analyze biomass allocation data from light x nitrogen x amf experiment
## full experiment led by Snehanjana Chatterjee
## this section is being led by Rea Farrow
## author: Nick Smith

## load packages
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(mvtnorm)

## load multcompView## load data
biomass_data <- read.csv('../data/lxnxamf_biomass.csv')
head(biomass_data)

## calculate metrics
biomass_data$leaf_weight <- biomass_data$focal_leaf_weight + 
  biomass_data$other_leaf_weight # total leaf biomass
biomass_data$total_weight <- biomass_data$leaf_weight + 
  biomass_data$stem_weight + biomass_data$root_weight # total plant biomass
biomass_data$leaf_allocation <- biomass_data$leaf_weight / biomass_data$total_weight # relative allocation to leaves
biomass_data$stem_allocation <- biomass_data$stem_weight / biomass_data$total_weight # relative allocation to stems
biomass_data$root_allocation <- biomass_data$root_weight / biomass_data$total_weight # relative allocation to roots

## analyze data

### total_weight
### hypothesis: weight should be greater with increased light, soil N, and in presence of AMF
#### check data distribution
hist(biomass_data$total_weight) # very skewed, is it better log transformed?
hist(log(biomass_data$total_weight)) # yes, much better. let's go with this
#### fit linear model
total_weight_lm <- lm(log(total_weight) ~ light_treatment * nitrogen_treatment * amf_treatment, data = biomass_data)
plot(resid(total_weight_lm) ~ fitted(total_weight_lm)) # look okay
summary(total_weight_lm)
Anova(total_weight_lm) # light treatment effect, N treatment effect
#### explore significant effects
emmeans(total_weight_lm, ~ light_treatment) # high light = more biomass
emmeans(total_weight_lm, ~ nitrogen_treatment) # high N = more biomass




