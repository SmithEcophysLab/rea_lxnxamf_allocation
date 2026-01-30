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

## load data
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
biomass_data$aboveground_allocation <- biomass_data$leaf_allocation + biomass_data$stem_allocation
biomass_data$rootshoot <- biomass_data$root_weight / (biomass_data$leaf_weight + 
                                                        biomass_data$stem_weight)

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

### leaf_allocation
### hypothesis: leaf allocation should be greater with increased soil N and in presence of AMF, but reduced with increased light
#### check data distribution
hist(biomass_data$leaf_allocation) # ok!
#### fit linear model
leaf_allocation_lm <- lm(leaf_allocation ~ light_treatment * nitrogen_treatment * amf_treatment, data = biomass_data)
plot(resid(leaf_allocation_lm) ~ fitted(leaf_allocation_lm)) # look okay
summary(leaf_allocation_lm)
Anova(leaf_allocation_lm) # light treatment effect, N treatment effect
#### explore significant effects
emmeans(leaf_allocation_lm, ~ light_treatment) # high light = less leaf allocation
emmeans(leaf_allocation_lm, ~ nitrogen_treatment) # high N = more leaf allocation

### stem_allocation
### hypothesis: stem allocation should be greater with increased soil N and in presence of AMF, but reduced with increased light
#### check data distribution
hist(biomass_data$stem_allocation) # ok!
#### fit linear model
stem_allocation_lm <- lm(stem_allocation ~ light_treatment * nitrogen_treatment * amf_treatment, data = biomass_data)
plot(resid(stem_allocation_lm) ~ fitted(stem_allocation_lm)) # look okay
summary(stem_allocation_lm)
Anova(stem_allocation_lm) # N treatment effect, AMF effect, NxAMF interaction
#### explore significant effects
emmeans(stem_allocation_lm, ~ nitrogen_treatment) # high N = more stem allocation
emmeans(stem_allocation_lm, ~ amf_treatment) # yes AMF = less stem allocation
emmeans(stem_allocation_lm, ~ nitrogen_treatment * amf_treatment) # lowest stem allocation in low N and yes AMF

### aboveground_allocation
### hypothesis: AG allocation should be greater with increased soil N and in presence of AMF, but reduced with increased light
#### check data distribution
hist(biomass_data$aboveground_allocation) # skewed, let's try log transformation...
hist(log(biomass_data$aboveground_allocation)) # better
#### fit linear model
aboveground_allocation_lm <- lm(log(aboveground_allocation) ~ light_treatment * nitrogen_treatment * amf_treatment, data = biomass_data)
plot(resid(aboveground_allocation_lm) ~ fitted(aboveground_allocation_lm)) # look okay
summary(aboveground_allocation_lm)
Anova(aboveground_allocation_lm) # light treatment effect, N treatment effect, AMF effect
#### explore significant effects
emmeans(aboveground_allocation_lm, ~ light_treatment) # high light = less AG allocation
emmeans(aboveground_allocation_lm, ~ nitrogen_treatment) # high N = more AG allocation
emmeans(aboveground_allocation_lm, ~ amf_treatment) # yes AMF = less AG allocation

### root_allocation
### hypothesis: root allocation should be greater with increased light, but reduced with soil N and in presence of AMF
#### check data distribution
hist(biomass_data$root_allocation) # skewed, let's try log transformation...
hist(log(biomass_data$root_allocation)) # better
#### fit linear model
root_allocation_lm <- lm(log(root_allocation) ~ light_treatment * nitrogen_treatment * amf_treatment, data = biomass_data)
plot(resid(root_allocation_lm) ~ fitted(root_allocation_lm)) # look okay
summary(root_allocation_lm)
Anova(root_allocation_lm) # light treatment effect, N treatment effect, AMF effect
#### explore significant effects
emmeans(root_allocation_lm, ~ light_treatment) # high light = more root allocation
emmeans(root_allocation_lm, ~ nitrogen_treatment) # high N = less root allocation
emmeans(root_allocation_lm, ~ amf_treatment) # yes AMF = more root allocation

### rootshoot
### hypothesis: root:shoot ratio should be greater with increased light, but reduced with soil N and in presence of AMF
#### check data distribution
hist(biomass_data$rootshoot) # skewed, let's try log transformation...
hist(log(biomass_data$rootshoot)) # better
#### fit linear model
rootshoot_lm <- lm(log(rootshoot) ~ light_treatment * nitrogen_treatment * amf_treatment, data = biomass_data)
plot(resid(rootshoot_lm) ~ fitted(rootshoot_lm)) # look okay
summary(rootshoot_lm)
Anova(rootshoot_lm) # light treatment effect, N treatment effect, AMF effect
#### explore significant effects
emmeans(rootshoot_lm, ~ light_treatment) # high light = greater root:shoot
emmeans(rootshoot_lm, ~ nitrogen_treatment) # high N = less root:shoot
emmeans(rootshoot_lm, ~ amf_treatment) # yes AMF = greater root:shoot

## plot data
### rootshoot
rootshoot_plot <- ggplot(aes(x = amf_treatment, y = rootshoot, alpha = light_treatment, fill = nitrogen_treatment),
                         data = biomass_data) +
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_boxplot() +
  scale_fill_manual(values = c('brown', 'darkgreen')) +
  scale_alpha_manual(values = c(1, 0.5)) +
  labs(alpha = 'Light treatment') +
  labs(fill = 'Nitrogen treatment') +
  xlab('AMF Treatment') +
  ylab('Root:Shoot') +
  ylim(c(0, 2))

jpeg(filename = "../results/plots/rea_rootshoot.jpeg", 
     width = 8, height = 8, units = 'in', res = 300)
plot(rootshoot_plot)
dev.off()
