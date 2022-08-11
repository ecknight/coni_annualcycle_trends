library(tidyverse)
library(usdm)
library(lme4)

#Adapted from https://github.com/crushing05/WOTH_retro/blob/master/single_pop_gvs.R

#Thought: Don't need to do per population because defined regions of connectivity already
#So how to reduce the models?
#And should I use the jags model or should I just use glmer?
#Should think about what I want to report
#Consider going back and adding p(use) to regions for covariate calculation via KDE and gaussian filter on the buffers.
#Make sure to reiterate the assumptions statement in Mike's paper

#1. Read in data & scale covs----
dat <- read.csv("Data/MonitoringData_Model.csv") %>% 
  mutate(years = (year-2002)/17,
         pdsis_breed = scale(pdsidt_breed),
         pdsis_fall = scale(pdsidt_fall),
         pdsis_winter = scale(pdsidt_winter),
         pdsis_spring = scale(pdsidt_spring),
         crops_breed = scale(cropdt_breed),
         crops_fall = scale(cropdt_fall),
         crops_winter = scale(cropdt_winter),
         crops_spring = scale(cropdt_spring),
         trees_breed = scale(treedt_breed),
         trees_fall = scale(treedt_fall),
         trees_winter = scale(treedt_winter),
         trees_spring = scale(treedt_spring))

#2. Check for covariance----
covs <- dat %>% 
  dplyr::select(pdsidt_breed:cropdt_winter)

cor(covs, use="na.or.complete")
#tree_breed and tree_winter are the only ones

vif(covs)
#I think leave it, everything is under 5

#2. MLE tests for univariate significance----

dat.pop <- dat %>% 
  group_by(method, pop, year) %>% 
  summarize(abun = mean(count),
            off = mean(p)) %>% 
  ungroup() %>% 
  left_join(dat %>% 
              dplyr::select(method, pop, year, pdsidt_breed:cropdt_winter) %>% 
              unique())
 
ggplot(pop_abun) +
  geom_histogram(aes(x=count))

pops <- unique(pop_abun$pop)

for(i in 1:length(pops)){
  
  dat.i <- dat %>% 
    dplyr::filter(pop==pops[i])
  
  fit1 <- glm(abun ~ year + pdsidt_breed, offset = off, family="quasipoisson", data=dat.i)
  fit1 <- glmer(count ~ years + pdsis_breed  + (1|route), offset = p, family="poisson", data=dat.i)
  fit2 <- glm(abun ~ year + pdsidt_winter, offset = off, family="quasipoisson", data=dat.i)
  fit3 <- glm(abun ~ year + pdsidt_spring, offset = off, family="quasipoisson", data=dat.i)
  fit4 <- glm(abun ~ year + cropdt_breed, offset = off, family="quasipoisson", data=dat.i)
  fit5 <- glm(abun ~ year + cropdt_winter, offset = off, family="quasipoisson", data=dat.i)
  fit6 <- glm(abun ~ year + cropdt_spring, offset = off, family="quasipoisson", data=dat.i)
  fit7 <- glm(abun ~ year + treedt_breed, offset = off, family="quasipoisson", data=dat.i)
  fit8 <- glm(abun ~ year + treedt_winter, offset = off, family="quasipoisson", data=dat.i)
  fit9 <- glm(abun ~ year + treedt_spring, offset = off, family="quasipoisson", data=dat.i)
  
  if(pops[i]!=16){
    fit10 <- glm(abun ~ year + pdsidt_fall, offset = off, family="quasipoisson", data=dat.i)
    fit11 <- glm(abun ~ year + cropdt_fall, offset = off, family="quasipoisson", data=dat.i)
    fit12 <- glm(abun ~ year + treedt_fall, offset = off, family="quasipoisson", data=dat.i)
    
    betas <- c(coef(fit1)[3],coef(fit2)[3],coef(fit3)[3],coef(fit4)[3],coef(fit3)[3],coef(fit6)[3],coef(fit7)[3],
               coef(fit8)[3],coef(fit9)[3],coef(fit10)[3],coef(fit11)[3],coef(fit12)[3])
    
    pop.lcl <- c(confint(fit1,level = .8)[3,1],confint(fit2,level = .8)[3,1],confint(fit3,level = .8)[3,1],
                 confint(fit4,level = .8)[3,1],confint(fit3,level = .8)[3,1],confint(fit6,level = .8)[3,1],
                 confint(fit7,level = .8)[3,1],confint(fit8,level = .8)[3,1],confint(fit9,level = .8)[3,1],
                 confint(fit10,level = .8)[3,1],confint(fit11,level = .8)[3,1],confint(fit12,level = .8)[3,1])
    
    pop.ucl <- c(confint(fit1,level = .8)[3,2],confint(fit2,level = .8)[3,2],confint(fit3,level = .8)[3,2],
                 confint(fit4,level = .8)[3,2],confint(fit3,level = .8)[3,2],confint(fit6,level = .8)[3,2],
                 confint(fit7,level = .8)[3,2],confint(fit8,level = .8)[3,2],confint(fit9,level = .8)[3,2],
                 confint(fit10,level = .8)[3,2],confint(fit11,level = .8)[3,2],confint(fit12,level = .8)[3,2])
    
    temp.df <- data.frame(pop=rep(pops[i],12), predictor = seq(1:12), est = betas, lcl=pop.lcl, ucl=pop.ucl)
    ifelse(i==1, pop.df <- temp.df, pop.df <- rbind(pop.df,temp.df))
    
  }
  else{
    
    betas <- c(coef(fit1)[3],coef(fit2)[3],coef(fit3)[3],coef(fit4)[3],coef(fit3)[3],coef(fit6)[3],coef(fit7)[3],
               coef(fit8)[3],coef(fit9)[3])
    
    pop.lcl <- c(confint(fit1,level = .8)[3,1],confint(fit2,level = .8)[3,1],confint(fit3,level = .8)[3,1],
                 confint(fit4,level = .8)[3,1],confint(fit3,level = .8)[3,1],confint(fit6,level = .8)[3,1],
                 confint(fit7,level = .8)[3,1],confint(fit8,level = .8)[3,1],confint(fit9,level = .8)[3,1])
    
    pop.ucl <- c(confint(fit1,level = .8)[3,2],confint(fit2,level = .8)[3,2],confint(fit3,level = .8)[3,2],
                 confint(fit4,level = .8)[3,2],confint(fit3,level = .8)[3,2],confint(fit6,level = .8)[3,2],
                 confint(fit7,level = .8)[3,2],confint(fit8,level = .8)[3,2],confint(fit9,level = .8)[3,2])
    
    temp.df <- data.frame(pop=rep(pops[i],9), predictor = seq(1:9), est = betas, lcl=pop.lcl, ucl=pop.ucl)
    ifelse(i==1, pop.df <- temp.df, pop.df <- rbind(pop.df,temp.df))
  }

}

pop.df$impt <- ifelse(pop.df$est > 0 & pop.df$lcl > 0, 1,
                      ifelse(pop.df$est < 0 & pop.df$ucl < 0, 1, 0))

View(pop.df)

#3. Try with a mixed model----
mod.all <- glmer(count ~ years + 
               pdsis_breed + pdsis_fall + pdsis_winter + pdsis_spring +
               crops_breed + crops_fall + crops_winter + crops_spring +
               trees_breed + trees_fall + trees_winter + trees_spring +
               (1|pop/route), offset=p, data=dat, family="poisson")

int.all <- confint.merMod(mod.all, level = 0.95, devtol = 1e-7)

#4. Compare to BBS only----
dat.bbs <- dat %>% 
  dplyr::filter(method=="BBS")

mod.bbs <- glmer(count ~ years + 
               pdsis_breed + pdsis_fall + pdsis_winter + pdsis_spring +
               crops_breed + crops_fall + crops_winter + crops_spring +
               trees_breed + trees_fall + trees_winter + trees_spring +
               (1|pop/route), offset=p, data=dat, family="poisson")

int.bbs <- confint.merMod(mod.bbs, level = 0.95)

mod.all <- mod
int.all <- mod.int
rm(mod)
rm(mod.int)
