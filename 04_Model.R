library(tidyverse)
library(usdm)
library(lme4)
library(R2jags)
library(sf)
library(performance)

options(scipen = 99999)

#Adapted from https://github.com/crushing05/WOTH_retro/blob/master/single_pop_gvs.R

#Thought: Don't need to do per population because defined regions of connectivity already
#So how to reduce the models?
#And should I use the jags model or should I just use glmer?
#Should think about what I want to report
#Consider going back and adding p(use) to regions for covariate calculation via KDE and gaussian filter on the buffers.
#Make sure to reiterate the assumptions statement in Mike's paper when writing 

#1. Read in data & scale covs----
dat <- read.csv("Data/MonitoringData_Model.csv") %>% 
  dplyr::filter(year!=2002) %>% 
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
  dplyr::select(pdsis_breed:trees_spring)

cor(covs, use="na.or.complete")
#tree_breed and tree_spring are 0.83

vif(covs)
#Should take out tree_breed

covs <- covs %>% 
  dplyr::select(-trees_breed)

vif(covs)

#2. MLE tests for univariate significance----
pops <- unique(dat$pop)
mod.df <- data.frame()
od.df <- data.frame()
for(i in 1:length(pops)){
  
  dat.i <- dat %>% 
    dplyr::filter(pop==pops[i])
  
  fit1 <- glmer.nb(count ~ years + pdsis_breed  + (1|route), offset = p,  data=dat.i)
  fit2 <- glmer.nb(count ~ years + pdsis_winter  + (1|route), offset = p,  data=dat.i)
  fit3 <- glmer.nb(count ~ years + pdsis_spring  + (1|route), offset = p,  data=dat.i)
  fit4 <- glmer.nb(count ~ years + crops_breed  + (1|route), offset = p,  data=dat.i)
  fit5 <- glmer.nb(count ~ years + crops_winter  + (1|route), offset = p,  data=dat.i)
  fit6 <- glmer.nb(count ~ years + crops_spring  + (1|route), offset = p,  data=dat.i)
  fit7 <- glmer.nb(count ~ years + trees_breed  + (1|route), offset = p,  data=dat.i)
  fit8 <- glmer.nb(count ~ years + trees_winter  + (1|route), offset = p,  data=dat.i)
  fit9 <- glmer.nb(count ~ years + trees_spring  + (1|route), offset = p,  data=dat.i)
  
  mods <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9)
  
  #set up if else because no fall points for population 16
  if(pops[i]!=16){
    fit10 <- glmer.nb(count ~ years + pdsis_fall  + (1|route), offset = p,  data=dat.i)
    fit11 <- glmer.nb(count ~ years + crops_fall  + (1|route), offset = p,  data=dat.i)
    fit12 <- glmer.nb(count ~ years + trees_fall  + (1|route), offset = p,  data=dat.i)
    
    mods <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12)
  }
  
  for(j in 1:length(mods)){
    od.df <- check_overdispersion(mods[[j]]) %>% 
      rbind(od.df)
    }
  
  for(j in 1:length(mods)){
    
    mod.df <- data.frame(summary(mods[[j]])$coefficients) %>% 
      mutate(pop = pops[i],
             cov = row.names(summary(mods[[j]])$coefficients)) %>% 
      dplyr::filter(!cov %in% c("(Intercept)", "years")) %>% 
      rbind(mod.df)
    
    # cl <- try(confint.merMod(mods[[j]], method="profile", level=0.8, devtol = 1e-5))
    # if(class(cl)=="matrix" & !is.na(mods[[j]]@beta[3])){
    #   mod.df <- data.frame(beta = mods[[j]]@beta[3],
    #                       lcl = cl[4,1],
    #                       ucl = cl[4,2],
    #                       pop = pops[i],
    #                       cov = names(mods[[j]]@frame)[3]) %>% 
    #     rbind(mod.df)
    # } else {
    #   mod.df <- data.frame(beta = mods[[j]]@beta[3],
    #                        lcl = NA,
    #                        ucl = NA,
    #                        pop = pops[i],
    #                        cov = names(mods[[j]]@frame)[3]) %>% 
    #     rbind(mod.df)
    # }
  }
  
  print(paste0("Finished population ", i, " of ", length(pops)))
}

#3. Calculate indicator values----
# mod.df$impt <- case_when(mod.df$beta > 0 & mod.df$lcl > 0 ~ 1,
#                          mod.df$beta < 0 & mod.df$ucl < 0 ~ 1,
#                          !is.na(mod.df$beta) ~ 0)
# mod.df$impt <- ifelse(is.na(mod.df$lcl), NA, mod.df$impt)

colnames(mod.df) <- c("est", "SE", "z", "p", "pop", "cov")
mod.df$impt <- ifelse(mod.df$p < 0.05, 1, 0)

View(mod.df)

mod.df %>% 
  group_by(cov) %>% 
  summarize(n=n(),
            impt = sum(impt, na.rm=TRUE),
            perc = impt/n) %>% arrange(-perc)

write.csv(mod.df, "Data/SingleMLModels_P.csv", row.names = FALSE)
mod.df <- read.csv("Data/SingleMLModels_P.csv")

#4. Full models----

#9a. All data----
pops <- unique(dat$pop)
coefs.all <- data.frame()
for(i in 1:length(pops)){
  
  dat.i <- dat %>% 
    dplyr::filter(pop==pops[i])
  
  if(pops[i]!=16){
    fit <- glmer.nb(count ~ years + 
                      pdsis_breed  + 
                      pdsis_fall  + 
                      pdsis_winter  + 
                      pdsis_spring  + 
                      crops_breed  + 
                      crops_fall  + 
                      crops_winter  + 
                      crops_spring  + 
                      trees_breed  + 
                      trees_fall  + 
                      trees_winter  + 
                      trees_spring  + 
                      (1|route), offset = p,  data=dat.i)
  }
  else{
    fit <- glmer.nb(count ~ years + 
                      pdsis_breed  + 
                      pdsis_winter  + 
                      pdsis_spring  + 
                      crops_breed  + 
                      crops_winter  + 
                      crops_spring  + 
                      trees_breed  + 
                      trees_winter  + 
                      trees_spring  + 
                      (1|route), offset = p,  data=dat.i)
  }
  
  coefs.all<- data.frame(summary(fit)[['coefficients']]) %>% 
    mutate(pop = pops[i],
           cov = row.names(summary(fit)[['coefficients']]),
           method = "All") %>% 
    rbind(coefs.all)
  
  print(paste0("Finished population ", i, " of ", length(pops)))
}

#9b. BBS only----
coefs.bbs <- data.frame()
for(i in 1:length(pops)){
  
  dat.i <- dat %>% 
    dplyr::filter(pop==pops[i]) %>% 
    dplyr::filter(method=="BBS")
  
  if(pops[i]!=16){
    fit <- glmer.nb(count ~ years + 
                   pdsis_breed  + 
                   pdsis_fall  + 
                   pdsis_winter  + 
                   pdsis_spring  + 
                   crops_breed  + 
                   crops_fall  + 
                   crops_winter  + 
                   crops_spring  + 
                   trees_breed  + 
                   trees_fall  + 
                   trees_winter  + 
                   trees_spring  + 
                   (1|route), offset = p,  data=dat.i)
  }
  else{
    fit <- glmer.nb(count ~ years + 
                   pdsis_breed  + 
                   pdsis_winter  + 
                   pdsis_spring  + 
                   crops_breed  + 
                   crops_winter  + 
                   crops_spring  + 
                   trees_breed  + 
                   trees_winter  + 
                   trees_spring  + 
                   (1|route), offset = p,  data=dat.i)
  }
  
  coefs.bbs<- data.frame(summary(fit)[['coefficients']]) %>% 
    mutate(pop = pops[i],
           cov = row.names(summary(fit)[['coefficients']]),
           method = "BBS") %>% 
    rbind(coefs.bbs)
  
  print(paste0("Finished population ", i, " of ", length(pops)))
}
coefs <- rbind(coefs.all, coefs.bbs)

colnames(coefs) <- c("est", "SE", "z", "p", "pop", "cov", "method")
write.csv(coefs, "MLFullModelResults.csv", row.names=FALSE)
coefs <- read.csv("MLFullModelResults.csv")

#10. Look at ML results----
coef.sum <- coefs %>% 
  mutate(sig = ifelse(p < 0.05, 1, 0)) %>% 
  group_by(cov, method) %>% 
  summarize(est = mean(est),
            p = mean(p),
            sig = sum(sig)/n()) %>% 
  ungroup() %>% 
  dplyr::filter(!cov %in% c("(Intercept)", "years")) %>% 
  separate(cov, into=c("cov", "season")) %>% 
  mutate(sig = ifelse(est < 0, -sig, sig))

write.csv(coef.sum, "MLFullModelResults_Summary.csv", row.names = FALSE)

ggplot(coef.sum) +
  geom_raster(aes(x=season, y=cov, fill=sig)) +
  facet_wrap(~method) +
  scale_fill_gradient2(low="red", mid="white", high = "blue", midpoint=0)

coef.grid <- coef.sum %>% 
  dplyr::select(method, cov, season, sig) %>% 
  pivot_wider(id_cols = c(method, cov), names_from=season, values_from=sig) %>% 
  arrange(method, cov)

regions <- read_sf("Data/PopulationPolygons.shp")

regions.breed <- regions %>% 
  dplyr::filter(stage=="breed") %>% 
  mutate(pop = as.integer(pop)) %>% 
  left_join(coefs %>% 
              dplyr::filter(cov=="years") %>% 
              mutate(sig = ifelse(p < 0.05, 1, 0)))

write_sf(regions.breed, "Data/PopulationPolygons_Results_Breed.shp")

regions.coef <- regions %>% 
  mutate(pop = as.integer(pop)) %>% 
  left_join(coefs %>% 
              dplyr::filter(!cov %in% c("(Intercept)", "years")) %>% 
              separate(cov, into=c("cov", "stage")) %>% 
              mutate(sig = ifelse(p < 0.05, 1, 0)))

write_sf(regions.coef, "Data/PopulationPolygons_Results.shp")

regions.breed.sig <- regions.breed %>% 
  dplyr::filter(sig==1)

regions.breed.nonsig <- regions.breed %>% 
  dplyr::filter(sig==0)

ggplot() +
  geom_sf(data=regions.breed.sig, aes(fill=est)) +
  geom_sf(data=regions.breed.nonsig, aes(fill=NA), fill="grey50") +
  facet_wrap(~method) +
  scale_fill_gradient2(low="red", mid="white", high = "blue", midpoint=0)
