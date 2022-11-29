#4. Define model----

sink("pop_model.jags")
cat("
    model {
    
    ######################## LIKELIHOOD #########################
    
    for (k in 1:nroutes){
      for (t in 1:nyears){
        noise[t,k] ~ dnorm(0, tau.noise)
        count[t,k] ~ dpois(lambda[t,k])
        log(lambda[t,k]) <- trend * year[t]   # Trend
            + beta[1]   * pdsis_breed[t,]      
            + beta[2]   * pdsis_fall[t,]        
            + beta[3]   * pdsis_winter[t,]        
            + beta[4]   * pdsis_spring[t,]        
            + beta[5]   * crops_breed[t,]        
            + beta[6]   * crops_fall[t,]      
            + beta[7]   * crops_winter[t,]      
            + beta[8]   * crops_spring[t,]      
            + beta[9]   * trees_breed[t,]      
            + beta[10]  * trees_fall[t,]      
            + beta[11]  * trees_winter[t,]   
            + beta[12]  * trees_spring[t,]   
            + alpha                           # Population intercept             
            + route[k]                        # Route effect
            + noise[t,k]                      # Noise
            + log(p[t,k])                     # Offset for availability of detection
            
    ### GOF statistics ###
        pres[t,k] <- (count[t,k] - lambda[t,k]) / sqrt(lambda[t,k])
        count.new[t,k] ~ dpois(lambda[t,k])
        pres.new[t,k] <- (count.new[t,k] - lambda[t,k]) / sqrt(lambda[t,k])
      }  # End t loop
    } # End k loop
    
    ####################### PRIORS ###############################
    
    ### Priors on population intercepts, trends, & covariate effects
    alpha   ~ dnorm(0, 0.001)
    trend   ~ dnorm(0, 0.001)
    for (j in 1:nvars){
      theta[j] ~ dnorm(0, 0.001)
    } # End j loop
    
    ### Priors on novice observer effect, noise sd, & pop intercept sd
    tau.noise  ~ dgamma(0.001, 0.001)
    sd.noise  <-  1/pow(tau.noise,0.5)
    tau.alpha <- 1/(sd.alpha*sd.alpha)
    sd.alpha  ~ dunif(0, 10)
    
    ### Beta coefficients
    for(j in 1:nvars) {
      g[j] ~ dbern(0.5*pop.ind[j])
      beta[j] <- theta[j]*g[j]
    }
    
    ### Priors on route effects
    for (r in 1:nroutes){
      route[r] ~ dnorm(0, tau.rte)T(-5,5)
    }
    tau.rte  ~ dgamma(0.001,0.001)
    sd.rte  <-  1/pow(tau.rte,0.5)
    
    ################ DERIVED PARAMETERS ######################
    
    fit <- sum(pres[,])
    fit.new <- sum(pres.new[,])
#    resid.tot <- sum(resid[,])
    
    for(t in 1:nyears){
      meanN[t] <- mean(lambda[t,])
    }
    
#    B.tot <- (meanN[nyears]/meanN[1])^(1/nyears)
#    trend.tot <- 100*(B.tot-1)
    } # End model
    ",fill = TRUE)
sink()

#5. Model settings----

# Initial values
inits <- function()list(alpha=runif(1,-1,1), 
                        trend=runif(1,-1,1), 
                        route = runif(ncol(count.pop), -1,1), 
                        tau.noise=10,
                        sd.alpha=1, 
                        tau.obs=1, tau.rte=10) 

# Parameters monitored
parameters.gvs <- c("sd.alpha", "beta", "g", "sd.noise", "fit", 
                    "fit.new", "sd.rte", "trend.tot", "meanN", "resid.tot")

#MCMC settings
ni <- 40000
nt <- 5
nb <- 10000
nc <- 2

#6. Set up loop----
pops <- unique(dat$pop)
for(i in 1:length(pops)){
  cat("Population", pops[i])
  
  #7. Format data----
  dat.i <- dat %>% 
    dplyr::filter(pop==pops[i]) %>% 
    group_by(pop, year, route) %>% 
    sample_n(1) %>% 
    ungroup() 
  #Randomly sample 1 per route - figure out why this is necessary later
  
  nroutes.i <- length(unique(dat.i$route))
  nyears.i <- length(unique(dat.i$year))
  
  count.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=count) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  
  p.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=p) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  
  pdsis_breed.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=pdsis_breed) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  pdsis_fall.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=pdsis_fall) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  pdsis_winter.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=pdsis_winter) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  pdsis_spring.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=pdsis_spring) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  crops_breed.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=crops_breed) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  crops_fall.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=crops_fall) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  crops_winter.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=crops_winter) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  crops_spring.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=crops_spring) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  trees_breed.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=trees_breed) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  trees_fall.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=trees_fall) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  trees_winter.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=trees_winter) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  trees_spring.i <- dat.i %>% 
    pivot_wider(id_cols=c(year), names_from = route, values_from=trees_spring) %>% 
    dplyr::select(-year) %>% 
    as.matrix()
  
  pop.data <- list(count=count.i, p=p.i,
                   nvars = 12, nroutes=nroutes.i,
                   nyears=nyears.i, year = dat.i$years,
                   pdsis_breed = pdsis_breed.i,
                   pdsis_fall = pdsis_fall.i,
                   pdsis_winter = pdsis_winter.i,
                   pdsis_spring = pdsis_spring.i,
                   crops_breed = crops_breed.i,
                   crops_fall = crops_fall.i,
                   crops_winter = crops_winter.i,
                   crops_spring = crops_spring.i,
                   trees_breed = trees_breed.i,
                   trees_fall = trees_fall.i,
                   trees_winter = trees_winter.i,
                   trees_spring = trees_spring.i,
                   pop.ind=rep(1, 12))
  
  
  #8. Call Jags----
  mod.name <- paste("fit_pop", pops[i], sep="")
  mod.file <- paste("fit2_pop", pops[i], ".csv", sep="")
  mod.name <- jags(pop.data, inits=NULL, parameters.gvs, "pop_model.jags", 
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
                   working.directory = getwd())
  
  pop_mcmc <- as.mcmc(mod.name$BUGSoutput$sims.matrix)
  write.csv(pop_mcmc,file=mod.file)
  
  print(mod.name, digits = 3)
}

