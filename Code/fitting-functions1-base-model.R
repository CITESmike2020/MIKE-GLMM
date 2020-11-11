# Model to fit the PIKE data
# Non-Hierarchical model for year effects to make the marginal means of year effect (on logit scale)
# We need to have a very wide prior on the effect to make the marginal means match


fit.pike <- function(pike, mike.pop.est,
                     IC.var   ="NumberOfIllegalCarcasses", # default variable in data frame
                     TC.var   ="TotalNumberOfCarcasses",
                     site.var ="MIKEsiteID",
                     year.var ="year",
                     subregion.var="SubregionName",
                     seed=NA
                     ){
# Fit a hierarchical Bayesian model to pike data
# Parametera are
#    pike - data frame with # of carcasses monitored (TC) and # illegally killed carcasses (IC)
#    mike.pop.est in long format 
  
# various checks on the data frame
if( !IC.var        %in% names(pike)){stop("Missing variable for number of intentionally killed carcasses")}
if( !TC.var        %in% names(pike)){stop("Missing variable for total number of carcasses")}
if( !site.var      %in% names(pike)){stop("Missing variable for site id in PIKE data")}
if( !year.var      %in% names(pike)){stop("Missing variable for year in PIKE data")}
if( !subregion.var %in% names(pike)){stop("Missing variable for subregion in PIKE data")}
    
if( is.na(seed)){stop("Must specific an initial seed")}  
# check that all of pike data has comparable mike.pop.est
if( length(setdiff(pike$MIKEsiteID, mike.pop.est$MIKEsiteID)) > 0){stop("PIKE sites not present in MIKE sites")}
if( length(setdiff(mike.pop.est$MIKEsiteID, pike$MIKEsiteID)) > 0){stop("MIKE sites not present in PIKE data")}
if( !site.var %in% names(mike.pop.est)){stop("Missing variable for site id in population data")}
if( !year.var %in% names(mike.pop.est)){stop("Missing variable for year in population data")}
 
# make sure there is no current site, year, or subregion index variable on the PIKE or population dataset
pike$site.index <- NULL
pike$year.index <- NULL
pike$subregion.index <- NULL

mike.pop.est$site.index <- NULL
mike.pop.est$year.index <- NULL
mike.pop.est$subregion.index <- NULL

# JAGS does not deal with character data. So we need to create an index number for each site
site.index <- unique(pike[,site.var, drop=FALSE])
site.index <- site.index[ order(site.index[,site.var]),,drop=FALSE]
site.index$site.index <- 1:nrow(site.index)

pike <- merge(pike, site.index, all.x=TRUE)

# Similarly, we need to create a subregion-index
subregion.index <- unique(pike[, subregion.var, drop=FALSE])
subregion.index <- subregion.index[ order(subregion.index[, subregion.var]),,drop=FALSE]
subregion.index$subregion.index <- 1:nrow(subregion.index)
pike <- merge(pike, subregion.index, all.x=TRUE)

# Set up the year index for convenience
year.index <- unique(pike[,year.var, drop=FALSE])
year.index <- year.index[ order(year.index[,year.var]),,drop=FALSE]
year.index$year.index <- 1:nrow(year.index)

pike <- merge(pike, year.index, all.x=TRUE)

#browser()
# get the subregion for each site in the right order
site.subregion = unique(pike[,c("site.index", "subregion.index")])
site.subregion <- site.subregion[ order(site.subregion$site.index),]
site.subregion.index <- site.subregion$subregion.index


# add the index number to the mike population data and create the population matrix
mike.pop.est <- merge(mike.pop.est, year.index)
mike.pop.est <- merge(mike.pop.est, site.index)
mike.pop <- matrix(NA, nrow=length(unique(mike.pop.est$year.index)), ncol=length(unique(mike.pop.est$site.index)))
mike.pop[ as.matrix(mike.pop.est[,c("year.index","site.index")])] <- mike.pop.est$population


# Save the BUGS code to a character variable
# The model file.
# The cat() command is used to save the model to the working directory.
# Notice that you CANNOT have any " (double quotes) in the bugs code
# between the start and end of the cat("...",) command.

model <- function() {
############################################################
    
  # model with 
  #  non-hierarchical model for year
  #  random effects (hierarchical model) for site and site-year  
  for(i in 1:Nsiteyear.obs){
     IC[i] ~ dbin(pi[i], TC[i])
     logit(pi[i]) <- Year.eff[ year.index.obs[i] ]+
                     Site.eff[ site.index.obs[i] ]+
                     Year.Site.eff[year.index.obs[i], site.index.obs[i] ]
  }
  
  # set up predicted Year.Site estimates of PIKE on logit scale for THESE particular sites and years
  for(y in 1:Nyear){
    for(s in 1:Nsite){
       Year.Site.est[y,s] <- Year.eff[ y ]+
                     Site.eff[ s ]+
                     Year.Site.eff[ y, s ]
    }
  }

  # Convert the Year.Site estimates from logit to [0,1] scale
  for(y in 1:Nyear){
     for(s in 1:Nsite){
        Year.SiteP.est[y,s] <- ilogit(Year.Site.est[y,s])
     }
  }
  
  
  # hierarchical model for site effects
  for(i in 1:Nsite){
    Site.eff[i] ~ dnorm(0, tau.site.eff)
  }
  sd.site.eff <- 1/sqrt(tau.site.eff)
  tau.site.eff       ~ dgamma(.001, .001)
  prior.tau.site.eff ~ dgamma(.001, .001)  # get a sample from the prior
  
  # hierarchical model for year.site effects
  for(s in 1:Nsite){
    for(y in 1:Nyear){
      Year.Site.eff[y,s] ~ dnorm(0, tau.year.site.eff)
    }
  }
  sd.year.site.eff <- 1/sqrt(tau.year.site.eff)
  tau.year.site.eff       ~ dgamma(.001, .001)
  prior.tau.year.site.eff ~ dgamma(.001, .001)  # get a sample from prior
  # non-hierarchical  effects for year effects
  # Year.eff2 is a "dummy" to make it consistent with hierarchical effects
  for(y in 1:Nyear){
     Year.eff      [y] ~ dnorm(0, .001)  
     prior.Year.eff[y] ~ dnorm(0, .001)  # get a sample from prior
     Year.eff2[y] <- Year.eff[y]
  }
  
  # also compute the year effect on the [0,1] scale
  for(y in 1:Nyear){
     YearP.eff[y] <- ilogit(Year.eff[y])    
  }
  
  
  # Unweighted marginal means (on logit scale) and on the p-scale
  # The unweighted marginal mean should match the Year.eff
  for(y in 1:Nyear){
     Year.est.MM [y] <- mean( Year.Site.est [y, 1:Nsite])
     YearP.est.MM[y] <- mean( Year.SiteP.est[y, 1:Nsite])
  }

  # Unweighted marginal means on the p-scale at the subregional level
  for(y in 1:Nyear){
     for(subregion in 1:Nsubregion){
        YearP.sub.est.MM[y,subregion]= sum( Year.SiteP.est[y, 1:Nsite]*(site.subregion.index[]==subregion))/
                                       sum(site.subregion.index[]==subregion)
     }
  }
   
  # Weighted marginal means on p-scale using the mike pop weights
  for(y in 1:Nyear){
     YearP.est.MMw[y] <- inprod( Year.SiteP.est[y, 1:Nsite], mike.pop[y, 1:Nsite])/ sum(mike.pop[y, 1:Nsite])
  }  
  

  # Compute the omnibus goodness of fit statistic for observed and simulated data
  # to perform the bayesian p-value computation
  for(i in 1:Nsiteyear.obs){
    # expected intentionally killed
    EIC[i] <- TC[i] * pi[i]
    # simulate number intentionally killed
    sim.IC[i] ~ dbin(pi[i], TC[i])
    # tukey freeman stats
    TF.stat.obs.indiv[i] <- (sqrt(IC[i])     -sqrt(EIC[i]))^2
    TF.stat.sim.indiv[i] <- (sqrt(sim.IC[i])- sqrt(EIC[i]))^2
  }
  # Get the Tukey-Freeman statistics
  TF.stat.obs <- sum( TF.stat.obs.indiv[])
  TF.stat.sim <- sum( TF.stat.sim.indiv[])
  
  # Count the number of simulated 0's
  Zero.sim <- sum( sim.IC[] == 0)
  
  # uses a bootstrap to sample from the generated site.year with replacement to
  # simulate a bayesian bootstrap sample.
  
  #Site.select[1:Nsite] ~ dmulti(rep(1,Nsite),Nsite)
  Site.select[1:Nsite] ~ ddirch(rep(1,Nsite))
  
  # compute a boot strap sample from the Year.site logit values (observed and imputed)
  # compute a boot strap sample from the Year.site anti.logit values (observed and imputed)
  # compute a boot strap sample using the weighted means 
  for(y in 1:Nyear){
     Year.est.MM.boot  [y] <- inprod( Site.select[1:Nsite], Year.Site.est [y, 1:Nsite])
     YearP.est.MM.boot [y] <- inprod( Site.select[1:Nsite], Year.SiteP.est[y, 1:Nsite])
  }
  # Weighted marginal means on p-scale using the mike pop weights
  for(y in 1:Nyear){
     YearP.est.MMw.boot[y] <- inprod( Site.select[1:Nsite], 
                                      Year.SiteP.est[y, 1:Nsite]* mike.pop[y, 1:Nsite]/
                                        inprod(Site.select[1:Nsite], mike.pop[y, 1:Nsite]))
  } 
} # End of the model

# set up functions
logit <- function(x){log(x/(1-x))}
expit <- function(x){1/(1+exp(-x))}

Nyear      <- nrow(year.index)
Nsite      <- nrow(site.index)
Nsubregion <- nrow(subregion.index)



# set up the data list
# The datalist will be passed to JAGS with the names of the data
# values.
data.list <- list(
  # observed site years
  Nsiteyear.obs      = nrow(pike),
  year.index.obs     = pike$year.index, 
  site.index.obs     = pike$site.index,
  IC        = pike[, IC.var], 
  TC        = pike[, TC.var],
  
  # site-subregion index in site index order
  site.subregion.index = site.subregion.index,
  
  # number of year and number of sites
  Nyear     = Nyear,
  Nsite     = Nsite,
  Nsubregion= Nsubregion,

  # MIKE population data
  mike.pop  = mike.pop   # matrix [year,site] of values
)

# initial value list
init.list <- list(
      list(Year.eff=logit(runif(Nyear, min=.1, max=.9))),
      list(Year.eff=logit(runif(Nyear, min=.1, max=.9))),
      list(Year.eff=logit(runif(Nyear, min=.1, max=.9)))
      )  # end of list of lists of initial values

## Next create the list of parameters to monitor.
# The deviance is automatically monitored.
# 
monitor.list <- c( "Year.eff","Year.eff2","YearP.eff",
                  "Year.est.MM","YearP.est.MM","YearP.est.MMw", # marginal means
                  "Site.eff","sd.site.eff","tau.site.eff",
                  "Year.Site.eff","sd.year.site.eff","tau.year.site.eff",
                  "Year.Site.est",  "Year.SiteP.est",# Year-Site estimates and marginal means
                  "YearP.sub.est.MM",  # marginal mean for subregions
                  "TF.stat.obs","TF.stat.sim",   # omnibus goodness of fit statistic
                  "Zero.sim",                    # number of simulated 0 counts
                  "prior.tau.site.eff", "prior.tau.year.site.eff","prior.Year.eff", # get samples of priors
                  "Year.est.MM.boot","YearP.est.MM.boot" ,"YearP.est.MMw.boot" # marginal means at using bootstrapping to mimic the Year.eff terms
                  ) # parameters to monitor
 


# Finally, the actual call to JAGS
set.seed(seed)  # intitalize seed for MCMC 

results <- jags( 
      data      =data.list,   # list of data variables
      inits     =init.list,   # list/function for initial values
      parameters=monitor.list,# list of parameters to monitor
      model.file=model,       # function with bugs model
      n.chains=3,
      n.iter  =5000,          # total iterations INCLUDING burn in
      n.burnin=2000,          # number of burning iterations
      n.thin=2,               # how much to thin
      DIC=TRUE,               # is DIC to be computed?
      working.dir=getwd()     # store temporary results in current working directory
      )

# return the fitted values
list(results=results,
     pike=pike,
     site.index=site.index, 
     year.index=year.index,
     subregion.index=subregion.index,
     IC.var=IC.var, 
     TC.var=TC.var,
     site.var=site.var,
     year.var=year.var,
     model=model,
     data.list,
     init.list,
     monitor.list,
     seed=seed
     )

}

