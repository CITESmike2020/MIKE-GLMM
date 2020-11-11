# Investigate the impact of dropping each site on the PIKE values
# This is run separately from the main document and the results are stored in a separate
# file for display by the main document


library(ggplot2)
library(ggforce)
library(plyr)
library(reshape2)
library(R2jags)
library(tidyverse)

# load the fitting functions
source("fitting-functions1-base-model.R")  # non-hierarchical model for year effect
source("fitting-functions2.R")  # extraction and other functions


# Get the data and restrict to Africa only
source("read.africa.data.R", echo=TRUE)


#####################################################################################
# Now we want to drop each site in turn and refit the Bayesian model.
# We should be able to parallelize this computation over the 50+ sites

doParallel <- FALSE # should I set up parallel processing of 
if(doParallel) {
   library(doMC) # for parallel model fitting
   # see http://viktoriawagner.weebly.com/blog/five-steps-to-
   detectCores()
   cl <- makeCluster(4)
   # Need to export some libraries to the cluster
   # see http://stackoverflow.com/questions/18981932/logging-
   #clusterEvalQ(cl, library(unmarked))
   registerDoMC(5)
}

site.names <- unique(pike$MIKEsiteID)

# now fit the model, dropping each site in turn
set.seed(4233423)
site.seeds <- runif(length(site.names),min=1, max=1000000000)

site.names <- data.frame(MIKEsiteID=site.names, seed=site.seeds)

drop.fits <- plyr::dlply(site.names,  "MIKEsiteID", function(x, pike, mike.pop.est){
   cat("\n\n*** Starting ", as.character(x$MIKEsiteID), "\n")
   # drop data corresponding to this site
   pike         <- pike        [ pike$MIKEsiteID         != x$MIKEsiteID,]
   mike.pop.est <- mike.pop.est[ mike.pop.est$MIKEsiteID != x$MIKEsiteID,]
   
   #print(setdiff(mike.pop.est$MIKEsiteID, pike$MIKEsiteID))
   #print(setdiff(pike$MIKEsiteID, mike.pop.est$MIKEsiteID))
   #browser()
   # do the fit
   fit <- fit.pike(pike, mike.pop.est,  seed=x$seed )

   # return the fit
   fit
}, pike=pike, mike.pop.est=mike.pop.est, .parallel=doParallel)


names(drop.fits)

# extract the values from the fits after dropping terms
# The [0,1] based on the estimated logit values
drop.eff <- plyr::ldply(drop.fits, function(x){
  effect  <- extract.effect(x,  effect.name="^YearP.eff\\[",    index.type="year", source="Drop")
  as.data.frame(effect)
})
head(drop.eff)

write.csv(drop.eff, file="africa.sens.drop.YearP.eff.csv", row.names=FALSE)


# The unweighted marginal mean of the $PIKE$ computed on the [0,1] scale is:
drop.eff <- plyr::ldply(drop.fits, function(x){
  effect  <- extract.effect(x,  effect.name="^YearP.est.MM\\[",    index.type="year", source="Drop")
  as.data.frame(effect)
})
head(drop.eff)
write.csv(drop.eff, file="africa.sens.drop.YearP.eff.MM.csv", row.names=FALSE)


# the weighted marginal mean of the $PIKE$ computed on the [0,1] scale is:
drop.eff <- plyr::ldply(drop.fits, function(x){
  effect  <- extract.effect(x,  effect.name="^YearP.est.MMw\\[",    index.type="year", source="Drop")
  as.data.frame(effect)
})
head(drop.eff)
write.csv(drop.eff, file="africa.sens.drop.YearP.eff.MMw.csv", row.names=FALSE)


if(doParallel) stopCluster(cl)
