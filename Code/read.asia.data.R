# read in the PIKE data, MIKE centroid data, MIKE population estimates, etc
# make list by MIKEsiteID of which site is in which data source.

UNRegion.select <- "Asia"


#----------------------------------------------
# Get the MIKE centroids and other details.
mike.centers <- readxl::read_excel(file.path("..","Data","190212_mike_sites_list_UpdatedTo2018.xlsx"),
                                   sheet="mike_sites_list_2018")
mike.centers <- plyr::rename(mike.centers, c("siteid"="MIKEsiteID",
                                             "name"="MIKEsiteName",
                                             "un_region"="UNRegion",
                                             "subregion"="SubregionName"))
data.source <- data.frame(MIKEsiteID=unique(mike.centers$MIKEsiteID[ mike.centers$UNRegion==UNRegion.select]), GIS=TRUE, stringsAsFactors=FALSE)






#----------------------------------------------
# read in the raw PIKE data
pike <- read.csv(file.path("..","Data","2020-05-14_PIKEStatsUpTo2019FusionTableFormat.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)

# only from 2003 onwards. Notice that no pdata available from SE Asia for 2018/2019 and it will be imputed by the GLMM
# based on past years.
cat("\n\n Exclude year 2002 from PIKE. Here is the number of records by year \n")
pike <- pike[ pike$year >= 2003,]
xtabs(~year, data=pike)

# select Africa data only
xtabs(~UNRegion, data=pike, exclude=NULL, na.action=na.pass)

cat("\n\nRestricting PIKE to those countries in ", UNRegion.select, "\n")
pike<- pike[ pike$UNRegion == UNRegion.select,]

pike.original <- pike
# exclude site-years with 0 carcasses reported as not useful for the analysis
select <- pike$TotalNumberOfCarcasses == 0
sum(select)
N.pike.site.years.0.carcasses <- sum(select)
pike <- pike[ !select,]

data.source <- merge(data.source, data.frame(MIKEsiteID=unique(pike$MIKEsiteID), pike=TRUE, stringsAsFactors=FALSE), all=TRUE)
data.source[ is.na(data.source)] <- FALSE
apply(data.source[,-1], 2, sum)

# find out total number of carcasses reported on
temp <- plyr::ddply(pike, "MIKEsiteID", plyr::summarize,
                    TC =sum(TotalNumberOfCarcasses))
temp[ temp$TC==0,] # any sites with 0 carcasses.
temp


#----------------------------------------------
# Get the MIKE population estimates 
# Currently there are only a few estimates for sites (about 1/3 are missing) and so we expand the
# present estimates to all years

mike.pop.est.wide <- readxl::read_excel(file.path("..","Data","MIKE.site.pop.ests-asia.xlsx"),
                                        sheet="Sheet1")
mike.pop.est.wide$MIKE.site.name <- NULL
temp <- plyr::adply(mike.pop.est.wide, 1, function(x){
    data.frame(x, year=min(pike$year):max(pike$year), stringsAsFactors=FALSE)
})
mike.pop.est.wide <- reshape2::dcast(temp, MIKEsiteID~year, value.var="population")

setdiff(mike.pop.est.wide$MIKEsiteID, mike.centers$MIKEsiteID)
setdiff(mike.centers$MIKEsiteID[mike.centers$UNRegion=="Africa"], mike.pop.est.wide$MIKEsiteID)


# Convert from wide to long and remove sites with all missing data
mike.pop.est <- reshape2::melt(mike.pop.est.wide,
                               id.var=c("MIKEsiteID"),
                               variable.name="year",
                               value.name="population")
mike.pop.est$year <- as.numeric(as.character(mike.pop.est$year))
# remove any missing values for population values
mike.pop.est <- mike.pop.est[ !is.na(mike.pop.est$population),]

data.source <- merge(data.source, data.frame(MIKEsiteID=unique(mike.pop.est$MIKEsiteID), Population=TRUE, stringsAsFactors=FALSE), all=TRUE)
data.source$Population[ is.na(data.source$Population)] <- FALSE

# Asia data doesn't have full population data yet, so I will use the mean for that subregion for missing population data for a placeholder
mike.pop.est <- reshape2::melt(mike.pop.est.wide,
                               id.var=c("MIKEsiteID"),
                               variable.name="year",
                               value.name="population")
mike.pop.est$year <- as.numeric(as.character(mike.pop.est$year))

# add the subregion information to the MIKE population estimates
mike.pop.est <- merge(mike.pop.est, mike.centers[,c("MIKEsiteID","SubregionName" )], all.x=TRUE)
cat("\n\nCheck that subregion names is available for all MIKE population sites\n")
mike.pop.est[ is.na(mike.pop.est$SubregionName),]

# use the subregional average abundance when abundance is missing.  THIS IS TEMPORARY
select <- is.na(mike.pop.est$population)
cat("Missing population information \n")
mike.pop.est[select,]
mike.pop.est <- plyr::ddply(mike.pop.est, "SubregionName", function(x){
    mean.pop <- mean(x$population, na.rm=TRUE)
    select <- is.na(x$population)
    x$population[ select] <- round(mean.pop)
    x$impute <- ""
    x$impute[select] <- "imputed mean"
    x
})
cat("Missing population estimates replaced by the mean population \n")
mike.pop.est[mike.pop.est$impute != "",]



mike.pop.est.original <- mike.pop.est
# add in the mike site name
mike.pop.est.original <- merge(mike.pop.est.original, mike.centers[,c("MIKEsiteID","MIKEsiteName")], all.x=TRUE)
mike.pop.est.original[ is.na(mike.pop.est.original$MIKEsiteName),]

# find out when the survey took place. This is indicated by bold face in the MIKE population workbook
# see the help for tidyxl on detail on how to extract this information

mike.pop.survey.years <- NULL

rm(mike.pop.cells, mike.pop.formats, bold_indices)



#----------------------------------------------
# Get the lsmeans for Asia. These are only available up to 2017

lsmeans.cont <- read.csv(file.path("..","Data","SC73_as_lsmns_90CL.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
lsmeans.cont <- plyr::rename(lsmeans.cont, c("yr"="year","estimate"="mean", "region"="Region"))
# compute 95% ci
lsmeans.cont$X2.5.  <- lsmeans.cont$mean - 1.96*lsmeans.cont$se
lsmeans.cont$X97.5. <- lsmeans.cont$mean + 1.96*lsmeans.cont$se
lsmeans.cont$Source <- "LSMeans"
lsmeans.cont$df       <- NULL
lsmeans.cont$t.stat   <- NULL
lsmeans.cont$p.value  <- NULL
lsmeans.cont$lower.ci <- NULL
lsmeans.cont$upper.ci <- NULL

pike$pike <- pike$NumberOfIllegalCarcasses/pike$TotalNumberOfCarcasses

pike$yearC <- as.factor(pike$year)
asia.fit <- lm( pike ~ SubregionID + yearC, data=pike, weight=TotalNumberOfCarcasses)
asia.emmo <- emmeans::emmeans(asia.fit, ~yearC)
asia.fit.lsmeans.cont <- summary(asia.emmo, infer=TRUE)

asia.fit.lsmeans.cont$year <- as.numeric(as.character(asia.fit.lsmeans.cont$yearC))
asia.fit.lsmeans.cont <- plyr::rename(asia.fit.lsmeans.cont, c("emmean"="mean",
                                             "SE"="se",
                                             "lower.CL"="X2.5.",
                                             "upper.CL"="X97.5."))
asia.fit.lsmeans.cont$df <- NULL
asia.fit.lsmeans.cont$t.ratio <- NULL
asia.fit.lsmeans.cont$p.value <- NULL
asia.fit.lsmeans.cont$yearC   <- NULL
asia.fit.lsmeans.cont$Source  <- "LSMeans"
asia.fit.lsmeans.cont$region  <- UNRegion.select

# subregional lsmeans computed similarly to Africa with a country effect
lsmeans.reg <- plyr::ddply(pike, "SubregionName", function(pike){
   fit <- lm( pike ~ CountryCode + yearC, data=pike, weight=TotalNumberOfCarcasses)
   fit.emmo <- emmeans::emmeans(fit, ~yearC)
   lsmeans  <- summary(fit.emmo, infer=TRUE)
   lsmeans
})
lsmeans.reg$year <- as.numeric(as.character(lsmeans.reg$yearC))
lsmeans.reg <- plyr::rename(lsmeans.reg, c("emmean"="mean",
                                             "SE"="se",
                                             "lower.CL"="X2.5.",
                                             "upper.CL"="X97.5."))
lsmeans.reg$df      <- NULL
lsmeans.reg$t.ratio <- NULL
lsmeans.reg$p.value <- NULL
lsmeans.reg$yearC   <- NULL
lsmeans.reg$Source  <- "LSMeans"




rm(asia.fit, asia.emmo, temp)


#----------------------------------------------
## Now include only those sites for which we have GIS and pike data.
# For Asia we don't worry about the population data because it is missing for so many site
# Problematic sites
select <- apply(!data.source[,c("GIS","pike")],1,any)
temp <- data.source[select,]
temp

select <- data.source$pike #& data.source$Population
MIKEsiteID.select <- as.character(data.source$MIKEsiteID[select])


# only those sites matching the population data
cat("\n\nRestricting PIKE to match sites with population data. These sites have no MIKE pop data\n")
select <- !pike$MIKEsiteID %in% MIKEsiteID.select
pike[select,]
pike <- pike[ !select,]

# restricting MIKE pop estimates to match PIKE data
cat("\n\nRestricting MIKE population estimates to match PIKE data. These data are excluded\n")
select <- !mike.pop.est$MIKEsiteID %in% MIKEsiteID.select
mike.pop.est[select,]

mike.pop.est <- mike.pop.est[ !select,]


# get the base map of Asia

# Geographic co-rdinate system is GCS_North_American_1983, projected to NAD_1983_BC_Environment_Albers 
# http://spatialreference.org/ref/sr-org/82/ has the projection string
proj4string = "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 

google.map <- ggmap::get_map(c(90,10), maptype="toner-lite",  source="stamen", zoom=4)
base.map <- ggmap(google.map)

google.map2 <- ggmap::get_map(c(left=60, bottom=-30, right=120, top =40), 
                              maptype="toner-lite",  source="stamen", zoom=4)
base.map2 <- ggmap(google.map2)



