# read in the PIKE data, MIKE centroid data, MIKE population estimates, etc
# make list by MIKEsiteID of which site is in which data source.

#----------------------------------------------
# Get the MIKE centroids and other details.
mike.centers <- readxl::read_excel(file.path("..","Data","190212_mike_sites_list_UpdatedTo2018.xlsx"),
                                   sheet="mike_sites_list_2018")
mike.centers <- plyr::rename(mike.centers, c("siteid"="MIKEsiteID",
                                             "name"="MIKEsiteName",
                                             "un_region"="UNRegion",
                                             "subregion"="SubregionName"))
data.source <- data.frame(MIKEsiteID=unique(mike.centers$MIKEsiteID[ mike.centers$UNRegion=="Africa"]), GIS=TRUE, stringsAsFactors=FALSE)


#----------------------------------------------
# Get the lsmeans for Africa 
lsmeans.cont <- read.csv(file.path("..","Data","SC73_af_lsmns_90CL.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
lsmeans.cont <- plyr::rename(lsmeans.cont, c("yr"="year","estimate"="mean"))
# compute 95% ci
lsmeans.cont$X2.5.  <- lsmeans.cont$mean - 1.96*lsmeans.cont$se
lsmeans.cont$X97.5. <- lsmeans.cont$mean + 1.96*lsmeans.cont$se
lsmeans.cont$Source <- "LSMeans"
lsmeans.cont$df       <- NULL
lsmeans.cont$t.stat   <- NULL
lsmeans.cont$p.value  <- NULL
lsmeans.cont$lower.ci <- NULL
lsmeans.cont$upper.ci <- NULL

# Get the lsmeans for Africa on the regional level
lsmeans.reg <- read.csv(file.path("..","Data","SC73_sr_lsmns_90CL.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
lsmeans.reg <- plyr::rename(lsmeans.reg, c("yr"="year","estimate"="mean","region"="SubregionName"))
# compute 95% ci
lsmeans.reg$X2.5.  <- lsmeans.reg$mean - 1.96*lsmeans.reg$se
lsmeans.reg$X97.5. <- lsmeans.reg$mean + 1.96*lsmeans.reg$se
lsmeans.reg$Source <- "LSMeans"
lsmeans.reg$df       <- NULL
lsmeans.reg$t.stat   <- NULL
lsmeans.reg$p.value  <- NULL
lsmeans.reg$lower.ci <- NULL
lsmeans.reg$upper.ci <- NULL




#----------------------------------------------
# Get the MIKE population estimates over time
# the 2019 values are currently just duplicated from the 2018 values.
mike.pop.est.wide <- readxl::read_excel(file.path("..","Data","MIKE site pop ests.xlsx"),
                                        sheet="Sheet1")
setdiff(mike.pop.est.wide$MIKEsiteID, mike.centers$MIKEsiteID)
setdiff(mike.centers$MIKEsiteID[mike.centers$UNRegion=="Africa"], mike.pop.est.wide$MIKEsiteID)


# Convert from wide to long and remove sites with all missing data
mike.pop.est.wide$CountryName <- NULL
mike.pop.est.wide$Long.MIKEsiteName <- NULL
mike.pop.est <- reshape2::melt(mike.pop.est.wide,
                               id.var=c("MIKEsiteID"),
                               variable.name="year",
                               value.name="population")
mike.pop.est$year <- as.numeric(as.character(mike.pop.est$year))
# remove any missing values for population values
mike.pop.est <- mike.pop.est[ !is.na(mike.pop.est$population),]

data.source <- merge(data.source, data.frame(MIKEsiteID=unique(mike.pop.est$MIKEsiteID), Population=TRUE, stringsAsFactors=FALSE), all=TRUE)

# add the subregion information to the MIKE population estimates
mike.pop.est <- merge(mike.pop.est, mike.centers[,c("MIKEsiteID","SubregionName" )], all.x=TRUE)
cat("\n\nCheck that subregion names is available for all MIKE population sites\n")
mike.pop.est[ is.na(mike.pop.est$SubregionName),]
mike.pop.est$SubregionName[ mike.pop.est$MIKEsiteID == "BBR"] <- "West Africa"
mike.pop.est$SubregionName[ mike.pop.est$MIKEsiteID == "ALE"] <- "West Africa"
mike.pop.est[ is.na(mike.pop.est$SubregionName),]


mike.pop.est.original <- mike.pop.est
# add in the mike site name
mike.pop.est.original <- merge(mike.pop.est.original, mike.centers[,c("MIKEsiteID","MIKEsiteName")], all.x=TRUE)
mike.pop.est.original[ is.na(mike.pop.est.original$MIKEsiteName),]
mike.pop.est.original$MIKEsiteName[ mike.pop.est$MIKEsiteID == "BBR"] <- "Babban Rafi Forest"
mike.pop.est.original$MIKEsiteName[ mike.pop.est$MIKEsiteID == "ALE"] <- "Monte Alen"
mike.pop.est.original[ is.na(mike.pop.est.original$MIKEsiteName),]

# gind out when the survey took place. This is indicated by bold face in the MIKE population workbook
# see the help for tidyxl on detail on how to extract this information

mike.pop.cells   <- tidyxl::xlsx_cells  (file.path("..","Data","MIKE site pop ests.xlsx"))
mike.pop.formats <- tidyxl::xlsx_formats(file.path("..","Data","MIKE site pop ests.xlsx"))

# find all cells in bold which indicates when surveys were conducted
bold_indices <- which(mike.pop.formats$local$font$bold)
mike.pop.survey.years <- mike.pop.cells[mike.pop.cells$local_format_id %in% bold_indices, ]
mike.pop.survey.years <- mike.pop.survey.years[ mike.pop.survey.years$sheet=="Sheet1",]  # only sheet1
mike.pop.survey.years$MIKEsiteID <- mike.pop.est.wide$MIKEsiteID[ mike.pop.survey.years$row-1]
mike.pop.survey.years$year       <- as.numeric(names(mike.pop.est.wide)[ mike.pop.survey.years$col-2])
mike.pop.survey.years <- mike.pop.survey.years[, c("MIKEsiteID","year")]
mike.pop.survey.years <- merge(mike.pop.survey.years, mike.pop.est.original[,c("MIKEsiteID","year","MIKEsiteName","SubregionName")], all.x=TRUE)

rm(mike.pop.cells, mike.pop.formats, bold_indices)
#----------------------------------------------
# read in the raw PIKE data
#pike <- read.csv(file.path("..","Data","190307_PikeStatsUpTo2018FusionTableFormat.csv"),     header=TRUE, as.is=TRUE, strip.white=TRUE)
 pike <- read.csv(file.path("..","Data","2020-05-14_PIKEStatsUpTo2019FusionTableFormat.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)

# only from 2003 onwards
cat("\n\n Exclude year 2002 from PIKE. Here is the number of records by year \n")
pike <- pike[ pike$year >= 2003,]
xtabs(~year, data=pike)

# select Africa data only
UNRegion.select <- "Africa"
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
## Now include only those sites for which we have population and pike data.
# Problematic sites
select <- apply(!data.source[,-1],1,any)
temp <- data.source[select,]
temp

select <- data.source$pike & data.source$Population
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


# get the base map of africa

# Geographic co-rdinate system is GCS_North_American_1983, projected to NAD_1983_BC_Environment_Albers 
# http://spatialreference.org/ref/sr-org/82/ has the projection string
proj4string = "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 

google.map <- ggmap::get_map(c(20,0), maptype="toner-lite",  source="stamen", zoom=3)
base.map <- ggmap(google.map)

google.map2 <- ggmap::get_map(c(20,0), maptype="toner-lite",  source="stamen", zoom=4)
google.map2 <- ggmap::get_map(c(left=-20, bottom=-30, right=50, top =20), 
                              maptype="toner-lite",  source="stamen", zoom=4)
base.map2 <- ggmap(google.map2)



