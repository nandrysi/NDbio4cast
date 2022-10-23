022 ND BIOS Biological Forecasting Class Repository Zhuoran Yu email: zyu3@nd.edu

Kayla Anderson email: kander42@nd.edu

 main
Stacy Mowry email: smowry@nd.edu

Nick Andrysiak email: nandrysi@nd.edu

## Download data (work in progress)
Currently downloads chla data, dissolved O2, phytoplankton biomass, and water temp, checks start and end date of data, resolution of data, and plots data time series


##load libraries
#install.packages("neonUtilities")
library(neonUtilities)
library(lubridate)


#Load Target Data
download_targets <- function(){
  readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)
}

target<-download_targets()

##subset target data
target<-subset(target, target$site_id == "PRPO")
target<-subset(target, target$variable == "chla")

##plot target data over time
plot(target$datetime,target$observation,type="p")

##Start and end date of data availability (2018-2022)
summary(target$datetime)

##resolution of data = DAILY
table(target$datetime)

##% NAs
na.perc<-sum(is.na(target$observation))/nrow(target)

## pull predictor variables 
## dissolved o2 is available monthly (2016-2022)
download_dissolved.O2<-function()
{
  dissolved.gasses <- loadByProduct(dpID="DP1.20097.001", 
                                    site=c("PRPO"),
                                    check.size = F)
  
  ## extract date and dissolved 02 data
  time.gas<-dissolved.gasses$sdg_fieldSuperParent$collectDate
  time.gas<-lubridate::as_datetime(time.gas)
  dissolved.02<-dissolved.gasses$sdg_fieldSuperParent$dissolvedOxygen
  dissolved.02.time<-cbind(time.gas,dissolved.02)
  
  ##plot time series
  plot(time.gas,dissolved.02, type="p")
  

  print(summary(time.gas))
  print(table(time.gas))

  
}


##Biomass 2014-2021, every few months
download_biomass<-function()
{
  biomass <- loadByProduct(dpID="DP1.20163.001", 
                           site=c("PRPO"),
                           check.size = F)
  
  ## extract date and biomass data
  time.biomass<-biomass$alg_fieldData$collectDate
  time.biomass<-lubridate::as_datetime(time.biomass)
  phyto<-biomass$alg_fieldData$phytoDepth1
  phyto.time<-cbind(time.biomass,phyto)
  
  ##plot time series
  plot(time.biomass, phyto,type="p")
  
  print(summary(time.biomass))
  print(table(time.biomass))
  
  
}

##water temperature available from 2018-22
##data resolution every 30 minutes

download_watertemp<-function()
{
  watertemp <- loadByProduct(dpID="DP1.20264.001", 
                             site=c("PRPO"),
                             check.size = F)
  
  ## extract date and water temp data 
  time.watertemp<-watertemp$TSD_30_min$startDateTime
  time.watertemp<-lubridate::as_datetime(time.watertemp)
  temp_depth<-watertemp$TSD_30_min$tsdWaterTempMean
  watertemp.time<-cbind(time.watertemp,temp_depth)
  
  ##plot time series
  plot(time.watertemp, temp_depth,type="p")
  
  print(summary(time.watertemp))
  print(table(time.watertemp))
 
}
