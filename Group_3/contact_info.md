022 ND BIOS Biological Forecasting Class Repository Zhuoran Yu email: zyu3@nd.edu

Kayla Anderson email: kander42@nd.edu

 main
Stacy Mowry email: smowry@nd.edu

Nick Andrysiak email: nandrysi@nd.edu

## Download data (work in progress)
Currently downloads chla data, dissolved O2, and phytoplankton biomass and plots data time series


##load libraries
#install.packages("neonUtilities")
library(neonUtilities)
library(lubridate)

##pull forecast target data
download_target<-function()
{
  water.quality <- loadByProduct(dpID="DP1.20288.001", 
                                 site=c("PRPO"),
                                 check.size = F)
  
  ## extract date and chlorophyll data from water quality data 
  time.water<-water.quality$waq_instantaneous$startDateTime
  time.water<-lubridate::as_datetime(time.water)
  chlorophyll<-water.quality$waq_instantaneous$chlorophyll
  chlorophyll.time<-cbind(time.water,chlorophyll)
  
  ##plot time series
  plot(time.water,chlorophyll, type="l")
  
  return(chlorophyll.time)
}

## pull predictor variables

download_dissolved.O2<-function()
{
  dissolved.gasses <- loadByProduct(dpID="DP1.20097.001", 
                                    site=c("PRPO"),
                                    check.size = F)
  
  ## extract date and chlorophyll data from water quality data 
  time.gas<-dissolved.gasses$sdg_fieldSuperParent$collectDate
  time.gas<-lubridate::as_datetime(time.gas)
  dissolved.02<-dissolved.gasses$sdg_fieldSuperParent$dissolvedOxygen
  dissolved.02.time<-cbind(time.gas,dissolved.02)
  
  ##plot time series
  plot(time.gas,dissolved.02, type="l")
  
  return(dissolved.02.time)
  
}



download_biomass<-function()
{
  biomass <- loadByProduct(dpID="DP1.20163.001", 
                           site=c("PRPO"),
                           check.size = F)
  
  ## extract date and chlorophyll data from water quality data 
  time.biomass<-biomass$alg_fieldData$collectDate
  time.biomass<-lubridate::as_datetime(time.biomass)
  phyto<-biomass$alg_fieldData$phytoDepth1
  phyto.time<-cbind(time.biomass,phyto)
  
  ##plot time series
  plot(time.biomass, phyto,type="l")
  
  return(phyto.time)
  
}

Stacy Mowry
email: smowry@nd.edu

Nick Andrysiak
email: nandrysi@nd.edu

