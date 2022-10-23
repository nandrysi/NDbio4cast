022 ND BIOS Biological Forecasting Class Repository Zhuoran Yu email: zyu3@nd.edu

Kayla Anderson email: kander42@nd.edu

 main
Stacy Mowry email: smowry@nd.edu

Nick Andrysiak email: nandrysi@nd.edu

Stacy Mowry
email: smowry@nd.edu

Nick Andrysiak
email: nandrysi@nd.edu

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

#Load Target Data
download_targets <- function(){
  readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)
}
target<-download_targets()


chla<-target$observation [target$site_id == "PRPO" & target$variable== "chla"]


#Downlodading chemical prop of surface water data
download_chem_prop<-function()
{
  chem_prop <- loadByProduct(dpID="DP1.20093.001", 
                                    site=c("PRPO"),
                                    check.size = F)
                                    
          #subset desired data
time<-chem_prop$swc_fieldSuperParent$collectDate
time<-lubridate::as_datetime(time)
oxygen<-chem_prop$swc_fieldSuperParent$dissolvedOxygen
temp<-chem_prop$swc_fieldSuperParent$waterTemp
chem_prop<-cbind(time,oxygen,temp)
  
  
  #return the dataframe
  chem_prop<<- as.data.frame(chem_prop)
  
  
  
}

oxygen<-chem_prop$oxygen
  plot (chla, chem_prop$oxygen)
  
  
  

