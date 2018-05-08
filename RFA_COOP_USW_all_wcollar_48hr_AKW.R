rm(list=ls())
#============================================================
# This code use the l-moments to calculate the 2, 5, 10, 25, 50, 100, 500 return periods for the 48-hr storms, for COOP and USW stations with at least 30 years of 90% complete data in IL and collar counties for 1957-2016. 
# Currently all data are using 'GEV' distribution for fitting of frequency curves

# Use 1.04 for 2-day to 48-hr conversion

#============================================================
#--------------------Work Directory--------------------------
# The working directory of this script is "\\swsatlas\Albedo\B70\DataAnalysis for Major GHCND\R"
setwd("\\\\swsatlas/Albedo/B70/DataAnalysis for Major GHCND/R/")
MainPath <- '\\\\swsatlas/Albedo/B70/DataAnalysis for Major GHCND/'
# setwd("C:/Users/kwang47/Documents/DataAnalysis for Major GHCND/R/")
# MainPath <- 'C:/Users/kwang47/Documents/DataAnalysis for Major GHCND/'
DataPath <- paste0(MainPath,'AMS_48hr/Stn_AMS_flagfixed_wMDPR_combined/')

library(lmom)
library(lmomRFA)
library(readr)
library(dplyr)
library(stringr)

distr <- c('glo','gev','gno','pe3','gpa')
nonexprob <- 1-1/c(2.54,5.52,10.51,25,50,100,500) # calculate non-exceedance probability, with conversion from AMS to PDS (Partial duration series)
returnperiod <- c('2-yr','5-yr','10-yr','25-yr','50-yr','100-yr','500-yr')
sta.d <- read_csv(paste0(MainPath,'Station_continuity/COOP_USW_forAMS_wcollar_AKW_d.csv')) # This station list name refer to the station list name in the AMS R code
regionlist <- unique(sta.d$REGION)
nyrs <- 60 # this code use 60 years of data
rn <- 1957:2016 # this code use data in 1957-2016, this is used later to name the rows in AMS dataframe

regstats <- data.frame() # dataframe for regional statistics for all IL and collar county stations, like l-moments and goodness-of-fit, etc.
reggrowthcurve <- data.frame() # dataframe for regional quantiles for all IL and collar county stations, not the real 'regional growth curve (non-dimentional)'
reg.L <- data.frame() # dataframe for lower CL of regional quantiles for all IL and collar county stations
reg.U <- data.frame() # dataframe for upper CL of regional quantiles for all IL and collar county stations
siteq.map <- data.frame() # dataframe for station quantiles for all IL and collar county stations

for (Region in regionlist) {
  
  lmompath <- paste0(MainPath,'RFA_48hr/',Region,'/')
  
  # Get station list
  stationlist  <- sta.d[sta.d$REGION==Region,]
  station <- str_replace(stationlist$NAME, "([,])", "")
  
  # ---------Get ams for all stations in the region-------------
  ams.total <- data.frame()
  for (i in 1:nrow(stationlist)) {
    data_name <- paste0(str_sub(stationlist$ID[i],7,17), ' ', str_replace(stationlist$NAME[i], "([,])", ""))
    data.file <- paste0(DataPath, data_name, '.csv')
    ams <- read_csv(data.file)
    # ams$NAME <- station[i]
    ams.total <- rbind(ams.total,ams)
  }
  
  # Set up ams matrix for RFA
  data <- data.frame(matrix(rep(NA,nyrs*length(station)),ncol=length(station)))
  rownames(data) <- rn
  colnames(data) <- station
  
  # Re-Write the AMS values to the input format of RFA analysis
  for (i in rn) {
    for (j in 1:length(station)) {
      if (i %in% ams.total$YEAR[ams.total$NAME==station[j]]) {
        data[which(rn==i),j] <- ams.total$VALUE[ams.total$NAME==station[j] & ams.total$YEAR==i]
      }
    }
  }
  data <- data[,colSums(!is.na(data))>=5] # remove those stations with less than 5 years of ams
  
  # ------------Start RFA analysis--------
  sitelmom<-regsamlmu(data) # Comput L-moments for each site
  sta <- regtst(sitelmom,nsim=500) # Compute discordancy, heterogeneity and goodness-of-fit measures for regional frequency analysis
  dis=sta$D
  het=sta$H
  Z=sta$Z
  bestfit <- distr[which.min(abs(Z))]
  
  #rfit=regfit(sitelmom,bestfit) # fit a frequency distribution of the best goodness-of-fit to a vector of regional average L-moments.
  rfit=regfit(sitelmom,'gev') # Manually select the distribution after checking the best-fit distributions in 'Regional l stats ...' csv file if necessary. Here use 'GEV' for all regions
  sq <- sitequant(nonexprob, rfit) # calculate station quantiles, before convertion
  regq <- regquant(nonexprob, rfit)*mean(sitelmom$l_1)*1.04 # get regional averaged quantiles, 1.04 is the convertion from daily measurement to 48hr precipitation
  reggrowthcurve <- rbind(reggrowthcurve,cbind(Region,matrix(regq,nrow=1))) # This is te actual quantiles, not growth curves
  
  correlmat=cor(data,use = "pairwise.complete.obs",method = "pearson")
  avecor=(sum(correlmat, na.rm=T)-length(station))/(2*0.5*length(station)*(length(station)-1))
  
  sim=regsimq(rfit$qfunc,cor=avecor, nrec=sitelmom$n,nrep=500, f=nonexprob,boundprob = c(0.025,0.975)) # Compute error bounds for a fitted regional frequency distribution
  sqb=sitequantbounds(sim,rfit) # Compute quantiles' CLs for individual stations, before convertion
  monte=sim$sim.rgcratio*1.04 #a matrix of dimension length(f)*nrep, containing the simulated values of the ratio of the estimated to the true regional growth curve for quantiles corresponding to probabilities in f
  matresult.est <-matrix(rep(NA,length(nonexprob)*length(station)),nrow = length(station))
  rownames(matresult.est) <- station
  colnames(matresult.est) <- returnperiod
  matresult.LCL <- matresult.est
  matresult.UCL <- matresult.est
  
  # Write station quantiles and CLs of this region to dataframes, with conversion from daily measurement to 48hr precipitation
  for (j in 1:length(nonexprob)) {
    for (i in 1:length(colnames(data))) {
      matresult.est[i,j] <- sq[i,j]*1.04 # 
      matresult.LCL[i,j] <- matrix(unlist(sqb[[i]][j,4]),nrow = 1,byrow = T)*1.04
      matresult.UCL[i,j] <- matrix(unlist(sqb[[i]][j,5]),nrow = 1,byrow = T)*1.04
    }
  }
  
  # combine regional values to dataframes of whole IL
  sitestats<-cbind(sitelmom,dis)
  regstats <- rbind(regstats,cbind(Region,het[1],t(Z),bestfit,avecor))
  reg.L <- rbind(reg.L,cbind(Region,matrix(colMeans(matresult.LCL),nrow=1)))
  reg.U <- rbind(reg.U,cbind(Region,matrix(colMeans(matresult.UCL),nrow=1)))
  siteq.map <- rbind(siteq.map,cbind(stationlist[,c(1,4,5)],matresult.est))
  
  write.csv(matresult.est,file = paste0(lmompath,'Frequency estimates_',Region,'_48hr_IL.csv'))
  write.csv(matresult.LCL,file = paste0(lmompath,'Frequency estimates Lower CL_',Region,'_48hr_IL.csv'))
  write.csv(matresult.UCL,file = paste0(lmompath,'Frequency estimates Upper CL_',Region,'_48hr_IL.csv'))
  write.csv(sitestats,file = paste0(lmompath,'Site l stats_',Region,'_48hr_IL.csv'),row.names = F)
  
  # the monte carlo results are not saved here, will be added later when applying weighted function
  
  
}

#--------------------------Get Reginal l stats-----------
colnames(regstats) <- c('Region','Heterogeneity',distr,'Bestfit','AverageCorrelation')
colnames(reggrowthcurve) <- c('Region',returnperiod)
colnames(reg.L) <- c('Region',returnperiod)
colnames(reg.U) <- c('Region',returnperiod)

write.csv(regstats,file = paste0(MainPath,'RFA_48hr/Regional l stats_IL.csv'),row.names = F)
write.csv(reggrowthcurve,file = paste0(MainPath,'RFA_48hr/Regional frequency estimate_48hr_IL.csv'),row.names = F)
write.csv(reg.L,file = paste0(MainPath,'RFA_48hr/Regional frequency estimate lower CL_48hr_IL.csv'),row.names = F)
write.csv(reg.U,file = paste0(MainPath,'RFA_48hr/Regional frequency estimate upper CL_48hr_IL.csv'),row.names = F)
write.csv(siteq.map,file = paste0(MainPath,'RFA_48hr/Site frequency estimates for map_48hr_IL.csv'),row.names = T)
