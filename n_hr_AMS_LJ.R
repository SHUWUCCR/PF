rm(list=ls())
#============================================================
# This code computes the AMS for different storm durations (1-hr, 2-hr, 3-hr, 6-hr, 12-hr, 18-hr, 24-hr, 5-d, 10-d) in 90% completed years in the period of 1948-2014 for 74 selected stations.
# The criteria of station selection in IL and collar county stations around IL with at least 20 years of 90% data + those for combinations.
#===========================================================
#--------------------Work Directory------------------------
# ==============edit your own directory==========
setwd("C:/Users/lujin2/N_hr_AMS_code_with_sample/")
MainPath <- 'C:/Users/lujin2/'

# =====Install and download packages used in this study====
# install packages before using
library("dplyr") # for using "%>%" operater
library("readr") # for read in csv file
library("stringr") # for using str_sub function
library("base") # for using rowSums function

#------------------------Read in station list and completed year list------------------------
# read in the station list 
sta <- read.csv('./Hourly_Stn_List_74.csv')
sta$ID <- as.character(sta$ID)
# read in the list of completed years in each station
year.list <- read.csv('./PRCP_Hourly_Station_Continuity_Table_74_stn.csv')
# set the list of 9 different time periods and number of hours
periods <- c('1hr', '2hr', '3hr', '6hr', '12hr', '18hr', '24hr', '5d', '10d')
num_hour <- c(1,2,3,6,12,18,24,5*24,10*24)
index <- 1 # Change from 1 to 9 to compute different durations of AMS
p <- periods[index] # time periods
n <- num_hour[index] # number of hours, match above

#----------------------loop through each station and compute the AMS-------------------------
# check i=1 in the loop, sample station 'COOP112417 DOWNS 2 NE, IL US'
for (i in 1:length(sta$ID)) {
  data.file <- read.csv(paste0('./', gsub(':','',sta$ID[i]), ' ',sta$NAME[i],'.csv')) # Read in station precipitation data
  # extract the completed year list for this station
  yr.temp <- as.character(year.list[year.list$ID==sta$ID[i],7])
  yr.list <- as.numeric(strsplit(yr.temp, " ")[[1]])
  # record the first recorded year 
  fst.yr <- data.file$YEAR[1]
  
  # Change all missing values in data set into 0
  data <- data.file
  data$NAME <- sta$NAME[i]
  data[(data$VALUE==99999)|(is.na(data$VALUE)),]$VALUE <- 0

 # --------------Start getting annual max------------
  this.annmax <- data.frame() # create an empty data frame for annual max series of this station
  # loop through each year in the completed year list
  for (y in 1:length(yr.list)) {
    yr <- yr.list[y]
    # check if add last n hours data from last year, list the data period we need
    if(yr!=fst.yr & n!=1){
      fst.day <- which(data$YEAR==yr & data$MON==1 & data$DAY==1 & data$TIME==0) #first hour of first day in a year
      yr.data <- rbind(data[(fst.day-n+1):(fst.day-1),],filter(data,YEAR==yr))
    }else{
      yr.data <- filter(data,YEAR==yr)
    }
    yr.val <- yr.data$VALUE
    m<-length(yr.val)
    # sum up n-hr PRCP
    add.data <- rowSums(sapply(1:n, function(x) yr.val[x:(x+m-n)]))
    row.index <- which.max(add.data) # find the row of the maximum precipitation in that year, the date of this row should be the end date of the n-hr
    thisyr <- data.frame(yr.data[(row.index+n-1), c(6,10,1:4,8:9)], max(add.data))
    colnames(thisyr) <- c('ID',	'NAME', 'YEAR', 'MON', 'DAY', 'TIME', 'MFLAG',	'QFLAG', 'VALUE') # notice that 'YEAR' and 'AMS_YR' are different, AMS for one year can happen in the perious year in some situations
    this.annmax <- rbind(this.annmax,thisyr) # bind AMS for this year to the total AMS df
  }
    
  # Write each station's annual max series to csv
  stationAM.file <- paste0('./AMS_',p,'_', gsub(':','',sta$ID[i]), ' ',sta$NAME[i],'.csv')
  write.csv(this.annmax, file=stationAM.file, row.names = FALSE)
}
