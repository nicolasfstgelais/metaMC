library(missForest)

data=read.csv("data/nico_function_input.csv",row.names=1,stringsAsFactors = F)

#impute missing values using a random forest
data=missForest(data);data=data$ximp

#import the sf for error 
dist=read.csv("data/dist.csv",row.names=1)


# function to run the metabolism function with error terms
# error terms (sd) are specified in the dist file, when sd=0 = no error
data.rand<-function(data,dist,n=1000){
  output=list()
  for(i in 1:nrow(data)){
    output[[rownames(data)[i]]]=metabolism(mapply(rnorm,mean=data[i,],n=1000,sd=dist[,1]))}
  
return(output)
}
 

#function to summarize the output from metabolism
# give the mean and sd for each output
sumStats<-function(l){
  
  for(i in names(l)){
    mean.res=apply((l[[i]]),2,mean)
    names(mean.res)=paste0(names(mean.res),".mean")
    
    sd.res=apply((l[[i]]),2,sd)
    names(sd.res)=paste0(names(sd.res),".sd")
    
    sumStats=c(mean.res,sd.res)
    sumStats=sumStats[order(names(sumStats))]
    sumStats=as.matrix(sumStats);colnames(sumStats)=i
    if(i==names(l)[1])output=sumStats
    if(i!=names(l)[1])output=cbind(output,sumStats)
    
  }
  
  return(t(output))
}


metabolism <- function(data){
  
  data=as.data.frame(data)
  
  #measured terms in mass balance from field and lab:
  data$doSat<-data$o2.g.m3 * 100/data$do.pct.sat #calculate dissolved oxygen concentration at equilibrium with atmosphere 
  data$ro18o2<-((data$delo18.o2/1000)*(0.00205))+(0.00205)#converting del value of DO-O18 back to O18:O16 ratio
  data$ro18h2o <- ((data$delo18.h2o/1000)*(0.00205))+(0.00205) #converting del value of H2O-O18 back to O18:O16 ratio
  data$o18o2 <- data$ro18o2 / (1 + data$ro18o2) #converting the ratio of O18:O16 in DO to an atomic fraction following Bogard et al. (2017)
  data$o18h2o <- data$ro18h2o / (1 + data$ro18h2o) #converting the ratio of O18:O16 in H2O to an atomic fraction following Bogard et al. (2017)
  
  #calculating the gas exchange coefficient for O2 empirically from lake area and wind speed:
  data$u10 <- data$wind.ms * (1+ (((0.0013^(0.5))/0.41) * (log(10/data$wind.height.ms)))) #calculating wind speed at 10 meter height following equation 3 Vachon & Prairie (2013)
  data$k600cmh <- 2.51 + 1.48*data$u10 + 0.39*data$u10*(log10(data$lake.area.km2)) #k600 in cm/h from table 2 equation B vachon & prairie 2013
  data$k600md <- data$k600cmh*24/100 #converting k600 to m/d
  data$sco2 <- 1800.6-(120.1*data$h2o.temp.c)+(3.7818*(data$h2o.temp.c^2))-(0.047608*(data$h2o.temp.c^3))#calculating schmidt number for oxygen from Jahne et al (1987)
  data$ko2md <- data$k600md*((data$sco2/600)^(-2/3)) #converting k600 to ko2 in m/d for use in mass balance
  #2/3 power used for wind speed less than 3.7 m/s following Vachon et al. (2010) and Guerin et al. (2007)
  
  
  #fixed terms in mass balance derived from literature:
  ro18air <- 0.00209895 #O18:O16 ratio of atmospheric O2 from 23.88 permil value of Barkan & Luz (2005) 
  o18air <- ro18air / (1 + ro18air) #converting the ratio of O18:O16 in aatmospheric O2 to an atomic fraction following Bogard et al. (2017)
  ffg <- 0.9972 # air-water exchange fractionation effect from Knox et al. (1992)*error*
  ffs <- 1.0007 #solubility fractionation effect from Benson & Kraus (1984)*error*
  ffp <- 1 # photosynthesis fractionation effect from Guy et al. (1993)*error*
  ffr <- 0.985 #respiration fractionation effect from Bogard et al. (2017)*error*
  

  #mass balance 
  #volumetric GPP as g DO/m3/d from equation 5 of Bogard et al (2017)
  gppv <- (data$ko2md/data$zmix.m) * ((data$o2.g.m3* (data$o18o2*ffg - data$o18o2*ffr)) - (data$doSat * (o18air*ffs*ffg - data$o18o2*ffr))) / (data$o18h2o*ffp - data$o18o2*ffr) 
  
  #areal GPP as g DO/m2/d for the mixed layer
  gppa <- gppv * data$zmix.m
  
  #R-volumetric as g DO/m3/d from equation 5 of Bogard et al (2017)
  resv <- (data$ko2md/data$zmix.m) * ((data$o2.g.m3 * (data$o18o2*ffg - data$o18h2o*ffp)) - (data$do.pct.sat* (o18air*ffs*ffg - data$o18h2o*ffp))) / (data$o18h2o*ffp - data$o18o2*ffr) 
  
  #R-areal as g DO/m2/d for the mixed layer
  resa <- resv * data$zmix.m
  
  #net ecosystem production - volumetric
  nepv <- gppv - resv
  
  #net ecosystem production - areal
  nepa <- nepv * data$zmix.m
  
  #The ratio of GPP to R 
  gpptor <- gppv/resv
  
  output<- data.frame(gppv=gppv,gppa=gppa,resv=resv,resa=resa,nepv=nepv,nepa=nepa,gpptor=gpptor)
  
  return(output)
}


#run the metabolism function with an associated error
res=data.rand(data,dist,1000)

# summary stats
sumStats(res)

write.csv(sumStats(res),"output/output.csv")



