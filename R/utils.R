#METHOD OF MOMENTS CONVERSIONS FOR WEIBULL DISTRIBUTION
.w_shape <- function(mean, sd){
  return((sd/mean)^-1.086)
}
.w_scale <- function(mean, sd){
  return(mean/gamma(1 + 1/.w_shape(mean, sd)))
}

#RETURNS THE PDF OF THE GIVEN DISTRIBUTION
.getDistPDF <- function(pars, dist){
  if (dist=="normal")
    return(function(x) dnorm(x, mean=pars["mean"], sd=pars["sd"]))
  else if (dist=="lognormal")
    return(function(x) dlnorm(x, meanlog=pars["meanlog"], sdlog=pars["sdlog"]))
  # 02FEB2025 NEW WEIBULL DISTRIBUTION #
  else if (dist=="Weibull")
    return(function(x) dweibull(x, .w_shape(pars["mean"],pars["sd"]),.w_scale(pars["mean"],pars["sd"])))
}

#RETURNS THE RANDOM NUMBER GENERATOR OF THE GIVEN DISTRIBUTION
.getDistRNG <- function(pars, dist){
  if (dist=="normal")
    return(function(x) rnorm(x, mean=pars["mean"], sd=pars["sd"]))
  else if (dist=="lognormal")
    return(function(x) rlnorm(x, meanlog=pars["meanlog"], sdlog=pars["sdlog"]))
  # 02FEB2025 NEW WEIBULL DISTRIBUTION #
  else if (dist=="Weibull")
    return(function(x) rweibull(x, .w_shape(pars["mean"],pars["sd"]),.w_scale(pars["mean"],pars["sd"])))
}

#RETURNS THE DESCRIPTION OF A GIVEN DISTRIBUTION
.getDistDesc <- function(pars, dist) {
  if (dist=="normal")
    return(paste("normal(mean=", round(pars["mean"], 2), ", ", "sd=", round(pars["sd"], 2),  ")", sep=""))
  else if (dist=="lognormal")
    return(paste("lognormal(meanlog=", round(pars["meanlog"], 2), ", ", "sdlog=", round(pars["sdlog"], 2), ")", sep=""))
  # 02FEB2025 NEW WEIBULL DISTRIBUTION #
  else if (dist=="Weibull")
    return(paste("Weibull(mean=", round(pars["mean"], 2), ", ", "sd=", round(pars["sd"], 2), ")", sep=""))
}

#RETRIEVES THE MEAN OF THE PRIOR FOR MOSQUITO BITING FREQUENCY
.getPriorMosqBitingFreqMean <- function(pars, dist) {
  if (dist=="normal")
    return(pars["mean"])
  else if (dist=="lognormal")
    return(exp(pars["meanlog"]+pars["sdlog"]^2/2))
}

#SMOOTHS A TIME SERIES WITH +- N STEPS AVERAGE PER POINT
.smoothUDSeries<- function(series, n){
  return(as.numeric(stats::filter(series, rep(1/n, n), sides=2, circular=TRUE)))
}

#CALCULATES MOSQ biting RATE COMPONENT DEPENDENT ON HUMIDITY IN TIME ##
.mvse_hum_effect_aV<- function(U,meanU){
  ef<- (U-meanU)/sqrt(1+(U-meanU)^2) ##Maintained the sign since the effect of VPD should increase biting per Brown et al., 2022 [DOI: 10.1111/ele.14228]
  return(ef)
  ##no need to trim the function for bio meaning
}

#CALCULATES MOSQ DEATH RATE COMPONENT DEPENDENT ON TEMPERATURE IN TIME
.mvse_temp_effect_muV_aedes <- function(temp){
  ef <- (0.8692-0.1590*temp+(0.01116*(temp^2))-0.0003408*(temp^3)+0.000003809*(temp^4))
  return(ef)
}
#NEW CULEX PIPIENS VERSION from 10.1007/s00436-014-4152-x
#EDITED TO RE-ENABLE TRIMMING FOR MAX LIFESPAN & PREVENT NEGATIVE
.mvse_temp_effect_muV <- function(temp){
  ef <- (-5.24*temp + 178.32)
  ef[which(ef > 120)] <- 120 # Max known lifespan of Culex pipiens from 10.1007/s00436-014-4152-x
  ef[which(ef <= 0)] <- .16 # fix negative numbers as bio -> trim to zero
  ef[which(temp > 34)] <- .16
  return(1/ef)
}

#CALCULATES MOSQ DEATH RATE COMPONENT DEPENDENT ON HUMIDITY IN TIME
.mvse_hum_effect_muV <- function(U, meanU){
  ef <- (U-meanU)/sqrt(1+(U-meanU)^2) ##Flipped the sign since the effect of VPD is opposite the effect of RH
  return(ef)
  ##no need to trim the function for bio meaning
}

#CALCULATES EPSVH (PROB TRANS VECTOR TO HUMAN) DEPENDENT ON TEMPERATURE IN TIME
.mvse_temp_epsVH <- function(temp){ ##testing
  ep <- 0.001044*temp*(temp-12.286)*((32.461-temp)^(0.5))
  ep[which(ep<0)] <- 0 # fix negative numbers as bio -> trim to zero
  ep[which(is.nan(ep))] <- 0
  # next two added to debug extreme values
    ep[which(ep>1)] <- 1
    ep[which(temp<0)] <- 0 #fix for bio meaning below zero there is no transmission
  return(ep)
}
##### Lambrechts et al. 10.1073/pnas.1101377108
#CULEX PIPIENS version from DOI: 10.1111/tbed.14513
# .mvse_temp_epsVH <- function(temp){ ## NEW ##
#   ep <- 0.00103*temp*(temp-14.02)*((34.76-temp)^(0.5)) #*.81 prob host is bird
#   ep[which(ep<0)] <- 0 # fix negative numbers as bio -> trim to zero
#   ep[which(is.nan(ep))] <- 0
#   ep[which(ep>1)] <- 1
#   ep[which(temp<0)] <- 0 #fix for bio meaning below zero there is no transmission
#   return(ep)
# }
# JAN 2025 Ultimately I decided not to use the new formulation for epsVH because the old one was already based on Culex mosquitos
# and the new formulation was for VC not for just epsVH. The old one used int he OG package is documented in Lambrechts et al., 2011

#CALCULATES EPSHV (PROB TRANS HUMAN TO VECTOR) DEPENDENT ON TEMPERATURE IN TIME
.mvse_temp_epsHV <- function(temp){ ##testing
  ep <- 0.0729*temp-0.9037
  ep[temp<12.4] <- 0
  ep[temp>26.1] <- 1
  return(ep)
}

#CALCULATES MOSQ INCUBATION TO INFECTION RATE DEPENDENT ON TEMPERATURE IN TIME
.mvse_temp_effect_gammaV <- function(temp){
  Tk <- temp+273.15
  R <- 1.987 # cal/(mol*K)
  ef <- (24.0*(0.003359*(Tk/298.)*exp((15000./R)*(1/298.-1./Tk))/(1.+ exp((6.203*(10^21)/R)*(1./(-2.176*(10^30))-1./Tk)))))
  ef[which(ef<0)]<- 0 #fix negative numbers as bio -> trim to zero
  return(ef)
}

##### 
## NEW FUNCTIONS ADDED BY HMJ ##

# Changelog:
# Replace .mvse_temp_effect_muV (Aedes death rate given temp) with a Culex pipiens version.
# Replace .mvse_temp_epsVH (Aedes V > H prob transmit) with a Culex pipiens version.
# Create new functions to ultimately calculate Vapor Pressure Deficit using Temp and RH
# Create .mvse_temp_effect_aV function to approximate the effect of temperature on biting rate
#   new function will need to be implemented deeper in mvsemodel where biting rate is calculated

## NEW AS OF 11 AUG 2024 - FUNCTIONS TO COMPUTE VAPOR PRESSURE DEFICIT as in DOI:10.1111/ele.14228
.tetens <- function(temp){return(.61078*exp(17.27*temp / (temp + 237.3)))}
.dewpoint <- function(temp, rh){return(temp - (100-rh)/5)} ##Using Lawrence 2005. http://dx.doi.org/10.1175/BAMS-86-2-225
# See: https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf
.vpd <- function(temp, rh){return(.tetens(temp) - .tetens(.dewpoint(temp, rh)))}
#Using avg daily temp for simplicity since the package already uses this variable.
#In future implementation, min and max temp could be used to improve accuracy.
#SEE https://www.fao.org/4/x0490e/x0490e07.htm FOR MORE DETAILS


## NEW AS OF 7 JULY 2024  -  FOR CULEX PIPIENS FROM https://doi.org/10.7554%2FeLife.58511
#CALCULATES MOSQ biting RATE COMPONENT DEPENDENT ON TEMPERATURE IN TIME
.mvse_temp_effect_aV <- function(temp){
  ef <- 1.7*10^(-4)*temp*(temp - 9.4)*(39.6-temp)^(1/2)
  ef[which(is.null(ef>39.6))] <- 0 # Setting biting rate to 0 past max threshold
  ef[which(temp>39.6)] <- 0 # Setting biting rate to 0 past max threshold
  ef[which(temp<9.4)] <- 0 # Setting biting rate to 0 past min threshold
  return(ef)
}
#####