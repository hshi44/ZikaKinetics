##Optional: Reading input
#argu <- as.numeric(commandArgs(trailingOnly = TRUE))
#ind <- argu+1
##Loading the packages##
library(deSolve)
library(minpack.lm)
#DVdir <- "ZikaModelParameterOptimization/ZikaKinetics/" ##Direcoty where data files are##
DVdir <-""

##Defining the model##
exp_delay_full_model <- function(time, y, parms){
  with(as.list(c(parms,y)),{
    dS <- k1*S*(1-(S+I1+I2)/N) - k2*S*V
    dI1 <- k2*S*V - k3*I1
    dI2 <- k3*I1 - k4*I2
    dV <- k5*I2 - k6*V- k2*S*V
    return(list(c(dS,dI1,dI2,dV)))
  })
}
##Defineing the function ssq to integrate and evaluate the prediction error of parameter set p##
ssq <- function(p){
  t <- c(seq(0,150,length=300),virus2$Time)
  t <- sort(unique(t))
  out <- ode(y0e,t,exp_delay_full_model,p)
  prd <- data.frame(out)
  prd <- prd[prd$time %in% virus2$Time,]
  res <- (prd$V-virus2$VR.Avg)/virus2$VR.Sdv
  ##Optional: Comparing the model predicted ratio of infected cells with experimetal measurements##
  #prds<-prd[prd$time>9&prd$time<68,]
  #sc = prds$S+prds$I1+prds$I2
  #pi = (prds$I1+prds$I2)/sc*100
  #pp = prds$I2/sc*100
  #cp <- c(as.numeric(!(pi>(spc$African_ExpressingPercentage_average+spc$African_ExpressingPercentage_average))),as.numeric(!(pp<(spc$African_ExpressingPercentage_average-spc$African_ExpressingPercentage_average))))
  #res <-c(res,5*cp)
  return(res)
}

##Reading experimetnal data from file(s)##
virus2 <- read.csv(paste0(DVdir,'CountsforLowMOIModel.csv'))
spc <- read.csv(paste0(DVdir,'LowMOIStainingPercentage.csv'))
##Time rage of the intergration##
t <- seq(0,150,0.1)
##Defining initial conditions##
y0e<-c(S=3.5e5,I1=0,I2=0,V=5.32e3)

put <- NULL#Defining a dataframe/vector to store the prediction error#

##Measured parameters##
k1 = 3.11e-2
N = 1.94e6
k6=0.05089

##Parameter fitting. Different strating values can be set by external input or for loops. res2 is the prediction error.##
#for(j in seq(1.50,1.56,0.002)){
#  for(l in seq(3.4,3.8,0.02)){
    parms<-c(k2=1.6e-7,k3=1.0e-6,k4=1e-5,k5=50)
    fit_af_cm <- nls.lm(parms,fn=ssq, lower = c (1e-10,1e-10,0,1e-4),upper=c(0.1,10,10,1e6),control = nls.lm.control(maxiter = 800))#Nonlinear Least Squares to minimize the prediction errors#
    ##Integration and calculation of the prediction error##
    out <- ode(y0e,t,exp_delay_full_model,fit_af_cm$par)
    prd <- data.frame(out)
    prd <- prd[prd$time %in% virus2$Time,]
    res1 <- prd$V/virus2$VR.Avg-1
    res2 <- (prd$V-virus2$VR.Avg)/virus2$VR.Sdv
    prds<-prd[prd$time>9&prd$time<68,]
    sc = prds$S+prds$I1+prds$I2
    pi = (prds$I1+prds$I2)/sc*100
    pp = prds$I2/sc*100
    cp <- c(as.numeric(!(pi>(spc$African_ExpressingPercentage_average+spc$African_ExpressingPercentage_average))),as.numeric(!(pp<(spc$African_ExpressingPercentage_average-spc$African_ExpressingPercentage_average))))
    put<-rbind(put,c(parms,fit_af_cm$par,res1%*%res1,res2%*%res2,cp%*%cp))
#  }
#}

##The parameters and the evaluation as the output##
write.csv(put,file = paste0(DVdir,'MSED_African.csv'))