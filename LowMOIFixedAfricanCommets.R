##Optional: Reading input
#argu <- as.numeric(commandArgs(trailingOnly = TRUE))
#ind <- argu+1
##Loading the packages##
library(deSolve)
library(minpack.lm)
DVdir <- "ZikaModelParameterOptimization/ZikaKinetics/" ##Direcoty where data files are##

##Defining the model##
fix_delay_ms <- function(time, y,parm){
  with(as.list(c(parm,y)),{
    dS <- k1*S*(1-(S+I1+I2)/N) - k2*S*V
    tlag <-time- tau
    if(tlag <=0){
      Ilag = 0}
    else{
      Ilag = lagvalue(tlag,1)*lagvalue(tlag,4)}
    dI1 <- k2*S*V - k2*Ilag
    dI2 <- k2*Ilag -k4*I2
    dV <- k5*I2 - k6*V- k2*S*V
    return(list(c(dS,dI1,dI2,dV)))
  })
}

##Defineing the function ssq_af to integrate and evaluate the prediction error of parameter set p##
ssq_af <- function(p){
  out <- ode(y0,t,fix_delay_ms,p)
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
t <- c(seq(0,150,length=300),virus2$Time)
t <- sort(unique(t))
##Defining initial conditions##
y0=c(S=3.5e5,I1=0,I2=0,V=5.32e3)

put2 <- NULL#Defining a dataframe/vector to store the prediction error#

##Defining parameters with input or by for loops##
k1 = 3.11e-2
N = 1.94e6
k6=0.05089


#for (j in 1:2){
  #for (l in 1:2){
    parms<-c(k2=1.2e-7,k4=1.0e-2,k5=5.4,tau=15.5)
    fit_af_msfd <- nls.lm(parms,fn=ssq_af, lower = c (1e-9,0,0.1,1),control = nls.lm.control(maxiter = 1000))#Nonlinear Least Squares to minimize the prediction errors#
    ##Integration and calculation of the prediction error##
    out <- ode(y0,virus2$Time,fix_delay_ms,fit_af_msfd$par)
    prd <- data.frame(out)
    prd <- prd[prd$time %in% virus2$Time,]
    res <- (prd$V-virus2$VR.Avg)/virus2$VR.Sdv
    prds<-prd[prd$time>9&prd$time<68,]
    sc = prds$S+prds$I1+prds$I2
    pi = (prds$I1+prds$I2)/sc*100
    pp = prds$I2/sc*100
    cp <- c(as.numeric(!(pi>(spc$African_ExpressingPercentage_average+spc$African_ExpressingPercentage_average))),as.numeric(!(pp<(spc$African_ExpressingPercentage_average-spc$African_ExpressingPercentage_average))))
    put2<-rbind(put2,c(parms,fit_af_msfd$par,res%*%res,cp%*%cp))
  #}
#}
##The parameters and the evaluation as the output##
write.csv(put2,file = paste0(DVdir,'MSFD/African.csv'))
#write.csv(put2,file = paste0(DVdir,'MSFD/Africann',ind,'.csv'))
