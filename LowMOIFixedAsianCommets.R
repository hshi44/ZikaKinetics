##Optional: Reading input
#argu <- as.numeric(commandArgs(trailingOnly = TRUE))
#ind <- argu+1
##Loading the packages##
library(deSolve)
library(minpack.lm)
##Direcoty where data files are##
#DVdir <- "ZikaModelParameterOptimization/ZikaKinetics/" 
DVdir <-""

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
ssq_as <- function(p){
  out <- dede(y0,t,fix_delay_ms,p)
  prd <- data.frame(out)
  prd <- prd[prd$time %in% virus2$Time,]
  res <- (prd$V-virus2$PR.Avg)/virus2$PR.Sdv
  ##Optional: Comparing the model predicted ratio of infected cells with experimetal measurements##
  #prds<-prd[prd$time>9&prd$time<68,]
  #sc = prds$S+prds$I1+prds$I2
  #pi = (prds$I1+prds$I2)/sc*100
  #pp = prds$I2/sc*100
  #cp <- c(as.numeric(!(pi>(spc$Asian_ExpressingPercentage_average+spc$Asian_ExpressingPercentage_average))),as.numeric(!(pp<(spc$Asian_ExpressingPercentage_average-spc$Asian_ExpressingPercentage_average))))
  #res <-c(res,5*cp)
  return(res)
}


##Reading experimetnal data from file(s)##
virus2 <- read.csv(paste0(DVdir,'CountsforLowMOIModel.csv'))
spc <- read.csv(paste0(DVdir,'LowMOIStainingPercentage.csv'))
##Time rage of the intergration##
t <- c(seq(0,150,length=300),virus2$Time)
t <- sort(unique(t))
##Defining a dataframe/vector to store the prediction error##
put2 <- NULL

##Defining initial conditions##
y0=c(S=3.5e5,I1=0,I2=0,V=682)

##Measured parameters##
k1 = 3.11e-2
N = 1.94e6
k6=0.06470

##Parameter fitting. Different strating values can be set by external input or for loops. res is the prediction error.##
#for (j in 1:2){
  #for (l in 1:2){
    parms<-c(k2=6.1e-6,k4=1.6e-2,k5=0.94,tau=8.5)
    fit_as_msfd <- nls.lm(parms,fn=ssq_as, lower = c (1e-9,1e-9,0.01,1),control = nls.lm.control(maxiter = 1000))#Nonlinear Least Squares to minimize the prediction errors#
    ##Integration and calculation of the prediction error##
    out <- ode(y0,t,fix_delay_ms,fit_as_msfd$par)
    prd <- data.frame(out)
    prd <- prd[prd$time %in% virus2$Time,]
    res <- (prd$V-virus2$PR.Avg)/virus2$PR.Sdv
    prds<-prd[prd$time>9&prd$time<68,]
    sc = prds$S+prds$I1+prds$I2
    pi = (prds$I1+prds$I2)/sc*100
    pp = prds$I2/sc*100
    cp <- c(as.numeric(!(pi>(spc$Asian_ExpressingPercentage_average+spc$Asian_ExpressingPercentage_average))),as.numeric(!(pp<(spc$Asian_ExpressingPercentage_average-spc$Asian_ExpressingPercentage_average))))
    put2<-rbind(put2,c(parms,fit_as_msfd$par,res%*%res,cp%*%cp))
  #}
#}
##The parameters and the evaluation as the output##
write.csv(put2,file = paste0(DVdir,'MSFD_Asian.csv'))
#write.csv(put2,file = paste0(DVdir,'MSFD/Asiann',ind,'.csv'))
