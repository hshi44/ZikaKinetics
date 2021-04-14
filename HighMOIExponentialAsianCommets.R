##Optional: Reading input
#argu <- as.numeric(commandArgs(trailingOnly = TRUE))
#ind <- argu+1
##Loading the packages##
library(deSolve)
library(minpack.lm)
#DVdir <- "ZikaModelParameterOptimization/ZikaKinetics/" ##Direcoty where data files are. Replace "ZikaModelParameterOptimization" with the pathname.##
DVdir <- ""

##Defining the model##
exp_delay_simple_model <- function(time, y, parms){
  with(as.list(c(parms,y)),{
    dI1 <-  - k3*I1
    dI2 <- k3*I1 - k4*I2
    dV <- k5*I2 - k6*V
    return(list(c(dI1,dI2,dV)))
  })
}

##Defineing the function ssq to integrate and evaluate the prediction error of parameter set p##
ssq_ed <- function(p){
  out <- ode(y0,t,exp_delay_simple_model,p)
  prd <- data.frame(out)
  prd <- prd[prd$time %in% virus$time,]
  res <- prd$V / virus$V-1
  return(res)
}

##Reading experimetnal data from file(s)##
counts <- read.csv(paste0(DVdir,'CountsforHighMOIModel.csv'))
virus <- data.frame(counts[,c(1,3)])
colnames(virus) <- c("time","V")
##Time rage of the intergration##
t <- c(seq(0,64,length=200),virus$time)
t <- sort(unique(t))
##Defining initial conditions##
y0 <- c(I1=4.28e5, I2=0, V=0.1)

##Defining a dataframe/vector to store the prediction error##
put <- NULL

##Measured parameter(s)##

k6=0.06470

##Parameter fitting. Different strating values can be set by external input or for loops##
#for(j in 1:10){
#  for(l in 1:10){
    parms<-c(k3=5.6e-5,k4=1e-8,k5=1.3)
    fit_as_ed <- nls.lm(par=parms,fn=ssq_ed,lower = c (1e-8,1e-9,0.1),control = nls.lm.control(maxiter = 300))#Nonlinear Least Squares to minimize the prediction errors#
    ##Integration and calculation of the prediction error##
    out <- ode(y0,virus$time,exp_delay_simple_model,fit_as_ed$par)
    prd <- data.frame(out)
    res <- prd$V / virus$V-1
    put <- rbind(put, c(parms,fit_as_ed$par,res%*%res))
#  }
#}

##The parameters and the evaluation as the output##
write.csv(put,file = paste0(DVdir,'OSED_Asian.csv'))