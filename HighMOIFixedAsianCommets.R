##Optional: Reading input
#argu <- as.numeric(commandArgs(trailingOnly = TRUE))
#ind <- argu+1
##Loading the packages##
library(deSolve)
library(minpack.lm)
#DVdir <- "ZikaModelParameterOptimization/ZikaKinetics/" ##Replace "ZikaModelParameterOptimization" with pathname of the directory##
DVdir <- ""

##Defining the model##
fixed_delay_model <- function(time, y, parms){
  with(as.list(c(parms,y)),{
    tlag <- time- tau
    if(tlag <0){
      dI = 0
      dV = - k6*V
    }
    
    else{
      dI = -k4*I
      dV <- k5*I - k6*V
    }
    return(list(c(dI,dV),tlag))
  })
}

##Defineing the function ssq to integrate and evaluate the prediction error of parameter set p##
ssq <- function(p){
  out <- dede(y0,t,fixed_delay_model,p)
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
y0 <- c(I=4.28e5, V=0)

##Defining a dataframe/vector to store the prediction error##
put <- NULL

##Measured parameter(s)##

k6=0.06470

##Parameter fitting. Different strating values can be set by external input or for loops##
#for(j in 1:10){
#  for(l in 1:10){
    parms<-c(k4=5.4e-2,k5=4.2,tau=11.3)
    fit_as_fd_new <- nls.lm(par=parms,fn=ssq,lower = c (1e-9,1e-4,1e-2),control = nls.lm.control(maxiter = 300))#Nonlinear Least Squares to minimize the prediction errors#
    ##Integration and calculation of the prediction error##
    out <- ode(y0,virus$time,fixed_delay_model,fit_as_fd_new$par)
    prd <- data.frame(out)
    res <- prd$V / virus$V-1
    put <- rbind(put,c(parms,fit_as_fd_new$par,res%*%res))
#  }
#}

##The parameters and the evaluation as the output##
write.csv(put,file = paste0(DVdir,'OSFD_Asian.csv'))