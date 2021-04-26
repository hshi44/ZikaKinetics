##Optional: Reading input
#argu <- as.numeric(commandArgs(trailingOnly = TRUE))
#ind <- argu+1
##Loading the packages##
library(deSolve)
library(minpack.lm)
#DVdir <- "ZikaModelParameterOptimization/ZikaKinetics/" ##Replace "ZikaModelParameterOptimization" with pathname of the directory##

##Defining the model##
exp_delay_simple_model <- function(time, y, parms){
  with(as.list(c(parms,y)),{
    dI1 <-  - k3*I1
    dI2 <- k3*I1 - k4*I2
    dV <- k5*I2 - k6*V
    return(list(c(dI1,dI2,dV)))
  })
}

##Defining the function ssq to integrate and evaluate the prediction error of parameter set p##
ssq_ed <- function(p){
  out <- ode(y0,t,exp_delay_simple_model,p)
  prd <- data.frame(out)
  prd <- prd[prd$time %in% virus$time,]
  res <- prd$V / virus$V-1
  return(res)
}

##Reading experimental data from file(s)##
#counts <- read.csv(paste0(DVdir,'CountsforHighMOIModel.csv'))
counts <- read.csv('CountsforHighMOIModel.csv')
virus <- data.frame(counts[,c(1,2)])
colnames(virus) <- c("time","V")
##Time rage of the integration##
t <- c(seq(0,64,length=200),virus$time)
t <- sort(unique(t))
##Defining initial conditions##
y0 <- c(I1=4.28e5, I2=0, V=0.1)

##Defining a dataframe/vector to store the prediction error##
put <- NULL

##Measured parameter(s)##

k6=0.05089

##Parameter fitting. Different starting values can be set by external input or for loops##
#for(j in 1:10){
#  for(l in 1:10){
    parms<-c(k3=1.0e-6,k4=1e-8,k5=50)
    fit_af_ed <- nls.lm(par=parms,fn=ssq_ed,lower = c (1e-8,1e-9,0.1),control = nls.lm.control(maxiter = 300))#Nonlinear Least Squares to minimize the prediction errors#
    ##Integration and calculation of the prediction error##
    out <- ode(y0,virus$time,exp_delay_simple_model,fit_af_ed$par)
    prd <- data.frame(out)
    res <- prd$V / virus$V-1
    put <- rbind(put, c(parms,fit_af_ed$par,res%*%res))
#  }
#}

##Plotting the predicted virus growth curve with optimized parameters. The dots represent experimental data. Write the output virus and cell kinetics as a file. Comment these lines out if the script is used for parameter optimization only.##
out <- ode(y0,t,exp_delay_simple_model,fit_af_ed$par)
prd <- data.frame(out)
plot(prd$time,prd$V,type = 'l',log="y",xlab = 'Time (h.p.i.)',ylab = 'Virus Titer (PFU/mL)',ylim = c(1e1,1e9))
points(virus)
write.csv(prd,file = 'OSED_African_modeloutput.csv')

##The parameters and the evaluation as the output##
#write.csv(put,file = paste0(DVdir,'OSED/African.csv'))
write.csv(put,file = 'OSED_African_parameters.csv')
