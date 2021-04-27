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
gamma_delay_model <- function(time, y, parms,nv){
  with(as.list(c(parms,y)),{
    der<-vector(mode = "numeric",length = length(y))
    C <-unname(sum(y)-y[2])
    der[1] <- k1*y[1]*(1-(C)/N) - k2*y[1]*y[2]
    der[3] <- k2*y[1]*y[2] - ne/taue*y[3]
    #der <-c(dS,dV,dE1)
    if (ne>1){
      for (i in 2:ne){
        der[2+i]<- ne/taue*(y[i+1]-y[i+2])
      }
    }
    der[3+ne] <-ne/taue*y[2+ne]-ni/taui*y[3+ne]
    for (i in 2:ni){
      der[2+ne+i]<- ni/taui*(y[1+ne+i]-y[2+ne+i])
    }
    der[2] <- sum(y[(3+ne):length(y)])*k5 - k6* y[2]-k2*y[1]*y[2]
    return(list(der))
  })
}

##Defining the function ssq_msas to integrate and evaluate the prediction error of parameter set p##
ssq_msas <- function(p){
  t <- c(seq(0,150,length=300),virus2$Time)
  t <- sort(unique(t))
  out <- ode(y0,t,gamma_delay_model,p)
  prd <- data.frame(out)
  prd <- prd[prd$time %in% virus2$Time,]
  res <- (prd$V-virus2$PR.Avg)/virus2$PR.Sdv
  ##Optional: Comparing the model predicted ratio of infected cells with experimental measurements##
  #prds<-prd[prd$time>9&prd$time<68,]
  #sc = prds$S+prds$I1+prds$I2
  #pi = (prds$I1+prds$I2)/sc*100
  #pp = prds$I2/sc*100
  #cp <- c(as.numeric(!(pi>(spc$Asian_ExpressingPercentage_average+spc$Asian_ExpressingPercentage_average))),as.numeric(!(pp<(spc$Asian_ExpressingPercentage_average-spc$Asian_ExpressingPercentage_average))))
  #res <-c(res,5*cp)
  return(res)
}


##Reading experimental data from file(s)##
virus2 <- read.csv(paste0(DVdir,'CountsforLowMOIModel.csv'))
spc <- read.csv(paste0(DVdir,'LowMOIStainingPercentage.csv'))
##Time rage of the integration##
t <- c(seq(0,150,length=300),virus2$Time)
t <- sort(unique(t))

##Measured parameters##
k1 = 3.11e-2
N = 1.94e6
k6=0.07045

##Defining ne and ni. Different ne's and ni's can be tested by external input or for loops##
#n_grid<-expand.grid(list(seq(40,60, by=2),1:5))
#ne<-n_grid[ind,1]
#ni<-n_grid[ind,2]
ne<-50
ni<-3
E<- replicate(ne,0)
I<- replicate(ni,0)
for (i in 1:ne){
  names(E)[i]<-paste("E", i, sep = "")
}
for (i in 1:ni){
  names(I)[i]<-paste("I", i, sep = "")
}
nv<-c(ne=ne,ni=ni)

put2 <- NULL#Defining a dataframe/vector to store the prediction error#

##Defining initial conditions##
y0=c(S = 3.5e5, V = 682, E,I)

##Parameter fitting. Different starting values can be set by external input or for loops##
#for (j in 1:2){
  #for (l in 1:2){
    parms <-c(k2 = 4e-6,taue=7.2,taui=100 ,k5 =0.68)
    fit_as_gd <- nls.lm(par=parms,fn=ssq_msas, lower = c (1e-9,0.5,0.5,1e-4),control = nls.lm.control(maxiter = 1000))#Nonlinear Least Squares to minimize the prediction errors#
    ##Integration and calculation of the prediction error##
    out <- ode(y0,t,gamma_delay_model,fit_as_gd$par)
    prd <- data.frame(out)
    prd <- prd[prd$time %in% virus2$Time,]
    res <- (prd$V-virus2$PR.Avg)/virus2$PR.Sdv
    prds<-prd[prd$time>9&prd$time<68,]
    sc = rowSums(prds[,c(2,4:ncol(prds))])
    pi = rowSums(prds[,4:ncol(prds)])/sc*100
    pp = rowSums(prds[,(4+ne):ncol(prds)])/sc*100
    cp <- c(as.numeric(!(pi>(spc$Asian_ExpressingPercentage_average+spc$Asian_ExpressingPercentage_average))),as.numeric(!(pp<(spc$Asian_ExpressingPercentage_average-spc$Asian_ExpressingPercentage_average))))
    put2 <- rbind(put2,c(parms,fit_as_gd$par,nv,res%*%res,cp%*%cp))
  #}
#}

##Plotting the predicted virus growth curve with optimized parameters. The dots represent experimental averages. The y-axis is on the log scale. Write the output virus and cell kinetics as a file. Comment these lines out if the script is used for parameter optimization only.##
prd <- data.frame(out)
plot(prd$time,prd$V,type = 'l',log="y",xlab = 'Time (h.p.i.)',ylab = 'Virus Titer (PFU/mL)',ylim = c(1e1,1e8))
points(virus2$Time,virus2$PR.Avg)
write.csv(put2,file = paste0(DVdir,'MSGD_Asian_modeloutput.csv'))

##The parameters and the evaluation as the output##
write.csv(put2,file = paste0(DVdir,'MSGD_Asian_parameters.csv'))
#write.csv(put2,file = paste0(DVdir,'MSGD/Asiann',ind,'.csv'))
