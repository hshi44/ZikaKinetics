##Optional: Reading input
#argu <- as.numeric(commandArgs(trailingOnly = TRUE))

##Loading the packages##
library(deSolve)
library(minpack.lm)
#DVdir <- "ZikaModelParameterOptimization/ZikaKinetics/" ##Replace "ZikaModelParameterOptimization" with pathname of the directory##
DVdir <-""

##Defining the model##
gamma_test <- function(time, y, parms,nv){
  with(as.list(c(parms,y)),{
    der<-vector(mode = "numeric",length = length(y))
    
    der[1] <- sum(y[(2+ne):length(y)])*k5 - k6* y[1]
    der[2] <- - ne/taue*y[2]
    #der <-c(dS,dV,dE1)
    if (ne>1){
      for (i in 2:ne){
        der[1+i]<- ne/taue*(y[i]-y[i+1])
      }
    }
    der[2+ne] <-ne/taue*y[1+ne]-ni/taui*y[2+ne]
    if (ni>1){
      for (i in 2:ni){
        der[1+ne+i]<- ni/taui*(y[ne+i]-y[1+ne+i])
      }
    }
    return(list(der))
  })
}

##Defineing the function ssq to integrate and evaluate the prediction error of parameter set p##
ssq_gd <- function(p){
  out <- ode(y0,t,gamma_test,p)
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

##Defining a dataframe/vector to store the prediction error##
put <- NULL

##Measured parameter(s)##

k6=0.06470

##Defineing ne and ni. Different ne's and ni's can be tested by external input or for loops##
#n_grid<-expand.grid(list(20:25,c(5*(10:20),125,200,500,1000)))
#ne<-n_grid[argu+1,1]
#ni<-n_grid[argu+1,2]
ne<-37
ni<-1

E<- replicate(ne,0)
I<- replicate(ni,0)
for (i in 1:ne){
  names(E)[i]<-paste("E", i, sep = "")
}
for (i in 1:ni){
  names(I)[i]<-paste("I", i, sep = "")
}

##Defining initial conditions##
nv<-c(ne=ne,ni=ni)
E[1]=4.28e5
y0=c(V = 0.1, E,I)


##Parameter fitting. Different strating values can be set by external input or for loops##
#for(j in 1:10){
#  for(l in 1:10){
    parms<-c(taue=12.0, taui=17,k5=4.5)
    fit_as_gd_os <- nls.lm(par=parms,fn=ssq_gd, lower = c (1e-4,1e-4,1e-5),control = nls.lm.control(maxiter = 500))#Nonlinear Least Squares to minimize the prediction errors#
    ##Integration and calculation of the prediction error##
    out <- ode(y0,virus$time,gamma_test,fit_as_gd_os$par)
    prd <- data.frame(out)
    res <- prd$V / virus$V-1
    put <- rbind(put,c(parms,fit_as_gd_os$par,nv,res%*%res))
#  }
#}

##The parameters and the evaluation as the output##
write.csv(put,file = paste0(DVdir,'OSGD_Asian.csv'))