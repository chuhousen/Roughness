
# R script for simultaneous estimation of roughness length (Z0) & zero-plane displacement height (d) 
#   from single-level eddy covariance measurements

# The code modifies and implements the 3 approaches used in Graf et al.,2014 Boundary-Layer Meteorol.
#   FP-RE-1: from logarithmic wind profile theory: WS/Ustr = 1/k*(ln((zm-d)/Z0)-beta*((zm-d)/Lm))

#   FP-RE-2: from logarithmic wind profile theory: WS = Ustr/k*(ln((zm-d)/Z0)-beta*((zm-d)/Lm))

#   FV-RE-1: from flux-variance similarity theory: std.Uz/Ustr = C1*(1-C2*(zm-d)/Lm)^(1/3)
#                                                  std.Uz = k*C1*WS/ln((zm-d)/Z0)                 

# [Desccription]
#   The goal is to estimate Z0 and d from continuous measurements of wind and turbulent statistics, i.e., WS, Ustr, std.Uz..... 
#   All input variables are originally measured at a 30-min time step, but are pre-filtered according to data quality & required meteorological conditions.
#   i.e., not every 30 min is used in model fitting.
#
#   While estimating Z0 & d, the time series is re-grouped to a daily step, i.e., treating each 30 min within a day as a sampling replicate.
#   Z0 and d are estimated at a daily time step, as they are both functions of canopy structures (height, leaf area...), assuming no changes within a day   
#   The time series of Z0 or d is assumed as an autoregressive process, and it's assumed there's correlation (cor) between the time series of Z0 & d
#   Z0[t+1] = Z0[t] + Z0.err ; Z0.err ~ N(0,sigma.Z0)
#   d[t+1] = d[t] + sigma.d/sigma.Z0*cor*(Z0[t+1]-Z0[t]) + d.err ; d.err ~ N(0,sigma.d.iid)
#
# [Input variables]  
#       WS: wind speed (m/s)
#       Ustr: friction velocity (m/s)
#       zm: measurement height (m)
#       Lm: Monin-Obukhov length (m)
#       std.Uz: standard deviation of vertical wind velocity (m/s)
#       hc: vegetation canopy height (m) *only used in filtering data*
#       ZL: Monin-Obukhov stability (unitless) *only used in filtering data*
#
# [Constants]      
#       k: 0.4 von karman constant
#       C1: 1.3
#       C2: 2.0 estimates for the universal constants, (Panofsky & Dutton 1984)
#       beta: 6.0
#
# [Output variables]
#       Z01[DOY]: roughness length estimated for each day (DOY) by using FP-RE-1 model
#       d1[DOY]: displacement height estimated for each day (DOY) by using FP-RE-1 model
#
#       Z02[DOY]: roughness length estimated for each day (DOY) by using FP-RE-2 model
#       d2[DOY]: displacement height estimated for each day (DOY) by using FP-RE-2 model
#
#       Z03[DOY]: roughness length estimated for each day (DOY) by using FV-RE-1 model
#       d3[DOY]: displacement height estimated for each day (DOY) by using FV-RE-1 model


rm(list=ls())
require(snow)
require(rjags)
require(dclone)
require(R2WinBUGS)
require(coda)

## control variables 
ver<-"4.5"
case<-"CL_2011"
output<-T

n.chains <- 3
n.adapt <- 5000
n.update <- 10000
n.iter <- 5000
n.keep <- 5000
thin <- max(c(1,floor((n.iter)*n.chains/n.keep)))

## Input data 
path<-"D:\\Housen\\Flux\\Data-exploring\\00_Synthesis_Fluxnet_Roughness\\"
setwd(paste(path,"Output\\",sep=""))

data.in<-read.table(paste(path,"Input\\2016-02-09_srt_",case,".csv",sep=""),
                    header=T,sep=",", na.string=c("-9999"))
data.in<-data.in[!is.na(data.in$Time.id),]

## pre-filtering data with thresholds 
data.in$Z0.max.L<-0.1*data.in$hc/data.in$Lm  
data.in$Z.max.L<-(data.in$zm-0.66*data.in$hc)/data.in$Lm
data.in[!is.na(data.in$Z0.max.L)&data.in$Z0.max.L>0.037,c("Ustr","WS")]<-NA
data.in[!is.na(data.in$Z0.max.L)&data.in$Z0.max.L<(-0.084),c("Ustr","WS")]<-NA
data.in[!is.na(data.in$WS)&data.in$WS>(5),c("std_Uz")]<-NA
data.in[!is.na(data.in$Ustr)&data.in$Ustr<(0.05),c("std_Uz")]<-NA
data.in[!is.na(data.in$Z.max.L)&data.in$Z.max.L<(-0.05),c("std_Uz")]<-NA
data.in[!is.na(data.in$Z.max.L)&data.in$Z.max.L>0.05,c("std_Uz")]<-NA
data.in<-na.omit(data.in)

## Model specification for JAGS
bugs.in<-function(infile=data.in){
  
  Ustr <- infile$Ustr
  WS <- infile$WS
  Lm <- infile$Lm
  std.Uz <- infile$std_Uz
  hc <- infile$hc
  zm <- infile$zm
  ZL <- infile$Z.max.L
  UzUstr <- infile$std_Uz/infile$Ustr
  WSUstr <- infile$WS/infile$Ustr
  DOY <- infile$DOY
  n <- dim(infile)[1] # length of input 30-min data
  n2 <- max(infile$DOY,na.rm=T)-min(infile$DOY,na.rm=T)+1  # number of days to be estimated
  n2.i <- min(infile$DOY,na.rm=T)  # the first day to be estimated 
  k <- 0.4
  beta <- 6.0
  C1 <- 1.3
  C2 <- 2.0
  max.zm <- max(infile$zm,na.rm=T)  # max measurement height, used as upper bound of d & Z0 estimation ( 0 < d < zm; 0 < Z0 < 0.1*zm )
  rng.WS <- max(infile$WS,na.rm=T)-min(infile$WS,na.rm=T)  # range of WS
  rng.WSUstr <- max(WSUstr,na.rm=T)-min(WSUstr,na.rm=T)    # range of WS/Ustr
  rng.UzUstr <- max(UzUstr,na.rm=T)-min(UzUstr,na.rm=T)    # range of std.Uz/Ustr
  rng.stdUz <- max(std.Uz,na.rm=T)-min(std.Uz,na.rm=T)     # range of std.Uz
  
  bugs.data <- list(Ustr=Ustr,WS=WS,Lm=Lm,std.Uz=std.Uz,hc=hc,zm=zm,UzUstr=UzUstr,WSUstr=WSUstr,DOY=DOY,
                    n=n,n2=n2,n2.i=n2.i,k=k,beta=beta,C1=C1,C2=C2,
                    max.zm=max.zm,rng.WS=rng.WS,rng.WSUstr=rng.WSUstr,rng.UzUstr=rng.UzUstr,rng.stdUz=rng.stdUz)

  bugs.ini <- list()
  for (i in 1:n.chains){
    bugs.ini[[i]] <- list(d1=runif(n2,0.001,max.zm),
                          d2=runif(n2,0.001,max.zm),
                          d3=runif(n2,0.001,max.zm),
                          Z01=runif(n2,0.001,max.zm*0.1),
                          Z02=runif(n2,0.001,max.zm*0.1),
                          Z03=runif(n2,0.001,max.zm*0.1),
                          cor1=runif(1,-1,1),
                          cor2=runif(1,-1,1),
                          cor3=runif(1,-1,1),
                          sigma1=runif(1,0.001,rng.WSUstr),
                          sigma2=runif(1,0.001,rng.WS),
                          sigma3=runif(1,0.001,rng.UzUstr),
                          sigma4=runif(1,0.001,rng.stdUz)
                          )
  }
  
  bugs.para <- c("d1","d2","d3","Z01","Z02","Z03",
                 "sigma.d1","sigma.d2","sigma.d3",
                 "sigma.Z01","sigma.Z02","sigma.Z03",
                 "cor1","cor2","cor3")
  
  return(list(data=bugs.data, 
              inits=bugs.ini, 
              para=bugs.para,
              n.chains=n.chains, 
              model=function(){
                for(i in 1:n){
                  ## FP-RE-1 model
                  WSUstr[i] ~ dnorm(mu1[i],pow(sigma1,-2));
                  mu1[i] <- log((zm[i]-d1[DOY[i]-n2.i+1])/Z01[DOY[i]-n2.i+1])/k+beta*(zm[i]-d1[DOY[i]-n2.i+1])/Lm[i]/k;
                  ## FP-RE-2 model
                  WS[i] ~ dnorm(mu2[i],pow(sigma2,-2));
                  mu2[i] <- log((zm[i]-d2[DOY[i]-n2.i+1])/Z02[DOY[i]-n2.i+1])/k*Ustr[i]+beta*(zm[i]-d2[DOY[i]-n2.i+1])/Lm[i]/k*Ustr[i];
                  ## FV-RE-1 model
                  UzUstr[i] ~ dnorm(mu3[i],pow(sigma3,-2));
                  mu3[i] <- C1*pow((1-C2*((zm[i]-d3[DOY[i]-n2.i+1])/Lm[i])),1/3);
                  std.Uz[i] ~ dnorm(mu4[i],pow(sigma4,-2));
                  mu4[i] <- k*C1*WS[i]/log((zm[i]-d3[DOY[i]-n2.i+1])/Z03[DOY[i]-n2.i+1]);
                }
                
                d1[1] ~ dnorm(0.01,pow(sigma.d1,-2));T(0,max.zm)
                sigma.d1 ~ dunif(0.001,max.zm);
                Z01[1] ~ dnorm(0.01,pow(sigma.Z01,-2));T(0,max.zm*0.1)
                sigma.Z01 ~ dunif(0.001,max.zm*0.1);
                cor1 ~ dunif(-1,1);
                tau.d1.iid <- 1/((1-cor1*cor1)*sigma.d1*sigma.d1);
                
                d2[1] ~ dnorm(0.01,pow(sigma.d2,-2));T(0,max.zm)
                sigma.d2 ~ dunif(0.001,max.zm);
                Z02[1] ~ dnorm(0.01,pow(sigma.Z02,-2));T(0,max.zm*0.1)
                sigma.Z02 ~ dunif(0.001,max.zm*0.1);
                cor2 ~ dunif(-1,1);
                tau.d2.iid <- 1/((1-cor2*cor2)*sigma.d2*sigma.d2);
                
                d3[1] ~ dnorm(0.01,pow(sigma.d3,-2));T(0,max.zm)
                sigma.d3 ~ dunif(0.001,max.zm);
                Z03[1] ~ dnorm(0.01,pow(sigma.Z03,-2));T(0,max.zm*0.1)
                sigma.Z03 ~ dunif(0.001,max.zm*0.1);
                cor3 ~ dunif(-1,1);
                tau.d3.iid <- 1/((1-cor3*cor3)*sigma.d3*sigma.d3);
                
                for (j in 2:n2){
                  Z01[j] ~ dnorm(Z01[j-1],pow(sigma.Z01,-2));T(0,max.zm*0.1)
                  d1[j] ~ dnorm(d1[j-1]+sigma.d1/sigma.Z01*cor1*(Z01[j]-Z01[j-1]),tau.d1.iid);T(0,max.zm)
                  Z02[j] ~ dnorm(Z02[j-1],pow(sigma.Z02,-2));T(0,max.zm*0.1)
                  d2[j] ~ dnorm(d2[j-1]+sigma.d2/sigma.Z02*cor2*(Z02[j]-Z02[j-1]),tau.d2.iid);T(0,max.zm)
                  Z03[j] ~ dnorm(Z03[j-1],pow(sigma.Z03,-2));T(0,max.zm*0.1)
                  d3[j] ~ dnorm(d3[j-1]+sigma.d3/sigma.Z03*cor3*(Z03[j]-Z03[j-1]),tau.d3.iid);T(0,max.zm)
                }
                
                sigma1 ~ dunif(0.001,rng.WSUstr);
                sigma2 ~ dunif(0.001,rng.WS);
                sigma3 ~ dunif(0.001,rng.UzUstr);
                sigma4 ~ dunif(0.001,rng.stdUz);
                
              }
              )
         )
}


### main workflow -- excute JAGS
input.to.bugs <- bugs.in()

date()
load.module("lecuyer")
cl <- makeCluster(n.chains,type="SOCK")
parLoadModule(cl,"lecuyer")
m2 <- jags.parfit(cl,data=input.to.bugs$data, 
                  params=input.to.bugs$para, 
                  model=input.to.bugs$model,
                  inits=input.to.bugs$inits,
                  n.adapt=n.adapt,
                  n.update=n.update,
                  n.iter=n.iter,
                  thin=thin,
                  n.chains=n.chains)
#nchain(m2)
unload.module("lecuyer")
parUnloadModule(cl,"lecuyer")
stopCluster(cl)
date()


## prepare output
jags.out <- NULL
for(i in 1:n.chains)
  jags.out <- rbind(jags.out,m2[[i]])

m2.sum<-summary(m2)
#para.scale<-gelman.diag(m2)
para.summary<-(cbind(m2.sum$statistics,m2.sum$quantiles))#,para.scale$psrf))

if(output){
  ## dignostic plot
  pdf(paste(getwd(),"/",substr(as.character(Sys.time()),1,10),"_",case,"_Z0d_",ver,"_diag.pdf",sep=""),
      width=9,height=9)
  plot(m2);
  dev.off()
  
  ## prediction time series plot
  pdf(paste(getwd(),"/",substr(as.character(Sys.time()),1,10),"_",case,"_Z0d_",ver,"_para.pdf",sep=""),
      width=8,height=10)
  n2 <- max(data.in$DOY)-min(data.in$DOY)+1
  n2.i <- min(data.in$DOY)
  par(mfrow=c(3,1),mar=c(4,4,0.5,0.5))
  plot(c(n2.i:(n2.i+n2-1)),para.summary[1:n2,7],col="red",type="l",lwd=3.5,ylim=c(0,0.2),ylab="Z0 (m)",xlab="DOY")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[1:n2,5],col="pink")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[1:n2,9],col="pink")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(n2+1):(2*n2),7],col="blue",lwd=3.5)
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(n2+1):(2*n2),5],col="lightblue")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(n2+1):(2*n2),9],col="lightblue")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(2*n2+1):(3*n2),7],col="forestgreen",lwd=3.5)
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(2*n2+1):(3*n2),5],col="green")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(2*n2+1):(3*n2),9],col="green")
  
  plot(c(n2.i:(n2.i+n2-1)),para.summary[(3*n2+4):(4*n2+3),7],col="red",type="l",lwd=3.5,ylim=c(0,2),ylab="d (m)",xlab="DOY")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(3*n2+4):(4*n2+3),5],col="pink")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(3*n2+4):(4*n2+3),9],col="pink")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(4*n2+4):(5*n2+3),7],col="blue",lwd=3.5)
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(4*n2+4):(5*n2+3),5],col="lightblue")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(4*n2+4):(5*n2+3),9],col="lightblue")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(5*n2+4):(6*n2+3),7],col="forestgreen",lwd=3.5)
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(5*n2+4):(6*n2+3),5],col="green")
  lines(c(n2.i:(n2.i+n2-1)),para.summary[(5*n2+4):(6*n2+3),9],col="green")
  
  plot(tapply(data.in$DOY,data.in$DOY,mean),
       tapply(data.in$DOY,data.in$DOY,length),type="s",lwd=1,ylim=c(0,48),ylab="Sample (n/day)",xlab="DOY")
  dev.off()
  
  ## model setup and posterior outputs
  write.table(jags.out,paste(getwd(),"/",substr(as.character(Sys.time()),1,10),"_",case,"_Z0d_",ver,"_para.csv",sep=""),sep=",",row.names=F)  
  write.table(para.summary,paste(getwd(),"/",substr(as.character(Sys.time()),1,10),"_",case,"_Z0d_",ver,"_para_sum_.csv",sep=""),sep=",",row.names=F)
  write.table(rownames(para.summary),paste(getwd(),"/",substr(as.character(Sys.time()),1,10),"_",case,"_Z0d_",ver,"_para_sum_name.txt",sep=""),sep="\t")
  sink((paste(getwd(),"/",substr(as.character(Sys.time()),1,10),"_",case,"_Z0d_",ver,"_mdl.txt",sep="")))
  lapply(input.to.bugs$model,print) 
  sink()
  sink((paste(getwd(),"/",substr(as.character(Sys.time()),1,10),"_",case,"_Z0d_",ver,"_init.txt",sep="")))
  lapply(input.to.bugs$inits,print) 
  sink()
}

