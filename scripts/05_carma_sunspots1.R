##########################
### A CARMA(2,1) model ###
##########################
library(ctsem)
#get data
sunspots<-sunspot.year
sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
id <- 1
time <- 1749:1924
datalong <- cbind(id, time, sunspots)
#setup model
ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1,
                   manifestNames='sunspots',
                   latentNames=c('ss_level', 'ss_velocity'),
                   LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
                   DRIFT=matrix(c(0, 'a21', 1, 'a22'), nrow=2, ncol=2),
                   MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
                   CINT=matrix(c(0, 0), nrow=2, ncol=1),
                   T0VAR=matrix(c(1,0,0,1), nrow=2, ncol=2), #Because single subject
                   DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))
ssmodel$pars$indvarying<-FALSE #Because single subject
ssmodel$pars$transform[14]<- '(param)*5+44 ' #Because not mean centered
ssmodel$pars$transform[4]<-'log(exp(-param*1.5)+1)' #To avoid multimodality
#fit
ssfit <- ctStanFit(datalong, ssmodel, iter=1000, chains=4)
#output
summary(ssfit, parmatrices = TRUE)
plot(ssfit)



