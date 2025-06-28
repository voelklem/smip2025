set.seed(1)
library(ctsem)

#Data generation (univariate process, 4 cohorts; T = 50)####

measurementErrorSD <- c(1, .5, .2, .2) #measurement error sd differs between the 4 cohorts
drift <- c(-.3, -.5, -.6, -.6) #drift coefficient differs between the 4 cohorts
                               #note that cohort 3 and 4 do not differ.  

for(cohorti in 1:4){
  
  gm <- ctModel(
    LAMBDA=diag(1), #single observed variable and single process
    Tpoints = 50,
    CINT=5,
    TRAITVAR=diag(1), #include individual differences in continuous intercept
    DRIFT=drift[cohorti], # different temporal dynamics
    DIFFUSION=.5,
    MANIFESTVAR=0) #measurement error sd fixed to zero so we can plot true scores (is added later on)
  
  #when generating data, free pars are set to 0
  d <- data.frame(ctGenerate(ctmodelobj = gm,n.subjects = 50,
                             burnin = 5,dtmean = .1),Cohort2=0,Cohort3=0,Cohort4=0, Cohort=cohorti)
  d$id <- paste0('Cohort',cohorti,'_id',d$id)
  
  d <- d[d$time < (cohorti) & d$time > (cohorti-1),]
  
  d$Yobs <- d$Y1 + rnorm(nrow(d),0,measurementErrorSD[cohorti]) #add measurement error that differs across cohorts
  
  if(cohorti == 1){
    dat <- d
  } else{
    d[[paste0('Cohort',cohorti)]] <- 1 #create three dummy variables Cohort2, Cohort3, Cohort4 (else baseline = cohort 1)
    dat <- rbind(dat,d)
  }
}

#Plotting the data####

library(ggplot2)

ggplot(dat, aes(y=Y1,x=time,colour=factor(Cohort),group=id))+
  geom_line()+theme_bw()+
  geom_line(mapping = aes(y=Yobs),alpha=.3)+
  guides(colour='none')

#setting the starting point (for all individuals) to the earliest observation in the data####


dat0 <- dat[!duplicated(dat$id),] #create extra dataframe containing one observation of each subject
dat0$time <- min(dat$time) #set the times of each observation to the time of earliest observation in our dataset
dat0$Y1 <- dat0$Yobs <- NA #set observed variables to NA (but keep the info for time independent predictors such as group dummies!)
dat <- merge(dat, dat0, all = TRUE)

#define ct model m1 that ignores cohort differences ####

m1 <- ctModel(
  type='stanct',
  CINT='cint', 
  MANIFESTMEANS = 0,
  LAMBDA=diag(1),
  manifestNames = 'Yobs')

f1 <- ctStanFit(datalong = dat, ctstanmodel = m1, cores=2)

#define ct model m2 accounting for cohort differences ####

m2 <- ctModel(
  type='stanct',
  CINT='cint', 
  TIpredNames = c('Cohort2','Cohort3','Cohort4'),
  MANIFESTMEANS = 0,
  LAMBDA=diag(1),
  manifestNames = 'Yobs')

#ctModelLatex(m2) #view the full model as a pdf

f2 <- ctStanFit(datalong = dat, ctstanmodel = m2, cores=2)

knitr::include_graphics('mdefaulttex.png')

plot<-ctKalman(f2, subjects = dat$id[!duplicated(dat$Cohort)], plot=T, 
               kalmanvec=c('yprior'),removeObs = T, polygonsteps=1)
print(plot+facet_wrap(vars(Subject)))

ctKalman(f2, subjects = dat$id[!duplicated(dat$Cohort)], plot=T, 
         kalmanvec=c('y','etasmooth'),facets=NA)

#Excursus: Define model m3 that accounts for cohort differences in random effects by ####
#using the state expansion trick

m3 <- ctModel(
  type='stanct',
  manifestNames = 'Yobs',
  latentNames=c('eta','intercept'), #we now have two processes -- main process eta, and extra process 'intercept'
  CINT=0, #we no longer estimate a parameter here, because we have moved it to the 2nd system process
  DRIFT=c( 
    'drift1', 1, #self feedback parameter, and effect of intercept process on main process fixed to 1
    0, 0), #no effect of main process on intercept, and intercept process does not change so auto effect is 0
  DIFFUSION=c(
    'diffusion1',0, #only the main process is subject to random fluctuations
    0,0), 
  T0VAR=c( #since we are now specifying individual differences by manually extending the system matrices, 
    't0v11',0, #stable individual differences come in via the initial variance (t0var) in starting points. 
    't0v21','t0v22'), #correlation parameters are specified in the lower off diagonal.
  T0MEANS=c(
    't0eta||FALSE', #freely estimate starting points for both processes but,
    't0intercept||FALSE'), # *disable* ctsem's random effects handling (set to FALSE)
  TIpredNames = c('Cohort2','Cohort3','Cohort4'),
  MANIFESTMEANS = 0,
  LAMBDA=matrix(c(1,0),nrow=1, ncol=2)) #intercept process does not directly affect observed variables)

f3 <- ctStanFit(datalong = dat, ctstanmodel = m3, cores=2)

ctKalman(f3, subjects = dat$id[!duplicated(dat$Cohort)],plot=T,kalmanvec='etaprior')

# Comparing models ####

print(data.frame(
  Model=c('m1','m2','m3'),
  AIC=c(summary(f1)$aic, summary(f2)$aic, summary(f3)$aic)))
##   Model      AIC
## 1    m1 4081.966
## 2    m2 3053.313
## 3    m3 3070.998

ctChisqTest(f2,f1)
## [1] 3.579979e-216
ctChisqTest(f3,f2)
## [1] 0.9999959

l1=ctLOO(f1,cores=1)
l2=ctLOO(f2,cores=1)
l3=ctLOO(f3,cores=1)

print(data.frame(
  Model=c('m1','m2','m3'),
  CV10FoldLogLik=c(l1$outsampleLogLik,l2$outsampleLogLik,l3$outsampleLogLik)))
##   Model CV10FoldLogLik
## 1    m1      -2045.691
## 2    m2      -1525.703
## 3    m3      -1534.040



















