library(ctsemOMX)

# Read in data
dataWide  	<- read.table(file ="5_sim_univariate_deltavarying.dat", header = FALSE)
names(dataWide)	<- c(c(paste("Y1_T",0:4, sep="")), c(paste("dT",1:4, sep="")))

# model specification
ct_model5a <- ctModel(type = "omx",
  n.manifest <- 1, 
  n.latent <- 1, 
  Tpoints <- 5, 
  LAMBDA=matrix(c(1),nrow=n.manifest, ncol=n.latent,byrow=TRUE),
  CINT = matrix(c(0),nrow=n.manifest, ncol=n.latent,byrow=TRUE),
  MANIFESTVAR = matrix(c(0),nrow=n.manifest, ncol=n.manifest, byrow=TRUE),
  MANIFESTMEANS = matrix(c(0),nrow=n.manifest, ncol=1, byrow=TRUE))

ct_model5afit <- ctFit(dataWide, ct_model5a)
#summary(ct_model5afit)
#summary(ct_model5afit)$omxsummary$Minus2LogLikelihood
#summary(ct_model5afit)$DRIFT
#summary(ct_model5afit)$DIFFUSION
#summary(ct_model5afit)$T0VAR
#summary(ct_model5afit)$T0MEANS



