library(ctsemOMX)

# Read in data
dataWide  	<- read.table(file ="4_sim_univariate_delta_1.dat", header = FALSE)
names(dataWide)	<- c(c(paste("Y1_T",0:4, sep="")), c(paste("dT",1:4, sep="")))

# model specification
ct_model4b <- ctModel(type = "omx",
        n.manifest <- 1, 
        n.latent <- 1, 
        Tpoints <- 5, 
        LAMBDA=matrix(c(1),nrow=n.manifest, ncol=n.latent,byrow=TRUE),
        CINT = matrix(c(0),nrow=n.manifest, ncol=n.latent,byrow=TRUE),
        MANIFESTVAR = matrix(c(0),nrow=n.manifest, ncol=n.manifest, byrow=TRUE),
        MANIFESTMEANS = matrix(c(0),nrow=n.manifest, ncol=1, byrow=TRUE))

ct_model4bfit <- ctFit(dataWide, ct_model4b)

summary(ct_model4bfit)
summary(ct_model4bfit)$omxsummary$Minus2LogLikelihood
summary(ct_model4bfit)$DRIFT
summary(ct_model4bfit)$DIFFUSION
summary(ct_model4bfit)$T0VAR
summary(ct_model4bfit)$T0MEANS


