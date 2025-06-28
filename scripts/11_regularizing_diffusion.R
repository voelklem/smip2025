set.seed(17046)

library(regCtsem)

## Population model:

ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)
# set the DIFFUSION matrix. Note that DIFFUSION eta_1_eta2 is set to 0 in the population.
DIFFUSION=matrix(c(.5,0,0,.5),2)

generatingModel<-ctsem::ctModel(Tpoints=10,n.latent=2,n.TDpred=0,
                                n.TIpred=0,n.manifest=2,
                                MANIFESTVAR=diag(0,2),
                                LAMBDA=diag(1,2),
                                DRIFT=ct_drift,
                                DIFFUSION=DIFFUSION,
                                CINT=matrix(c(0,0),nrow=2),
                                T0MEANS=matrix(0,ncol=1,nrow=2),
                                T0VAR=diag(1,2), type = "omx")

# simulate a training data set
dat <- ctsem::ctGenerate(generatingModel, n.subjects = 100, wide = TRUE)

## Build the analysis model. Note that DIFFUSION eta1_eta2 is freely estimated
# although it is 0 in the population.
myModel <- ctsem::ctModel(Tpoints=10,n.latent=2,n.TDpred=0,
                          n.TIpred=0,n.manifest=2,
                          LAMBDA=diag(1,2),
                          MANIFESTVAR=diag(0,2),
                          CINT=matrix(c(0,0),nrow=2),
                          # Note: we estimate the diffusion:
                          DIFFUSION=matrix(c('eta1_eta1','eta1_eta2',0,'eta2_eta2'),2),
                          T0MEANS=matrix(0,ncol=1,nrow=2),
                          T0VAR="auto", type = "omx")

# fit the model using ctsemOMX:
fit_myModel <- ctsemOMX::ctFit(dat, myModel)

# select diffusion values for regularization:
# Note: If you are unsure what the parameters are called in
# your model, check: showParameters(fit_myModel)
showParameters(fit_myModel)

# regularize the off-diagonal diffusion value:
regIndicators <- c("eta1_eta2")

# Optimize model using GIST with lasso penalty
regModel <- regCtsem(ctsemObject = fit_myModel,
                     dataset = dat,
                     regIndicators = regIndicators,
                     lambdas = "auto",
                     lambdasAutoLength = 20)

# estimates for the off-diagonal diffusion value:
regModel$parameterEstimatesRaw[regIndicators,]
# final estimate:
regModel$parameterEstimatesRaw[regIndicators,which.min(regModel$fit["BIC",])]

# The problem lies in the complex transformation. 
# While we can regularized the lower-triangular diffusion-base values,
# it is not obvious for the user, how these relate to the covariances that we 
# are actually interested in
varianceChol <- matrix(c(1  , 0, 0,
                         -.2, 1, 0,
                         0  ,.3, 1),3,3,byrow = T)
varianceValues <- varianceChol %*% t(varianceChol)
# here, the zero-parameter in the Cholesky matrix results
# in a zero-parameter in the variance:
any(varianceValues == 0)

varianceChol <- matrix(c(1  , 0, 0,
                         -.3, 1, 0,
                         .5 , 0, 1),3,3,byrow = T)
varianceValues <- varianceChol %*% t(varianceChol)
# here, the zero-parameter in the Cholesky matrix does not result
# in a zero-parameter in the variance:
any(varianceValues == 0)

# Now, let's assume we have a sparse covariance matrix:
cov_mat <- matrix(c(1  , -.3,.5,
                    -.3,   1, 0,
                    .5 ,   0, 1),3,3,byrow = T)
# the Cholesky of this covariance matrix is non-sparse:
chol(cov_mat)

# while the Cholesky of the following matrix:
cov_mat <- matrix(c(1  ,-.3,  0,
                    -.3,  1, .5,
                    0  , .5,  1),3,3,byrow = T)
# is sparse:
chol(cov_mat)
