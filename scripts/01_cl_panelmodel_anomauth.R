###################################################
### cross-lagged panel model in continuous time ###
###################################################

#install.packages("ctsem")
library(ctsem)
library(ctsemOMX)

#Bivariate CT model 
data(AnomAuth) 
AnomAuthmodel <- ctModel(	 Tpoints = 5, n.latent = 2, n.manifest = 2,
                           LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
                           MANIFESTVAR=diag(0, 2), 
                           TRAITVAR = NULL) 
AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel)
summary(AnomAuthfit)
summary(AnomAuthfit, verbose=T)
plot(AnomAuthfit)

#Examine the RAM specification
attributes(AnomAuthfit)
AnomAuthfit$mxobj$matrices$A$labels
AnomAuthfit$mxobj$matrices$M$labels
AnomAuthfit$mxobj$matrices$F$values

#Examine the parameter constraints
AnomAuthfit$mxobj$algebras
AnomAuthfit$mxobj$algebras$discreteDRIFT_T1


#1) LR test for testing whether the effect of anomia (eta1) on authoritarianism (eta2) is significantly different from zero
AnomAuthmodel_restricted <- AnomAuthmodel
AnomAuthmodel_restricted$DRIFT[2,1] <- 0
AnomAuthfit_restricted <- ctFit(AnomAuth, AnomAuthmodel_restricted)

summary(AnomAuthfit)$DRIFT
summary(AnomAuthfit_restricted)$DRIFT
chi2 <- summary(AnomAuthfit_restricted)$omxsummary$Minus2LogLikelihood-summary(AnomAuthfit)$omxsummary$Minus2LogLikelihood
pchisq(chi2, 1, lower.tail=F) #p-value
mxCompare(AnomAuthfit$mxobj, AnomAuthfit_restricted$mxobj)

#2) are the two cross effects significantly different from each other?
AnomAuthmodel_restricted2 <- AnomAuthmodel
AnomAuthmodel_restricted2$DRIFT[2,1] <- "cross"
AnomAuthmodel_restricted2$DRIFT[1,2] <- "cross"
AnomAuthfit_restricted2 <- ctFit(AnomAuth, AnomAuthmodel_restricted2)
chi2 <- summary(AnomAuthfit_restricted2)$omxsummary$Minus2LogLikelihood-summary(AnomAuthfit)$omxsummary$Minus2LogLikelihood
pchisq(chi2, 1, lower.tail=F) #p-value
mxCompare(AnomAuthfit$mxobj, AnomAuthfit_restricted2$mxobj)






