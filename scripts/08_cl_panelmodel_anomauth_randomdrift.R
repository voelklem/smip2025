################################################################
### Hierarchical cross-lagged panel model in continuous time ###
################################################################

#install.packages("ctsem")
library(ctsem)

data(AnomAuth)
#Transforming data from wide to long
AnomAuth_intervals <- ctWideToLong(datawide = AnomAuth, Tpoints=5, n.manifest=2)   #convert wide to long format
AnomAuthlong <- ctDeintervalise(datalong = AnomAuth_intervals, id='id', dT='dT')   #convert intervals to absolute time

#Bivariate CT model -- standard AnomAuth model using stan
AnomAuthmodelstan <- ctModel(	 Tpoints = 5, n.latent = 2, n.manifest = 2,
                           LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
                           MANIFESTVAR=diag(0, 2), 
                           TRAITVAR = NULL, type ="stanct") 
AnomAuthmodelstan$pars[c(1:2,19:20),7] <- FALSE #overwrite default of random T0means and random manifestmeans

AnomAuthfitstan <- ctStanFit(datalong=AnomAuthlong, AnomAuthmodelstan)
summary(AnomAuthfitstan, verbose=T)
ctStanDiscretePars(AnomAuthfitstan, plot=T)


#Allow drift coefficients to vary across individuals
#Disclaimer: This is just to sketch the implementation // takes a while...
AnomAuthmodelstan_rdrift <- AnomAuthmodelstan
AnomAuthmodelstan_rdrift$pars[c(7:10),7] <- T
AnomAuthmodelstan_rdrift_fit <- ctStanFit(datalong=AnomAuthlong, AnomAuthmodelstan_rdrift)
summary(AnomAuthmodelstan_rdrift_fit, verbose=T)
ctModelLatex(AnomAuthmodelstan)
