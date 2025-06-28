#############################
### modeling oscillations ###
#############################

library(ctsem)
library(ctsemOMX)
data(Oscillating)

#interfacing to OpenMx
inits <- c(-39, -.3, 1.01, 10.01, .1, 10.01, 0.05, .9, 0)
names(inits) <- c("crosseffect","autoeffect", "diffusion",
                  "T0var11", "T0var21", "T0var22","m1", "m2", 'manifestmean')

oscillatingm <- ctModel(n.latent = 2, n.manifest = 1, Tpoints = 11,
                        MANIFESTVAR = matrix(c(0), nrow = 1, ncol = 1),
                        LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
                        T0MEANS = matrix(c('m1', 'm2'), nrow = 2, ncol = 1),
                        T0VAR = matrix(c("T0var11", "T0var21", 0, "T0var22"), nrow = 2, ncol = 2),
                        DRIFT = matrix(c(0, "crosseffect", 1, "autoeffect"), nrow = 2, ncol = 2),
                        CINT = matrix(0, ncol = 1, nrow = 2),
                        MANIFESTMEANS = matrix('manifestmean', nrow = 1, ncol = 1),
                        DIFFUSION = matrix(c(0, 0, 0, "diffusion"), nrow = 2, ncol = 2),
                        startValues=inits, type="omx")

oscillating_fit <- ctFit(Oscillating, oscillatingm)
summary(oscillating_fit)
plot(oscillating_fit)


#interfacing to Stan (long-data format)
oscillationlong_intervals <- ctWideToLong(datawide = Oscillating, Tpoints=11, n.manifest=1) #convert wide to long format
oscillationlong <- ctDeintervalise(datalong = oscillationlong_intervals, id='id', dT='dT')   #convert intervals to absolute time
#hist(oscillationlong[,2], ylab = "frequency", xlab="interval length", main="", col="gray")


oscillatingmodel_stan <-ctModel(n.latent = 2, n.manifest=1, Tpoints=11,
                                MANIFESTVAR=matrix(c(0), nrow=1, ncol=1),
                                LAMBDA=matrix(c(1, 0), nrow=1, ncol=2),
                                DRIFT=matrix(c(0, "cross", 1, "auto"), nrow=2, ncol=2),
                                CINT=matrix(c(0,0), ncol=1, nrow=2, ),
                                DIFFUSION=matrix(c(0, 0, 0, "diffusion22"), nrow=2, ncol=2),
                                #startValues = inits,
                                type="stanct")
oscillatingmodel_stan$pars$indvarying <- FALSE # no individual differences
oscillatingfit_stan <-ctStanFit(oscillationlong, oscillatingmodel_stan)
summary(oscillatingfit_stan, verbose=T)
plot(oscillatingfit_stan)







