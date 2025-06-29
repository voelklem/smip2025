##############################################
### continous time dynamic factor analysis ###
##############################################

library(ctsem)
library(ctsemOMX)
data(ctExample3)

#interfacing to OpenMx (wide-data format)
DFAmodel <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 100,
                    T0MEANS = matrix(0,1,1),
                    LAMBDA = matrix(c(1, 'lambda2', 'lambda3'), nrow = 3, ncol = 1),
                    MANIFESTMEANS = matrix(c('manifestmean1', 'manifestmean2', 'manifestmean3'), nrow = 3, ncol = 1),
                    T0VAR = diag(1),
                    type = "omx")
DFAmodelfit <- ctFit(ctExample3, ctmodelobj = DFAmodel)
summary(DFAmodelfit, verbose=T)
plot(DFAmodelfit)

#Transforming data from wide to long
ctExample3long_intervals <- ctWideToLong(datawide = ctExample3, Tpoints=100, n.manifest=3) #convert wide to long format
ctExample3long <- ctDeintervalise(datalong = ctExample3long_intervals, id='id', dT='dT')   #convert intervals to absolute time

#interfacing to Stan (long-data format)
DFAmodelv2 <- ctModel(n.latent = 1, n.manifest = 3,
                      manifestNames = c('Y1', 'Y2', 'Y3'),
                      latentNames = "eta",
                      id = "id",
                      time = "time",
                      LAMBDA = matrix(c(1, 'lambda2', 'lambda3'), nrow = 3, ncol = 1),
                      MANIFESTMEANS = matrix(c('manifestmean1', 'manifestmean2', 'manifestmean3'), nrow = 3, ncol = 1),
                      T0VAR = diag(1),
                      T0MEANS = matrix(0,1,1),
                      CINT = matrix(0,1,1),
                      type = "ct")

DFAmodelv2fit <- ctStanFit(datalong=ctExample3long, ctstanmodel=DFAmodelv2, indvarying=FALSE, optimize = TRUE, priors = FALSE)
summary(DFAmodelv2fit)
plot(DFAmodelv2fit)

#data visualization
matplot(ctExample3long[,2], ctExample3long[,3:5], type="b",pch=1, col = 1:4, ylab = "outcome", xlab="time", lwd=5)
hist(ctExample3long_intervals[,2], ylab = "frequency", xlab="interval length", main="", col="gray")



