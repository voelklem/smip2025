##########################################################################
### cross-lagged panel model with random intercepts in continuous time ###
##########################################################################

library(ctsem)
library(ctsemOMX)
data(ctExample1)

#Ignoring unobserved heterogeneity:
example1model <- ctModel(n.latent = 2, n.manifest = 2, Tpoints = 6,
                         manifestNames = c('LeisureTime', 'Happiness'),
                         latentNames = c('LeisureTime', 'Happiness'), LAMBDA = diag(2))
example1fit <- ctFit(dat = ctExample1, dataform = "wide", ctmodelobj = example1model)
summary_example1fit <- summary(example1fit, verbose=T)
summary_example1fit
#plot(example1fit)

#Accounting for unobserved heterogeneity:
traitmodel <- ctModel(n.manifest = 2, n.latent = 2, Tpoints = 6,
                         manifestNames = c('LeisureTime', 'Happiness'),
                         latentNames = c('LeisureTime', 'Happiness'), LAMBDA = diag(2),
                         TRAITVAR = "auto")
traitfit   <- ctFit(dat = ctExample1, dataform = "wide", ctmodelobj = traitmodel)
summary_traitfit <- summary(traitfit, verbose=T)
summary_traitfit
#plot(traitfit)

#Comparing the drift matrices:
summary_example1fit$DRIFT
summary_traitfit$DRIFT

summary_example1fit$discreteDRIFT
summary_traitfit$discreteDRIFT





