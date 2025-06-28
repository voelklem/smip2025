#####################
### Input effects ###
#####################

library(ctsem)
library(ctsemOMX)

#Time independent predictors
data(ctExample1TIpred)
tipredmodel <- ctModel(n.manifest = 2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
                       n.TIpred = 1,
                       manifestNames = c('LeisureTime', 'Happiness'),
                       latentNames = c('LeisureTime', 'Happiness'),
                       TIpredNames = 'NumFriends',
                       TRAITVAR = "auto")
tipredfit <- ctFit(dat = ctExample1TIpred, dataform="wide", ctmodelobj = tipredmodel)
summary(tipredfit)['TIPREDEFFECT']

#Time dependent predictors
data(ctExample2)
tdpredmodel <- ctModel(n.manifest = 2, n.latent=2, Tpoints=8, LAMBDA=diag(2), 
                       n.TDpred = 1,
                       manifestNames = c('LeisureTime', 'Happiness'),
                       TDpredNames = 'MoneyInt', latentNames = c('LeisureTime', 'Happiness'),
                       TRAITVAR = "auto")
tdpredfit <- ctFit(dat = ctExample2, dataform="wide", ctmodelobj = tdpredmodel)
summary(tdpredfit)['TDPREDEFFECT']
plot(tdpredfit)
