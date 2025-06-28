#############################################
### Nonlinear dynamics in continuous time ###
#############################################

library(ctsem)

sunspots<-sunspot.year
sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
id <- 1
time <- 1749:1924
datalong <- data.frame(id, time, sunspots)

m <- ctModel(type='stanct',
             manifestNames='sunspots',
             latentNames=c('ss_level', 'ss_velocity'),
             LAMBDA=c( 1, 'ma1|log(1+exp(param))'),
             DRIFT=c(0, 1,
                     '-log1p_exp(freqintercept + freqbylevel * ss_level)','a22'),
             MANIFESTMEANS=c('m1|param * 10 + 44'),
             MANIFESTVAR=diag(0,1), #As per original spec
             CINT=0,
             DIFFUSION=c(0, 0,
                         0, "diffusion"),
             PARS=c('freqintercept', 'freqbylevel'))
#may take a few minutes...
ssfitnl <- ctStanFit(datalong, m) # if you get a warning regarding g++ (not an error), just ignore it:-)
summary(ssfitnl)

ctModelLatex(m, compile=F)


#when sampling
#library(shinystan)
#ssfitnl_shiny <- as.shinystan(ssfitnl$stanfit)
#launch_shinystan(ssfitnl_shiny)