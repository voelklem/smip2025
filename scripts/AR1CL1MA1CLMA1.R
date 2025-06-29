###Fitting a CLPM
library(lavaan)
datain <- read.table("../Gallup2011_missingNA_T6.csv", sep=",", header=T)


#AR1CL1MA1CLMA1 model from zyphur (2020)

ar1cl1ma1macl1model <- 
  '
#traits
ksi_SWB =~ lswb1*SWB1 + lswb2*SWB2 + lswb3*SWB3 + lswb4*SWB4 + lswb5*SWB5 + 1*SWB6
ksi_INC =~ linc1*INC1 + linc2*INC2 + linc3*INC3 + linc4*INC4 + linc5*INC5 + 1*INC6

ksi_SWB ~~ ksi_SWB #cov matrix traits
ksi_SWB ~~ ksi_INC 
ksi_INC ~~ ksi_INC 

SWB1d =~ 1*SWB1 # explicitly define latent error terms (to define MA effects later on)
SWB2d =~ 1*SWB2 
SWB3d =~ 1*SWB3
SWB4d =~ 1*SWB4
SWB5d =~ 1*SWB5
SWB6d =~ 1*SWB6

INC1d =~ 1*INC1
INC2d =~ 1*INC2 
INC3d =~ 1*INC3
INC4d =~ 1*INC4
INC5d =~ 1*INC5
INC6d =~ 1*INC6

SWB1 ~ muSWB1*1 #saturated mean structure
SWB2 ~ muSWB2*1
SWB3 ~ muSWB3*1
SWB4 ~ muSWB4*1 
SWB5 ~ muSWB5*1
SWB6 ~ muSWB6*1

INC1 ~ muINC1*1 
INC2 ~ muINC2*1
INC3 ~ muINC3*1
INC4 ~ muINC4*1 
INC5 ~ muINC5*1
INC6 ~ muINC6*1

SWB1d ~~ vars1*SWB1d  #initial cov matrix
INC1d ~~ vari1*INC1d
SWB1d ~~ varsi1*INC1d

SWB2 ~ beta_ss*SWB1 + beta_si*INC1 # ARCL structure
SWB3 ~ beta_ss*SWB2 + beta_si*INC2
SWB4 ~ beta_ss*SWB3 + beta_si*INC3
SWB5 ~ beta_ss*SWB4 + beta_si*INC4
SWB6 ~ beta_ss*SWB5 + beta_si*INC5

INC2 ~ beta_is*SWB1 + beta_ii*INC1
INC3 ~ beta_is*SWB2 + beta_ii*INC2
INC4 ~ beta_is*SWB3 + beta_ii*INC3
INC5 ~ beta_is*SWB4 + beta_ii*INC4
INC6 ~ beta_is*SWB5 + beta_ii*INC5

SWB2 ~ delta_ss*SWB1d + delta_si*INC1d   # MA(1) MACL(1) structure
SWB3 ~ delta_ss*SWB2d + delta_si*INC2d  
SWB4 ~ delta_ss*SWB3d + delta_si*INC3d
SWB5 ~ delta_ss*SWB4d + delta_si*INC4d
SWB6 ~ delta_ss*SWB5d + delta_si*INC5d

INC2 ~ delta_is*SWB1d + delta_ii*INC1d
INC3 ~ delta_is*SWB2d + delta_ii*INC2d
INC4 ~ delta_is*SWB3d + delta_ii*INC3d
INC5 ~ delta_is*SWB4d + delta_ii*INC4d
INC6 ~ delta_is*SWB5d + delta_ii*INC5d

SWB2d ~~ vards2*SWB2d #var/cov of disturbance terms 
SWB3d ~~ vards3*SWB3d
SWB4d ~~ vards4*SWB4d
SWB5d ~~ vards5*SWB5d
SWB6d ~~ vards6*SWB6d

INC2d ~~ vardi2*INC2d
INC3d ~~ vardi3*INC3d
INC4d ~~ vardi4*INC4d
INC5d ~~ vardi5*INC5d
INC6d ~~ vardi6*INC6d

SWB2d ~~ vardsi2*INC2d #contemporaneous effects
SWB3d ~~ vardsi3*INC3d
SWB4d ~~ vardsi4*INC4d
SWB5d ~~ vardsi5*INC5d
SWB6d ~~ vardsi6*INC6d
  '

ar1cl1ma1macl1model_fit <- lavaan(ar1cl1ma1macl1model, data = datain,
              missing = 'ML',                    
              estimator = 'MLR',
              int.ov.free = F,
              int.lv.free = F,
              auto.fix.first = F,
              auto.fix.single = F,
              auto.cov.lv.x = F,
              auto.cov.y = F,
              auto.var = F)
summary(ar1cl1ma1macl1model_fit, standardized = T)
fitMeasures(ar1cl1ma1macl1model_fit)

#parameterestimates(ar1cl1ma1macl1model_fit) #no abbreviations
#library(lavaanPlot) #don't use it. Most confusing plot ever;-)
#lavaanPlot(model = ar1cl1ma1macl1model_fit, node_options = list(shape = "box", fontname = "Helvetica"), edge_options = list(color = "grey"), coefs = TRUE)

