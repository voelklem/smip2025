# Estimating a proportional latent change score model for T = 5 with measurement error AND dynamic error via OpenMx and Lavaan
require(OpenMx)
require(lavaan)

# 1. read in data and assign names
data_in			<-read.table("3_data_DCSmodel.dat")
names(data_in)	<- c("V1", "V2", "V3", "V4", "V5", "I1", "I2", "I3", "I4")	


# 2. specify and estimate LGC model in OpenMx
LCSmodelMX <- 	mxModel("LCS",
				mxData(data_in, type="raw"),

    A <- mxMatrix("Full", 16, 16, 
           values=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   
					0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					1, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					1, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					1, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
					1, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
					0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
					0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
			free =c(F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, T, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F),
			labels=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,"b",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,"b",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,"b",NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,"b",NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA), byrow=TRUE, name="A"),
				
    S <- mxMatrix("Full", 16, 16,
           values=c(1, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
					0.5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
			free =c(T, T, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					T, T, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, T, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, 
					F, F, F, F, F, F, F, F, F, F, F, F, T, F, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, T, F, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, T, F,
					F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, T),
			labels=	c("var(S)","cov(1S)",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					"cov(1S)","var1",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					NA,NA,NA,"var(z)",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					NA,NA,NA,NA,NA,"var(z)",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					NA,NA,NA,NA,NA,NA,NA,"var(z)",NA,NA,NA,NA,NA,NA,NA,NA,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,"var(z)",NA,NA,NA,NA,NA,NA,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"var(e)",NA,NA,NA,NA, 
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"var(e)",NA,NA,NA,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"var(e)",NA,NA,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"var(e)",NA,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"var(e)"), byrow=TRUE, name="S"),

	Filter <- mxMatrix("Full", 5, 16,
           values=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
			free =	FALSE,
			labels=	NA,
			dimnames=list(c("V1", "V2", "V3", "V4", "V5"),c("S", "I", "eta1", "Deta1", "eta2", "Deta2", "eta3", "Deta3", "eta4", "Deta4", "eta5", "V1", "V2", "V3", "V4", "V5")),			
			byrow=TRUE, name="F"),

	M <- mxMatrix(type="Full", nrow=1, ncol=16,		
			values	=c(2,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
			free	=c(T, T, F, F, F, F, F, F, F, F, F, F, F, F, F, F),
			labels	=c("mS","meta1",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
			dimnames=list(1, c("S", "I", "eta1", "Deta1", "eta2", "Deta2", "eta3", "Deta3", "eta4", "Deta4", "eta5", "V1", "V2", "V3", "V4", "V5")),
			byrow=FALSE, name="M"),
			
	mxExpectationRAM("A","S","F","M"),
	mxFitFunctionML()
)

fitMX	<- mxRun(LCSmodelMX)
summary(fitMX)


# 3. cross-check with lavaan
LCSmodel 	<- "
				#Measurement model
				eta1 =~	1*V1
				eta2 =~	1*V2
				eta3 =~	1*V3
				eta4 =~	1*V4
				eta5 =~	1*V5
								
				#Latent change scores
				Deta1 =~ 1*eta2 			
				Deta2 =~ 1*eta3
				Deta3 =~ 1*eta4	
				Deta4 =~ 1*eta5				

				Deta1 ~ a*eta1
				Deta2 ~ a*eta2
				Deta3 ~ a*eta3	
				Deta4 ~ a*eta4
				
				eta2 ~ 1*eta1 
				eta3 ~ 1*eta2 
				eta4 ~ 1*eta3 
				eta5 ~ 1*eta4 
				
				#Latent slope
				
				S =~ 1*Deta1 + 1*Deta2 + 1*Deta3 + 1*Deta4
				S ~ 1
				
				eta1 ~~ 0*eta1	
				eta2 ~~ 0*eta2 
				eta3 ~~ 0*eta3
				eta4 ~~ 0*eta4
				eta5 ~~ 0*eta5
				
				#Dynamic error
				Deta1 ~~ u*Deta1
				Deta2 ~~ u*Deta2
				Deta3 ~~ u*Deta3
				Deta4 ~~ u*Deta4

				#Measurement error
				V1 ~~ e*V1
				V2 ~~ e*V2
				V3 ~~ e*V3
				V4 ~~ e*V4
				V5 ~~ e*V5
				
				#Initial intercept and variance			
				Int  =~ 1*eta1
				Int ~ 1
				Int ~~ Int			

				#Covariance Int and Slope
				
				Int ~~ S
				
				#Set remaining intercepts to zero
				eta1 ~ 0*1
				eta2 ~ 0*1
				eta3 ~ 0*1
				eta4 ~ 0*1
				eta5 ~ 0*1
				
				Deta1 ~ 0*1
				Deta2 ~ 0*1
				Deta3 ~ 0*1
				Deta4 ~ 0*1				
				"
				
fit 		<- lavaan(LCSmodel, auto.var = TRUE, data=data_in[c("V1", "V2", "V3", "V4", "V5")])
summary(fit, fit.measures = TRUE)	

