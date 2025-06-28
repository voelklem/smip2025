#------- CT and cohorts  #
##########################

##################################################################################################
### Generating and fitting data of a cross-sequential design under two different conditions of ###
### cohort equivalence and cohort non-equivalence.                                             ###
##################################################################################################

library(ctsemOMX)
set.seed(818)

### Step 1. Generate data of single group model---------------------------------------------------
complete_model <-ctModel(Tpoints=7, n.latent=1, n.manifest=1, LAMBDA=matrix(c(1),ncol=1),
                 DRIFT=matrix(c(-0.1),nrow=1),MANIFESTVAR=diag(1),DIFFUSION=matrix(c(sqrt(.5)),1,1),
                 CINT=matrix(c(2),ncol=1),T0MEANS=matrix(c(5),ncol=1,nrow=1), T0VAR=matrix(c(sqrt(2)),nrow=1))

N <- 600 #select total N
complete_data  <-ctGenerate(complete_model,n.subjects=N,burnin=0, wide=T)
complete_data  <-data.frame(complete_data)

#plot A----
N_sel <- 100 # select number of individuals to plot
matplot(t(matrix(rep(0:6,N_sel),N_sel,7,byrow=T)), t(complete_data[sample(1:600, N_sel),1:7]), type="l", ylim=c(0,15), ylab="", xlab="")
par(new=T)
plot(0:6, colMeans(complete_data[,1:7]), type="b", lwd=5, ylim=c(0,15), ylab="Dependent variable", xlab="Time")
summary(complete_data)
#----------

### Step 2. Estimate full (single group) model----------------------------------------------------
inits        <- c(-0.1, 5, 4, 0.5, 2, 1.1) #set starting values (not necessary, but speeds up the estimation)
names(inits) <-c("a", "T0means", "T0var", "Q", "cint", "vare")
complete_model_specification <-ctModel(Tpoints=7, n.latent=1, n.manifest=1, LAMBDA=matrix(c(1),ncol=1),
                         DRIFT=matrix(c("a"),nrow=1),
                         MANIFESTVAR=matrix(c("vare"),nrow=1),
                         DIFFUSION=matrix(c("Q"),1,1),CINT=matrix(c("cint"),ncol=1),T0MEANS=matrix(c("T0means"),ncol=1,nrow=1),
                         T0VAR=matrix(c("T0var"),nrow=1),
                         startValues = inits)

complete_model_out <- ctFit(dat=complete_data, ctmodelobj=complete_model_specification, showInits =TRUE, stationary=NULL)
summary(complete_model_out)

### Step 3. Split sample in three cohorts---------------------------------------------------------
mg_data       <- complete_data
mg_data$group <- c(rep("g1",(N/3)), rep("g2",(N/3)), rep("g3",(N/3)))

mg_data[mg_data$group == "g1",][,4:7] <- NA             #delete data group 1
mg_data[mg_data$group == "g2",][,c(1:2,6:7)] <- NA      #delete data group 2
mg_data[mg_data$group == "g3",][,c(1:4)] <- NA          #delete data group 3

### Step 4. Estimate model under the assumption of cohort equivalence-----------------------------
mg_model1_out <- ctFit(dat=mg_data[,-14], ctmodelobj=complete_model_specification, stationary=NULL, showInits =TRUE)
summary(mg_model1_out)

### Step 5. Estimate model under the assumption of cohort equivalence using multiple group sem (equivalent to Step 4)
mg_constrained <- complete_model_specification
mg_constrained$DRIFT        <- "groupfixed" # constrain all parameter to equality across cohorts
mg_constrained$MANIFESTVAR  <- "groupfixed"
mg_constrained$DIFFUSION    <- "groupfixed"
mg_constrained$T0VAR        <- "groupfixed"
mg_constrained$T0MEANS      <- "groupfixed"
mg_constrained$CINT         <- "groupfixed"
mg_constrained_out <-  ctMultigroupFit(dat=mg_data[,-14], groupings=mg_data$group,
                       ctmodelobj=complete_model_specification,
                       fixedmodel=mg_constrained)
summary(mg_constrained_out)


#plot B----
N_sel=100
matplot(t(matrix(rep(0:6,N_sel),N_sel,7,byrow=T)), t(mg_data[sample(1:N, N_sel),1:7]), type="l", ylim=c(0,15), ylab="", xlab="")
par(new=T)
plot(0:6, colMeans(mg_data[,1:7],na.rm = T), type="b", lwd=5, ylim=c(0,15), ylab="Dependent variable", xlab="Time")
#----------

### Step 6. Test hypothesis of differing T0MEANS across cohorts-----------------------------------
mg_T0MEANSfree              <- mg_constrained
mg_T0MEANSfree$T0MEANS      <- "groupfree"
mg_T0MEANSfree_out1          <-  ctMultigroupFit(dat=mg_data[,-14], groupings=mg_data$group,
                                       ctmodelobj=complete_model_specification,
                                       fixedmodel=mg_T0MEANSfree)
summary(mg_T0MEANSfree_out1)
mxCompare(mg_T0MEANSfree_out1$mxobj, mg_constrained_out$mxobj) # Likelihood ratio test (n.s.)

### Step 7 Generate three cohorts that differ in T0MEANS (g2:-1; g3:+2)----------------------------
g1_model <-ctModel(Tpoints=7, n.latent=1, n.manifest=1, LAMBDA=matrix(c(1),ncol=1),
                   DRIFT=matrix(c(-0.1),nrow=1),MANIFESTVAR=diag(1),DIFFUSION=matrix(c(sqrt(0.5)),1,1),
                   CINT=matrix(c(2),ncol=1),T0MEANS=matrix(c(5),ncol=1,nrow=1), T0VAR=matrix(c(sqrt(2)),nrow=1))
g2_model <-ctModel(Tpoints=7, n.latent=1, n.manifest=1, LAMBDA=matrix(c(1),ncol=1),
                   DRIFT=matrix(c(-0.1),nrow=1),MANIFESTVAR=diag(1),DIFFUSION=matrix(c(sqrt(0.5)),1,1),
                   CINT=matrix(c(2),ncol=1),T0MEANS=matrix(c(4),ncol=1,nrow=1), T0VAR=matrix(c(sqrt(2)),nrow=1))
g3_model <-ctModel(Tpoints=7, n.latent=1, n.manifest=1, LAMBDA=matrix(c(1),ncol=1),
                   DRIFT=matrix(c(-0.1),nrow=1),MANIFESTVAR=diag(1),DIFFUSION=matrix(c(sqrt(0.5)),1,1),
                   CINT=matrix(c(2),ncol=1),T0MEANS=matrix(c(7),ncol=1,nrow=1), T0VAR=matrix(c(sqrt(2)),nrow=1))

g1_data  <-ctGenerate(g1_model,n.subjects=(N/3),burnin=0, wide=T)
g2_data  <-ctGenerate(g2_model,n.subjects=(N/3),burnin=0, wide=T)
g3_data  <-ctGenerate(g3_model,n.subjects=(N/3),burnin=0, wide=T)
row.names(g1_data) <- NULL
row.names(g2_data) <- NULL
row.names(g3_data) <- NULL
mg_data_leveldifferences <- data.frame(rbind(g1_data,g2_data,g3_data))
mg_data_leveldifferences$group <- c(rep("g1",(N/3)), rep("g2",(N/3)), rep("g3",(N/3)))
mg_data_leveldifferences[mg_data_leveldifferences$group == "g1",][,4:7] <- NA             #delete data group 1
mg_data_leveldifferences[mg_data_leveldifferences$group == "g2",][,c(1:2,6:7)] <- NA      #delete data group 2
mg_data_leveldifferences[mg_data_leveldifferences$group == "g3",][,c(1:4)] <- NA          #delete data group 3

#plot C----
N_sel=100
matplot(t(matrix(rep(0:6,N_sel),N_sel,7,byrow=T)), t(mg_data_leveldifferences[sample(1:N, N_sel),1:7]), type="l", ylim=c(0,15), ylab="", xlab="")
par(new=T)
plot(0:6, colMeans(mg_data_leveldifferences[,1:7],na.rm = T), type="b", lwd=5, ylim=c(0,15), ylab="Dependent variable", xlab="Time")
#----------

### Step 8. Test hypothesis of differing T0MEANS across cohorts-----------------------------------
mg_T0MEANSfixed            <- mg_constrained
mg_T0MEANSfixed_out        <- ctMultigroupFit(dat=mg_data_leveldifferences[,-14], groupings=mg_data_leveldifferences$group,
                                                ctmodelobj=complete_model_specification,
                                                fixedmodel=mg_T0MEANSfixed)
summary(mg_T0MEANSfixed_out)

mg_T0MEANSfree              <- mg_T0MEANSfixed 
mg_T0MEANSfree$T0MEANS      <- "groupfree"
mg_T0MEANSfree_out2         <-  ctMultigroupFit(dat=mg_data_leveldifferences[,-14], groupings=mg_data_leveldifferences$group,
                                                ctmodelobj=complete_model_specification,
                                                fixedmodel=mg_T0MEANSfree, showInits=TRUE)
summary(mg_T0MEANSfree_out2)
mxCompare(mg_T0MEANSfree_out2$mxobj, mg_T0MEANSfixed_out$mxobj) # Likelihood ratio test (significant)

### Step 9. Test hypothesis of differing T0MEANS and differing DRIFT across cohorts---------------
mg_T0MEANS_DRIFT_free       <- mg_T0MEANSfree
mg_T0MEANS_DRIFT_free$DRIFT <- "groupfree"
mg_T0MEANS_DRIFT_free_out   <-  ctMultigroupFit(dat=mg_data_leveldifferences[,-14], groupings=mg_data_leveldifferences$group,
                                                ctmodelobj=complete_model_specification,
                                                fixedmodel=mg_T0MEANS_DRIFT_free)
summary(mg_T0MEANS_DRIFT_free_out)
mxCompare(mg_T0MEANS_DRIFT_free_out$mxobj, mg_T0MEANSfree_out2$mxobj) # Likelihood ratio test (n.s.)
