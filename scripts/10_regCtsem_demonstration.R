#if(!require(devtools))install.packages("devtools")
#devtools::install_github("jhorzek/regCtsem") #install regCtsem

library(regCtsem)
set.seed(123)

#### PART 1: Data Simulation ####

# Population parameter values:
DRIFT <- matrix(c(-0.973, 0, 0.434, 
                  0.1, -0.795, 0, 
                  0.264, 0, -2.065),3,3, TRUE)
DIFFUSIONchol <- matrix(c(1.275, 0,0,
                          0.367, 1.177, 0,
                          .806, -.153, 1.414),3,3,TRUE)

generatingModel <- ctModel(LAMBDA = diag(3), n.manifest = 3, n.latent = 3,
                           T0VAR = diag(3), T0MEANS = 0, MANIFESTMEANS = 0,
                           MANIFESTVAR = 0, DRIFT = DRIFT,
                           DIFFUSION = DIFFUSIONchol, TRAITVAR = NULL,
                           Tpoints = 10)

simulatedData <- ctGenerate(ctmodelobj = generatingModel,
                                n.subjects = 100,
                                burnin = 100,
                                wide = T)

#### PART 2: Specify & estimate an unregularized CTSEM ####

DiffusionEstim <- matrix(paste0("Diff",rep(1:3,each = 3), paste0("_Diff", 1:3)),
                         nrow = 3,
                         ncol = 3,
                         byrow = T)
DiffusionEstim[upper.tri(DiffusionEstim)] <- "0"

T0VAREstim <- matrix(paste0("T0VAR", rep(1:3, each = 3), rep(1:3)),3,3, T)
T0VAREstim[upper.tri(T0VAREstim)] <- "0"

analysisModel <- ctModel(LAMBDA = diag(3),
                         n.manifest = 3,
                         n.latent = 3,
                         T0VAR = T0VAREstim,
                         T0MEANS = paste0("T0MEANS", 1:3),
                         MANIFESTMEANS = 0,
                         MANIFESTVAR = 0,
                         DRIFT = "auto",
                         DIFFUSION = DiffusionEstim,
                         TRAITVAR = NULL,
                         Tpoints = 10)

fit.Model <- ctFit(dat = simulatedData, ctmodelobj = analysisModel)

#### PART 3: Specify & estimate a regularized CTSEM ####

# Which parameters do we want to regularize?
regIndicators <- fit.Model$mxobj$DRIFT$labels[!diag(T, 3)] # all cross-effects
print(regIndicators)

# regularization
regModel <- try(regCtsem::regCtsem(ctsemObject = fit.Model,
                                   dataset = simulatedData,
                                   regIndicators = regIndicators,
                                   lambdas = "auto", # the maximally required lambda
                                   # will be computed automatically
                                   lambdasAutoLength  = 20, # note: in practice, we should
                                   # use as many lambdas as possible; we restrict this here to reduce
                                   # the runtime.
                                   penalty = "adaptiveLasso"))

# Plot results and extract best estimates:
plot(regModel)
plot(regModel, what = "fit")
getFinalParameters(regCtsemObject = regModel, criterion = "BIC")

