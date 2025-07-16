rm(list=ls())
# remotes::install_github("ampl-psych/EMC2@pp_conditional",dependencies=TRUE, Ncpus=8)
library(EMC2)
wd <- './'
source(file.path(wd, 'plotting_utils.R'))
wgm <- EMC2:::loadRData(file.path(wd, 'datasets/wagenmakers2004_CS.RData'))
wgm <- EMC2:::add_trials(wgm)
ADmat <- matrix(c(-.5,.5), ncol=1, dimnames=list(NULL,'d'))


# In an RDM ---------------------------------------------------------------
# In the paper, we used  u = u_0 + w_AM*(1-Q_AM), with Q_AM learning *accuracy*, not *error rate*.
# The (1-Q_AM) parametrisation is currently not possible in EMC2, but we can reverse the logic, and estimate *error rate* directly
isError <- function(d) d$S!=d$R
# and now we can apply u = u_0 + w_AM*(Q_AM), with Q_AM tracking errors.
trend_AM=make_trend(kernel='delta', base='lin', cov_names ='error', par_names='v', premap = TRUE, pretransform = FALSE, filter_lR=TRUE)
trend_AM$v$trend_pnames <- c('w_AM', 'q0_AM', 'alpha_AM')

ADmat <- matrix(c(-.5,.5), ncol=1, dimnames=list(NULL,'d'))
design <- design(model=RDM,
                 data=wgm,
                 contrast=list(lM=ADmat),
                 covariates=c('error'),
                 functions=list(error=isError),              # this *needs* to be passed as a function for posterior predictives to work (correctly)
                 matchfun=function(d) d$S==d$lR,
                 formula=list(B ~ 1, v ~ lM, t0 ~ 1, s~lM),
                 constants=c('q0_AM'=0, 's'=log(1)),         # don't try to estimate Q0 -- won't work. Assume it starts at 0 (no errors)
                 trend=trend_AM)
samplers <- make_emc(wgm, design=design, compress=FALSE)

# samplers <- fit(samplers, iter=1000, cores_per_chain=6, cores_for_chains=3, fileName=file.path(wd, 'samples/wgm_AM.RData'))
samplers <- EMC2:::loadRData(file.path(wd, 'samples/wgm_AM.RData'))

check(samplers)    # learning rate alpha tends to be difficult to estimate
plot_pars(samplers)


# Here, we introduce something new. If you are adapting to a quantity that is an output of the model (e.g., RT, accuracy),
# you need to simulate posterior predictives by looping over trials. 
# After each trial t, you update based on the *simulated* outcome, not based on the outcome observed in the empirical data on trial t.
# That is, the updating should *NOT* be conditional on the empirical data, but conditional on the *simulated* data.
# We can toggle this behavior by setting the argument conditional_on_data to FALSE.
# Note that, since we cannot use a vectorised simulation anymore, it is much slower to simulate new data.
pp <- predict(samplers, n_cores=20, conditional_on_data=FALSE)
pp$accuracy <- 1-pp$error

data_PES <- getErrorEffects(wgm)
pp_PES <- getErrorEffects(pp, mc.cores=10)
plotPES(data_PES=data_PES$average, pp_PES=pp_PES$average, 
        mean_rt=mean(aggregate(rt~subjects,wgm,mean)[,2]),
        main='Error-related effects')


# As an exercise: Compare the above plot with a plot obtained by simulating conditional on the observed data, 
# and to try to understand the difference
pp2 <- predict(samplers, n_cores=20, conditional_on_data=TRUE)
pp2$accuracy <- 1-pp2$error

pp2_PES <- getErrorEffects(pp2, mc.cores=10)
plotPES(data_PES=data_PES$average, pp_PES=pp2_PES$average, 
        mean_rt=mean(aggregate(rt~subjects,wgm,mean)[,2]),
        main='Error-related effects')



# Or in a DDM -------------------------------------------------------------
## Make design
trend_AM_DDM=make_trend(kernel='delta', base='lin', cov_names ='error', par_names='a', premap = TRUE, pretransform = FALSE)
trend_AM_DDM$a$trend_pnames <- c('w_AM', 'q0_AM', 'alpha_AM')

Smat <- matrix(c(-1,1), nrow = 2,dimnames=list(NULL,"dif"))
design_DDM <- design(model=DDM,
                     data=wgm,
                     covariates=c('error'),
                     functions=list(error=isError),   # this *needs* to be passed as a function for posterior predictives to work (correctly)
                     contrasts=list(S=Smat),
                     formula=list(Z ~ 1, v ~ S, a~1, t0 ~ 1),
                     constants=c('q0_AM'=0, 'Z'=qnorm(0.5)),
                     trend=trend_AM_DDM)
sampled_pars(design_DDM)

samplers <- make_emc(wgm, design=design_DDM, compress=FALSE)
# samplers <- fit(samplers, cores_per_chain=6, cores_for_chains=3, fileName=file.path(wd, './samples/wgm_AM_DDM.RData'))
samplers <- EMC2:::loadRData(file.path(wd, 'samples/wgm_AM_DDM.RData'))


pp <- predict(samplers, n_cores=30, conditional_on_data=FALSE)
pp$accuracy <- 1-pp$error

data_PES <- getErrorEffects(wgm)
pp_PES <- getErrorEffects(pp, mc.cores=10)
plotPES(data_PES=data_PES$average, pp_PES=pp_PES$average, 
        mean_rt=mean(aggregate(rt~subjects,wgm,mean)[,2]),
        main='Error-related effects')


## Which model wins?
compare(sList=list(rdm=EMC2:::loadRData(file.path(wd, 'samples/wgm_AM.RData')),
                   ddm=EMC2:::loadRData(file.path(wd, 'samples/wgm_AM_DDM.RData'))), BayesFactor = FALSE)
## 