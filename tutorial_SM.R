rm(list=ls())
# remotes::install_github("ampl-psych/EMC2@pp_conditional",dependencies=TRUE, Ncpus=8)
library(EMC2)
wd <- getwd()
source(file.path(wd, 'plotting_utils.R'))  # utilities for plotting
wgm <- EMC2:::loadRData(file.path(wd, 'datasets/wagenmakers2004_CS.RData'))
wgm <- EMC2:::add_trials(wgm)  # I just like to have trial numbers already in the data


# RDM + SM ----------------------------------------------------------------------
# Recode into -1 (left) and 1 (right)
wgm$S01 <- NA
wgm$S01[wgm$S==levels(wgm$S)[1]] <- -1
wgm$S01[wgm$S==levels(wgm$S)[2]] <- 1

## Make design
trend_SM=make_trend(kernel='delta', base='lin', cov_names ='S01', par_names='B_lRd', premap = TRUE, pretransform = FALSE, filter_lR=TRUE)
trend_SM$B_lRd$trend_pnames <- c('w_SM', 'q0_SM', 'alpha_SM')
# Note that we recoded 'stimulus' to [-1, 1]; so the mapping can be linear B = B_0 + w*Q_SM.
# Alternatively, you could code S as [0, 1] and the mapping could be centered: B = B_0 + w*(Q_SM - 0.5). In that case, w will be twice as large, as Q_SM is halved.

# Mean-difference parametrisation for v and B
ADmat <- matrix(c(-.5,.5), ncol=1, dimnames=list(NULL,'d'))
design_RDM <- design(model=RDM,
                 data=wgm,
                 contrast=list(lM=ADmat, lR=ADmat),
                 covariates=c('S01'),
                 matchfun=function(d) d$S==d$lR,
                 formula=list(B ~ lR, v ~ lM, t0 ~ 1),
                 transform=list(func=c(B='identity')),    # see note below
                 constants=c('q0_SM'=0, 'B_lRd'=0),    # don't try to estimate Q0 -- won't work. Assume it starts at 0 (unbiased)
                 trend=trend_SM)
sampled_pars(design_RDM)

# Think about why we're *not* doing a exp-transform on threshold here
# This is because we first apply the design matrix to get an accumulator-wise threshold in the dadm
# Only *then* we apply the transformations. So suppose we would sample B on the logscale
# That would mean that we allow B_lR to vary with Q_SM on the log scale. Addition on the logscale is multiplication on the natural scale.
# To prevent allow B_lR to vary with Q_SM on the natural scale, simply estimate B on the natural scale.
# *HOWEVER*: remember your priors! An N(0,1) prior on the natural scale is not a great choice for thresholds.
# For estimation, this isn't much of an issue, but for Bayes Factors, it might be.


# Fit
samplers <- make_emc(wgm, design=design_RDM, compress=FALSE)
samplers <- fit(samplers, iter=1000, cores_per_chain=6, cores_for_chains=3, fileName=file.path(wd, './samples/wgm_SM.RData'))
samplers <- EMC2:::loadRData(file.path(wd, './samples/wgm_SM.RData'))
## You may notice that sampling is not as fast as it was without trends. This is hard stuff to sample!

check(samplers)    # learning rate alpha tends to be difficult to estimate
plot_pars(samplers)

## Visualize estimated Q-values. E.g., take subject 1
all_pars <- parameters(samplers, selection='alpha')
median_pars <- aggregate(.~subjects, all_pars, median)
p_vector <- t(t(median_pars[median_pars$subjects==1,-1]))
updated <- EMC2:::get_pars_matrix(p_vector=p_vector, dadm=samplers[[1]]$data[[1]], model=samplers[[1]]$model())
par(mfrow=c(1,2))
plot(updated[seq(1, nrow(updated),2),'q0_SM'], xlim=c(1,20), type='l')  # Q-values
plot(updated[seq(1, nrow(updated),2),'B']-updated[seq(2, nrow(updated),2),'B'], xlim=c(1,20), type='l')  # threshold bias


## Posterior predictives: What does this look like?
pp <- predict(samplers, n_cores=20)
plot_cdf(wgm, pp, factors = c('S'))  # Overall fit?

## Convenience function to plot the stimulus history effects
plot_history_effects(wgm,pp)
# this does a fairly good job at capturing the stimulus memory effects on choices, and also RTs with some remaining misfit

# Additionally, we can ask predict() to return the updated covariates
pp <- predict(samplers, n_cores=20, return_covariates=TRUE)
covariates <- attr(pp, 'covariates')
# Note that covariates are in the length of the dadm, so in this case, we have *two* rows per trial (one per accumulator).
# To plot the evolution of the covariate
plot(covariates[seq(1, nrow(covariates), 2),'q0_SM'], xlim=c(1, 30))



# DDM + SM -------------------------------------------------------------
## Let's do the same thing, but in a DDM. Here, we want the start point to vary with SM.
trend_SM_DDM=make_trend(kernel='delta', base='lin', cov_names ='S01', par_names='Z', premap = TRUE, pretransform = FALSE)
trend_SM_DDM$Z$trend_pnames <- c('w_SM', 'q0_SM', 'alpha_SM')

Smat <- matrix(c(-1,1), nrow = 2,dimnames=list(NULL,"dif"))
design_DDM <- design(model=DDM,
                 data=wgm,
                 covariates=c('S01'),
                 contrasts=list(S=Smat),
                 formula=list(Z ~ 1, v ~ S, a~1, t0 ~ 1),
                 constants=c('q0_SM'=0, 'Z'=qnorm(0.5)),
                 trend=trend_SM_DDM)
sampled_pars(design_DDM)

samplers <- make_emc(wgm, design=design_DDM, compress=FALSE)
samplers <- fit(samplers, cores_per_chain=6, cores_for_chains=3, fileName=file.path(wd, './samples/wgm_SM_DDM.RData'))
samplers <- EMC2:::loadRData(file.path(wd, './samples/wgm_SM_DDM.RData'))

pp <- predict(samplers, n_cores=20)
plot_cdf(wgm, pp, factors = c('S'))  # Overall fit?

plot_history_effects(wgm,pp)

# What's the better model in this case, RDM or DDM?
compare(sList=list(rdm=EMC2:::loadRData(file.path(wd, 'samples/wgm_SM.RData')),
                   ddm=EMC2:::loadRData(file.path(wd, 'samples/wgm_SM_DDM.RData'))), BayesFactor = FALSE)



# Exercise ----------------------------------------------------------------
# So far, we've assumed that SM biases start points (or threshold differences). An alternative hypothesis could be that
# SM affects drift rates -- e.g., participants selectively attend to the option that was most often present in the
# recent trial history. Implement an RDM or DDM where drift rates vary with SM.

