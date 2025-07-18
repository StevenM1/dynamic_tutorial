rm(list=ls())
# remotes::install_github("ampl-psych/EMC2@pp_conditional",dependencies=TRUE, Ncpus=8)
library(EMC2)
wd <- './'
source(file.path(wd, 'plotting_utils.R'))
wgm <- EMC2:::loadRData(file.path(wd, 'datasets/wagenmakers2004_CS.RData'))
wgm <- EMC2:::add_trials(wgm)
ADmat <- matrix(c(-.5,.5), ncol=1, dimnames=list(NULL,'d'))
wgm$S01 <- NA
wgm$S01[wgm$S==levels(wgm$S)[1]] <- -1
wgm$S01[wgm$S==levels(wgm$S)[2]] <- 1


# MS3 ---------------------------------------------------------------------
# Now we can simply combine the previous trends in one model, and we have the MS3!
trend_MS3=make_trend(kernel=c('delta', 'delta', 'deltab'), base=c('lin','lin', 'lin'),
                     cov_names =c('S01', 'error', 'rt'),
                     par_names=c('B_lRd', 'v', 'B'),
                     premap = TRUE, pretransform = FALSE, filter_lR=TRUE)

design <- design(model=RDM,
                 data=wgm,
                 contrast=list(lM=ADmat, lR=ADmat),
                 covariates=c('S01'),
                 functions=list(error=function(x) x$S!=x$R), # this *needs* to be passed as a function for posterior predictives to work (correctly)
                 matchfun=function(d) d$S==d$lR,
                 formula=list(B ~ lR, v ~ lM, t0 ~ 1, s~lM),
                 transform=list(func=c(B='identity'),
                                # see paper for explanation  -- enforce variance in Q_FM, to prevent Q_FM becoming a second intercept),
                                lower=c(B.alpha=0.05,
                                        v.alpha=0.01)),
                 constants=c('v.q0'=0, 'B_lRd.q0'=0, 'B.q0'=3,
                             's'=log(1), 'B_lRd'=0),
                 trend=trend_MS3)
samplers <- make_emc(wgm, design=design, compress=FALSE)

fn <- file.path(wd, 'samples/wgm_MS3.RData')
samplers <- fit(samplers, iter=1000, cores_per_chain=6, cores_for_chains=3, fileName=fn)
# This one takes quite a bit longer to converge -- this is a hard model to sample, and 
# the data are not super informative for AM. But still, it convergences with decent
# chains
# took some 40 minutes here

plot(samplers)
plot_pars(samplers)

pp <- predict(samplers, n_cores=30, n_post=100, conditional_on_data=FALSE)

## plot all effects: 1. stimulus history
plot_history_effects(wgm, pp)

## 2. PES. Note that there's a bit of pre-error speeding in the predictions now
pp$accuracy <- 1-pp$error
wgm <- EMC2:::add_trials(wgm)
data_PES <- getErrorEffects(wgm)
pp_PES <- getErrorEffects(pp)
plotPES(data_PES=data_PES$average, pp_PES_CI=pp_PES$average,
        mean_rt=mean(aggregate(rt~subjects,wgm,mean)[,2]),
        main='Error-related effects')

## 3. Fourier spectra
trial_duration = mean(wgm$rt) +.5 # mean RT plus iti etc
plotSpectrum(dat=wgm, pp=pp, trial_duration=trial_duration)


# And add some trends! ----------------------------------------------------
# We can go nuts and add a 4-th order polynomial on the quality of evidence.
# E.g., people might grow tired over time and their efficiency of information sampling could decrease
# Vice versa, it could also increase due to practice. In either case,
# there might be a trend.
trend_MS3t=make_trend(kernel=c('delta', 'delta', 'deltab', 'poly4'), base=c('lin','lin', 'lin', 'add'),
                      cov_names =c('S01', 'error', 'rt', 'trials2'),
                      par_names=c('B_lRd', 'v', 'B', 'v_lMd'),
                      premap = TRUE, pretransform = FALSE, filter_lR=TRUE)

design <- design(model=RDM,
                 data=wgm,
                 contrast=list(lM=ADmat, lR=ADmat),
                 covariates=c('S01'),
                 functions=list(error=function(x) {x$S!=x$R}, # this *needs* to be passed as a function for posterior predictives to work (correctly)
                                trials2=function(x) x$trials/1000),
                 matchfun=function(d) d$S==d$lR,
                 formula=list(B ~ lR, v ~ lM, t0 ~ 1, s~lM),
                 transform=list(func=c(B='identity'),
                                # see paper for explanation  -- enforce variance in Q_FM, to prevent Q_FM becoming a second intercept),
                                lower=c(B.alpha=0.05,
                                        v.alpha=0.01)),
                 constants=c('v.q0'=0, 'B_lRd.q0'=0, 'B.q0'=3,
                             's'=log(1), 'B_lRd'=0),
                 trend=trend_MS3t)
priors <- prior(design, mu_sd=c('v_lMd.d1'=0.1, 'v_lMd.d2'=0.1,
                                'v_lMd.d3'=0.1, 'v_lMd.d4'=0.1))
samplers <- make_emc(wgm, design=design, compress=FALSE, prior_list=priors)

fn <- file.path(wd, 'samples/wgm_MS3T.RData')
#samplers <- fit(samplers, iter=1000, cores_per_chain=6, cores_for_chains=3, fileName=fn)
# This is even harder to sample, of course. Yet, again, give it enough time and it will converge.
# took some 53 minutes here.
samplers <- EMC2:::loadRData(fn)

plot_pars(samplers)

pp <- predict(samplers, n_cores=20, conditional_on_data=FALSE)
plotSpectrum(dat=wgm, pp=pp)


# And the other effects
plot_history_effects(wgm, pp)
pp$accuracy <- 1-pp$error
wgm <- EMC2:::add_trials(wgm)
data_PES <- getErrorEffects(wgm)
pp_PES <- getErrorEffects(pp)
plotPES(data_PES=data_PES$average, pp_PES_CI=pp_PES$average,
        mean_rt=mean(aggregate(rt~subjects,wgm,mean)[,2]),
        main='Error-related effects')
