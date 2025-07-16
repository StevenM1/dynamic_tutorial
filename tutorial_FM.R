rm(list=ls())
# remotes::install_github("ampl-psych/EMC2@pp_conditional",dependencies=TRUE, Ncpus=8)
library(EMC2)
wd <- './'
source(file.path(wd, 'plotting_utils.R'))
wgm <- EMC2:::loadRData(file.path(wd, 'datasets/wagenmakers2004_CS.RData'))
wgm <- EMC2:::add_trials(wgm)
ADmat <- matrix(c(-.5,.5), ncol=1, dimnames=list(NULL,'d'))

# FM ----------------------------------------------------------------------
trend_FM=make_trend(kernel='deltab', base='lin',
                    cov_names ='rt',
                    par_names='B', premap = TRUE, pretransform = FALSE, filter_lR=TRUE)
# trend_FM
ADmat <- matrix(c(-.5,.5), ncol=1, dimnames=list(NULL,'d'))
design <- design(model=RDM,
                 data=wgm,
                 contrast=list(lM=ADmat),
                 matchfun=function(d) d$S==d$lR,
                 formula=list(B ~ 1, v ~ lM, t0 ~ 1, s~lM),
                 transform=list(func=c(B='identity'),          # should also update prior in this case
                                lower=c(B.alpha=0.05)),        
                 constants=c('B.q0'=3, 's'=log(1)),            # don't try to estimate Q0 -- won't work. Assume it starts at 0 (no errors)
                 trend=trend_FM)
samplers <- make_emc(wgm, design=design, compress=FALSE)

samplers <- fit(samplers, iter=1000, cores_per_chain=6, cores_for_chains=3, fileName=file.path(wd, 'samples/wgm_FM.RData'))

# Simulate again
pp <- predict(samplers, n_post=100, n_cores=30, conditional_on_data=FALSE)

plotSpectrum(dat=wgm, pp=pp)


# FM in a dataset with objective difficulty -------------------------------
# The wagenmakers 2004 dataset does not have any objective measure of difficulty, so we cannot inspect
# carry-over effects of difficulty. However, the first block of the second experiment in Miletic & Van Maanen (2019) does.
mvm1 <- EMC2:::loadRData(file.path(wd, 'datasets/mileticvanmaanen2019exp2block1.RData'))
mvm1 <- EMC2:::add_trials(mvm1)
ADmat <- matrix(c(-.5,.5), ncol=1, dimnames=list(NULL,'d'))

design <- design(model=RDM,
                 data=mvm1,
                 contrast=list(lM=ADmat),
                 matchfun=function(d) d$S==d$lR,
                 formula=list(B ~ 1, v ~ lM*difficulty, t0 ~ 1, s~lM),
                 transform=list(func=c(B='identity'),          # should also update prior in this case
                                lower=c(B.alpha=0.05)),        
                 constants=c('B.q0'=3, 's'=log(1),
                             v_difficultyD2=0, v_difficultyD3=0, v_difficultyD4=0,v_difficultyD5=0),            # don't try to estimate Q0 -- won't work. Assume it starts at 0 (no errors)
                 trend=trend_FM)
samplers <- make_emc(mvm1, design=design, compress=FALSE)

# note that we now have 56 subjects instead of 6 (albeit with way fewer trials per subject, this is still a substantial difference).
# It should be a bit slower. This was run on a pretty decent server so we can just increase the compute to compensate in this case.
samplers <- fit(samplers, iter=1000, cores_per_chain=15, cores_for_chains=3, fileName=file.path(wd, 'samples/mvm1_FM.RData'))
# This took ~10 minutes on this server.

samplers <- EMC2:::loadRData(file.path(wd, 'samples/mvm1_FM.RData'))
# Simulate again
pp <- predict(samplers, n_post=100, n_cores=30, conditional_on_data=FALSE)

# In some regimes, dynamical systems can be unstable. E.g., a speed-up
# can lead to an increase in fluency so a decrease in threshold #
# followed by another speed-up, etc etc. In humans, this would lead to a 
# brief set of very fast trials, followed by some cognitive break on that
# In our models, this leads to an Inf response.
sum(is.infinite(pp$rt)) # 4
# This is quite rare though
mean(is.infinite(pp$rt))

# Excluding these trials leads to a similarly-looking spectrum
plotSpectrum(dat=mvm1, pp=pp[!is.infinite(pp$rt),])   # ?

# but what about those carry-over effects?
mvm1$prevD <- Hmisc::Lag(mvm1$difficulty,1)
pp$prevD <- Hmisc::Lag(pp$difficulty,1)
agg_dat <- aggregate(rt~prevD, aggregate(rt~prevD*subjects,mvm1[mvm1$trials>1,],mean), mean)
agg_pp <- aggregate(rt~prevD, aggregate(rt~prevD*postn, aggregate(rt~prevD*postn*subjects,pp[pp$trials>1,],mean), mean), quantile, c(0.025, .975))
plot(1:5, agg_dat$rt, xlab='Difficulty previous trial (1=hard)', ylab='RT (s)', ylim=range(c(agg_pp$rt, agg_dat$rt)))
arrows(1:5, y0=agg_pp$rt[,1], y1=agg_pp$rt[,2], angle=90, code=3)
