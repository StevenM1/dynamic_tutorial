rm(list=ls())
## Fit SM, AM
#source('./utils.R')
library(EMC2)
dat <- EMC2:::loadRData('./datasets/wagenmakers2004_CS.RData')
head(dat)

# SM ----------------------------------------------------------------------
# Recode into -1 (left) and 1 (right)
dat$S01 <- -1
dat$S01[dat$S=='odd'] <- 1

dat <- EMC2:::add_trials(dat)
dat$trials2 <- dat$trials/max(dat$trials)
dat$trials2 <- dat$trials2-.5


## Make design
#trend_help()
trend_poly=make_trend(kernel='poly2', base='add', cov_names ='trials2', par_names='v_lM', premap = TRUE, pretransform = FALSE)

# Mean-difference parametrisation for v and B
ADmat <- matrix(c(-.5,.5), ncol=1, dimnames=list(NULL,'d'))
design <- design(model=RDM, 
                 data=dat,
                 contrast=list(lM=ADmat, lR=ADmat),
                 covariates=c('S01'),
                 matchfun=function(d) d$S==d$lR,
                 formula=list(B ~ 1, v ~ lM, t0 ~ 1),
                 #transform=list(func=c(B='identity')),    # see note below
                 #constants=c('B_lRd.q0'=0, 'B_lRd'=0),    # don't try to estimate Q0 -- won't work. Assume it starts at 0 (unbiased)
                 trend=trend_poly)
design$model()$transform
sampled_pars(design)

##
samplers <- make_emc(dat, design=design, compress=FALSE)
samplers <- fit(samplers, iter=500, cores_per_chain=3, cores_for_chains=3, fileName='./samples/wagenmakers2004_CS_poly2.RData')


## alternatively, perhaps the threshold changes over time -- people grow more/less cautious?
trend_poly_b=make_trend(kernel='poly2', base='lin', cov_names ='trials2', par_names='b', premap = TRUE, pretransform = FALSE)

# Mean-difference parametrisation for v and B
ADmat <- matrix(c(-.5,.5), ncol=1, dimnames=list(NULL,'d'))
design <- design(model=RDM, 
                 data=dat,
                 contrast=list(lM=ADmat, lR=ADmat),
                 covariates=c('trials2'),
                 matchfun=function(d) d$S==d$lR,
                 formula=list(B ~ 1, v ~ lM, t0 ~ 1),
                 transform=list(func=c(B='identity')),    # see note below
                 #constants=c('B_lRd.q0'=0, 'B_lRd'=0),    # don't try to estimate Q0 -- won't work. Assume it starts at 0 (unbiased)
                 trend=trend_poly_b)
design$model()$transform
sampled_pars(design)

##
samplers <- make_emc(dat, design=design, compress=FALSE)
samplers <- fit(samplers, iter=500, cores_per_chain=3, cores_for_chains=3, fileName='./samples/wagenmakers2004_CS_poly2_b.RData')

