### RL in EMC2 with trend?
rm(list=ls())
# remotes::install_github("ampl-psych/EMC2@RL",dependencies=TRUE, Ncpus=8);.rs.restartR()
library(EMC2)

wd <- '.'
source(file.path(wd, 'RL_utility_functions.R'))
source(file.path(wd, 'RL_plotting_utils.R'))

# Dataset 1 ---------------------------------------------------------------
dat <- EMC2:::loadRData(file.path(wd, 'datasets/data_exp1.RData'))
dat <- dat[!dat$excl,]  # exclude subjects based on criteria of Miletic et al 2021
dat$subjects <- as.factor(dat$pp)
dat$s_left <- dat$appear_left
dat$s_right <- dat$appear_right
dat$p_left <- ifelse(dat$s_left==dat$low_stim, dat$low_stim_prob, dat$high_stim_prob)
dat$p_right <- ifelse(dat$s_right==dat$low_stim, dat$low_stim_prob, dat$high_stim_prob)
dat$R <- factor(dat$choiceDirection, levels=c('left', 'right'))
dat$S <- factor(ifelse(dat$p_left>dat$p_right, 'left', 'right'), levels=c('left', 'right'))
dat <- dat[!is.na(dat$R),] # remove NA-responses
dat <- dat[dat$rt>.15,] # remove too fast responses
dat$trials <- dat$TrialNumber
dat$reward <- dat$outcome/100
dat <- dat[,c('subjects', 'trials', 'S', 'R', 'rt', 'reward', 's_left', 's_right', 'p_left', 'p_right', 'trialBin')]
head(dat)
samplers_fn <- file.path(wd, './samples/dat_exp1_full.RData')


# function to generate 'rewards'
feedback_generator <- function(d) {
  p_chosen <- ifelse(d$R=='left', d$p_left, d$p_right)
  rbinom(1:nrow(d), 1, p_chosen)
  # alternatively, map symbol to reward probability or reward function here.
  # can change trialwise/subjectwise etc
}

all_symbols <- unique(c(dat$s_left, dat$s_right))
make_covariate_columns <- setNames(lapply(all_symbols, covariate_column_generator), all_symbols)

## Make trend
trend_RL <- make_trend(kernel='delta', base='lin', cov_names=list(all_symbols), par_names='v',
                       premap = FALSE, pretransform = FALSE, filter_lR=FALSE)
trend_RL$v$feedback_generator <- feedback_generator

## Make design
design_RDM <- design(model=RDM,
                     data=dat,
                     matchfun=function(d) d$S == d$lR,
                     functions=c(list(RS=RS, Scorrect=Scorrect, lS=lS), make_covariate_columns),
                     formula=list(B ~ 1, v ~ 1, t0 ~ 1),
                     constants=c('v.q0'=0),    # don't try to estimate Q0 -- won't work. Assume it starts at 0 (unbiased)
                     trend=trend_RL)

samplers <- make_emc(dat, design=design_RDM, compress=FALSE)
head(samplers[[1]]$data[[1]])

samplers <- fit(samplers, iter=1000, cores_per_chain=10, cores_for_chains=3, fileName=samplers_fn)
# Time difference of 8.635182 mins
samplers <- EMC2:::loadRData(samplers_fn)
plot_pars(samplers)

ppC <- predict(samplers, n_cores=30, n_post = 100,
              conditional_on_data=TRUE, return_covariates=TRUE, return_trialwise_parameters=FALSE)
ppU <- predict(samplers, n_cores=30, n_post = 100,
              conditional_on_data=FALSE, return_covariates=TRUE, return_trialwise_parameters=FALSE)

plot_exp1(dat=dat,pp=ppC)
plot_exp1(dat=dat,pp=ppU)


# Dataset 3 = Reversal learning paradigm ---------------------------------------------------------------
dat <- EMC2:::loadRData(file.path(wd, 'datasets/dataset-trondheim_task-revl.RData'))
dat$S <- NA
dat$RS <- RS(dat)  # figures out which symbol was chosen
dat[dat$p_left>dat$p_right,'S'] <- 'left'
dat[dat$p_left<dat$p_right,'S'] <- 'right'
dat$S <- factor(dat$S, levels=levels(dat$R))
dat$chosen_symbol_was_correct_prereversal <- Smatch_prereversal(dat)
samplers_fn <- file.path(wd, './samples/dat_revl_full.RData')

# This is some example code to find all reversals per subject -- it's a bit brute force but does the trick
for(subject in unique(dat$subjects)) {
  tmp <- dat[dat$subjects==subject,]
  all_stimuli <- unique(c(tmp$s_left, tmp$s_right))
  all_stimuli <- setNames(rep(NA, length(all_stimuli)), all_stimuli)
  for(trial in 1:nrow(tmp)) {
    sym_left <- tmp[trial,'s_left']
    sym_right <- tmp[trial,'s_right']
    p_left <- tmp[trial,'p_left']
    p_right <- tmp[trial,'p_right']
    if(is.na(all_stimuli[[sym_left]])) all_stimuli[[sym_left]] <- p_left
    if(is.na(all_stimuli[[sym_right]])) all_stimuli[[sym_right]] <- p_right
    tmp[trial,'hasReversed'] <- all_stimuli[[sym_left]] != p_left
    all_stimuli[[sym_left]] <- p_left
    all_stimuli[[sym_right]] <- p_right
  }
  n_reversal_per_set <- mean(aggregate(hasReversed~s_left, tmp, sum)[,2])*2
  for(stim in names(all_stimuli)) {
    for(reversal_n in 1:n_reversal_per_set) {
      tmp2 <- tmp[tmp$s_left==stim | tmp$s_right==stim,]
      reversal_trial <- which(tmp2$hasReversed)[reversal_n]
      dat[dat$subjects==subject & (dat$s_left==stim | dat$s_right==stim), paste0('trialNrelativetoreversal',reversal_n)] <- 1:nrow(tmp2)-reversal_trial
    }
  }
}
## include additional columns for additional reversals if needed
dat <- dat[,c('subjects', 'S', 'R', 'rt', 'reward', 's_left', 's_right', 'p_left', 'p_right', 'trialNrelativetoreversal1')]

all_symbols <- unique(c(dat$s_left, dat$s_right))
make_covariate_columns <- setNames(lapply(all_symbols, covariate_column_generator), all_symbols)

## Make trend
trend_RL <- make_trend(kernel='delta', base='lin', cov_names=list(all_symbols), par_names='v',
                       premap = FALSE, pretransform = FALSE, filter_lR=FALSE)
trend_RL$v$feedback_generator <- feedback_generator

## Make design
design_RDM <- design(model=RDM,
                     data=dat,
                     matchfun=function(d) d$S == d$lR,
                     functions=c(list(RS=RS, Scorrect=Scorrect, lS=lS), make_covariate_columns),
                     formula=list(B ~ 1, v ~ 1, t0 ~ 1),
                     constants=c('v.q0'=0),    # don't try to estimate Q0 -- won't work. Assume it starts at 0 (unbiased)
                     trend=trend_RL)

samplers <- make_emc(dat, design=design_RDM, compress=FALSE)
head(samplers[[1]]$data[[1]])

samplers <- fit(samplers, iter=1000, cores_per_chain=10, cores_for_chains=3, fileName=samplers_fn)
samplers <- EMC2:::loadRData(samplers_fn)
plot_pars(samplers)

ppC <- predict(samplers, n_cores=30, n_post = 100,
               conditional_on_data=TRUE, return_covariates=TRUE, return_trialwise_parameters=FALSE)
ppU <- predict(samplers, n_cores=30, n_post = 100,
               conditional_on_data=FALSE, return_covariates=TRUE, return_trialwise_parameters=FALSE)

debug(plot_revl)
plot_revl(dat=dat,pp=ppC)

plot_revl(dat=dat,pp=ppU)
plot_revl(dat=dat,pp=ppU,plot_all_RT_quantiles=FALSE)  # focus only on median RTs?



# RL-DDM ------------------------------------------------------------------
## to be implemented later



