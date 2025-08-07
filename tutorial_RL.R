### RL in EMC2 with trend?
rm(list=ls())
#remotes::install_github("ampl-psych/EMC2@RL",dependencies=TRUE, Ncpus=8)
#library(EMC2)
#

wd <- '~/Desktop/dynamic_tutorial/'
# #dat <- EMC2:::loadRData('~/Desktop/dynamic_tutorial/datasets/dataset-trondheim_task-revl.RData')
# dat <- EMC2:::loadRData('~/Desktop/dynamic_tutorial/datasets/data_exp1.RData')
# dat <- dat[!dat$excl,]  # exclude subjects based on criteria of Miletic et al 2021
# dat$subjects <- as.factor(dat$pp)
# dat$s_left <- dat$appear_left
# dat$s_right <- dat$appear_right
# dat$p_left <- ifelse(dat$s_left==dat$low_stim, dat$low_stim_prob, dat$high_stim_prob)
# dat$p_right <- ifelse(dat$s_right==dat$low_stim, dat$low_stim_prob, dat$high_stim_prob)
# dat$R <- factor(dat$choiceDirection, levels=c('left', 'right'))
# dat$S <- factor(ifelse(dat$p_left>dat$p_right, 'left', 'right'), levels=c('left', 'right'))
# dat <- dat[!is.na(dat$R),] # remove NA-responses
# dat <- dat[dat$rt>.15,] # remove too fast responses
# dat$trials <- dat$TrialNumber
# dat$reward <- dat$outcome/100
# dat <- dat[,c('subjects', 'trials', 'S', 'R', 'rt', 'reward', 's_left', 's_right', 'p_left', 'p_right', 'trialBin')]
# head(dat)

## alternatively, revl
# for reversal learning paradigms
Smatch_prereversal <- function(d) {
  # this figures out which responses correspond to the symbol that had the high pay-off probability at the start of the experiment.
  tmp <- d[!duplicated(d[,c('s_left', 'p_left')]),]
  tmp$correct_symbols <- ifelse(tmp$p_left>tmp$p_right, tmp$s_left, tmp$s_right)
  correct_symbols_prereversal <- tmp[!duplicated(tmp$s_left),'correct_symbols']
  return(d$Rs %in% correct_symbols_prereversal)
}

dat <- EMC2:::loadRData('~/Desktop/dynamic_tutorial/datasets/dataset-trondheim_task-revl.RData')
dat$S <- NA
dat[dat$p_left>dat$p_right,'S'] <- 'left'
dat[dat$p_left<dat$p_right,'S'] <- 'right'
dat$S <- factor(dat$S, levels=levels(dat$R))
dat$Smatchprerev <- Smatch_prereversal(dat)
samplers_fn <- file.path(wd, './samples/dat_revl_full.RData')

# Getting the data in the right format ------------------------------------
# For RL paradigms, we need to use a combination of response-coding of accumulators (left/right)
# and *stimulus/symbol-coding* of accumulators (mapping each accumulator onto a symbol that represents each choice option),
RS <- function(d) {
  ## which stimulus *symbol* was chosen?
  ifelse(d$R=='left', d$s_left, d$s_right)
}

lRS <- function(d) {
  ## For race models: In dadms, there's a column 'lR', which codes which latent response (left/right) each accumulator corresponds to.
  ## In RL paradigms we also need to know to which latent response *symbol* the accumulator corresponds. That's lRS
  factor(d[cbind(1:nrow(d), match(paste0('s_', d$lR), colnames(d)))])
}

Smatch <- function(d) {
  ## which stimulus *symbol* was correct?
  ifelse(d$p_left > d$p_right, d$s_left, d$s_right)
}



# this is a function generator that generates one column per symbol, with the appropriate structure in the dadm for updating.
covariate_column_generator <- function(col_name) {
  function(d) {
    d[,col_name] <- NA
    d[d$RS == col_name, col_name] <- d[d$RS == col_name, 'reward']
    d[d$RS != d$lRS, col_name] <- NA
    d[[col_name]]
  }
}
all_symbols <- unique(c(dat$s_left, dat$s_right))
make_covariate_columns <- setNames(lapply(all_symbols, covariate_column_generator), all_symbols)

## Make trend
trend_RL <- make_trend(kernel='delta', base='lin', cov_names=list(all_symbols), par_names='v',
                       premap = FALSE, pretransform = FALSE)

## Make design
design_RDM <- design(model=RDM,
                     data=dat,
                     matchfun=function(d) d$S == d$lR,
                     functions=c(list(RS=RS, Smatch=Smatch, lRS=lRS), make_covariate_columns),
                     formula=list(B ~ 1, v ~ 1, t0 ~ 1),
                     constants=c('v.q0'=0),    # don't try to estimate Q0 -- won't work. Assume it starts at 0 (unbiased)
                     trend=trend_RL)

samplers <- make_emc(dat, design=design_RDM, compress=FALSE)
head(samplers[[1]]$data[[1]])

samplers <- fit(samplers, iter=1000, cores_per_chain=2, cores_for_chains=3, fileName=samplers_fn)
samplers <- EMC2:::loadRData(samplers_fn)
plot_pars(samplers)

# for debugging
undebug(make_data)
debug(EMC2:::run_trend)
pp <- predict(samplers, n_cores=1, conditional_on_data=TRUE, return_covariates=TRUE, return_trialwise_parameters=TRUE)

# dadm_full <- do.call(rbind, samplers[[1]]$data)
# covs<-attr(pp, 'covariates')
# twise_pars<-attr(pp, 'trialwise_parameters')

dat$accuracy <- dat$S==dat$R
pp$accuracy <- pp$S==pp$R

## Aggregations for Exp 1
dat$bin <- as.numeric(cut(dat$trials, breaks=10))
pp$bin <- as.numeric(cut(pp$trials, breaks=10))

# Part 1. Plot fit
aggAccS <- aggregate(accuracy~subjects*bin, dat, mean)
aggAccG <- aggregate(accuracy~bin, aggAccS, mean)

aggRTS <- aggregate(rt~subjects*bin*accuracy, dat,quantile, c(0.1,.5,.9))
aggRTG <- aggregate(rt~bin*accuracy, aggRTS, mean)

# pp
ppaggAccS <- aggregate(accuracy~subjects*bin*postn, pp, mean)
ppaggAccG <- aggregate(accuracy~bin*postn, pp, mean)
ppaggAcc <- aggregate(accuracy~bin, ppaggAccG, quantile, c(0.025, 0.5, 0.975))

ppaggRTS <- aggregate(rt~subjects*bin*accuracy*postn, pp, quantile, c(0.1,.5,.9))
ppaggRTG <- aggregate(rt~bin*accuracy*postn, ppaggRTS, mean)
ppaggRT <- aggregate(cbind(`10%`,`50%`,`90%`)~bin*accuracy, ppaggRTG, quantile, c(0.025, 0.5, 0.975))

## plot: 1. accuracy
par(mfrow=c(1,3))
plot(0,0,type='n', xlim=c(1,10), ylim=c(0.4,.9), ylab='', xlab='Trial bin', main='')#, xaxt=ifelse(condition_=='SPD', 's', 'n'))
abline(h=seq(0,1,.1), col='lightgray', lty=2)
polygon(c(1:10, 10:1), c(ppaggAcc$accuracy[,1],rev(ppaggAcc$accuracy[,3])),col=adjustcolor(2, alpha.f=.3), border = FALSE)
lines(aggAccG$bin, aggAccG$accuracy, lwd=1.5)
points(aggAccG$bin, aggAccG$accuracy, pch=19, lwd=1.5)

# 2. RT (correct)
plot(0,0,type='n', xlim=c(1,10), ylim=range(c(ppaggRT[,3:5],aggRTG[,3:5])), xlab='Trial bin', ylab='RT (s)', main='')#, xaxt=ifelse(condition_=='SPD', 's', 'n'))
abline(h=seq(0,2,.1), col='lightgray', lty=2)
for(quantile_ in c('10%', '50%', '90%')) {
  polygon(c(1:10, 10:1), c(ppaggRT[ppaggRT$accuracy==1,quantile_][,'2.5%'],
                           rev(ppaggRT[ppaggRT$accuracy==1,quantile_][,'97.5%'])),
          col=adjustcolor(2, alpha.f=.3), border = FALSE)

  lines(aggRTG$bin[aggRTG$accuracy==1], aggRTG[aggRTG$accuracy==1, quantile_], lwd=1.5) # data
  points(aggRTG$bin[aggRTG$accuracy==1], aggRTG[aggRTG$accuracy==1, quantile_], pch=19, lwd=1.5) # data
}

plot(0,0,type='n', xlim=c(1,10), ylim=range(c(ppaggRT[,3:5],aggRTG[,3:5])), xlab='Trial bin', ylab='RT (s)', main='')#, xaxt=ifelse(condition_=='SPD', 's', 'n'))
abline(h=seq(0,2,.1), col='lightgray', lty=2)
for(quantile_ in c('10%', '50%', '90%')) {
  polygon(c(1:10, 10:1), c(ppaggRT[ppaggRT$accuracy==0,quantile_][,'2.5%'],
                           rev(ppaggRT[ppaggRT$accuracy==0,quantile_][,'97.5%'])),
          col=adjustcolor(2, alpha.f=.3), border = FALSE)

  lines(aggRTG$bin[aggRTG$accuracy==0], aggRTG[aggRTG$accuracy==0, quantile_], lwd=1.5) # data
  points(aggRTG$bin[aggRTG$accuracy==0], aggRTG[aggRTG$accuracy==0, quantile_], pch=19, lwd=1.5) # data
}


# aggregations for REVL ------------------------------------------------------------
pp$RS <- RS(pp) #pp$Rs==pp$Smatchprerev
pp$acc <- pp$S==pp$R
pp$Racc <- pp$acc==pp$Smatchprerev
#pp$Racc <- pp$RS==pp$Scorrectprerev

aggRT <- aggregate(rt~trialNreversal,dat,quantile, c(.1, .5,.9))
aggRTpp <- aggregate(rt~trialNreversal, aggregate(rt~trialNreversal*postn, pp, quantile, c(.1, .5, .9)), quantile, c(.025, .5,.975))

aggChoice <- aggregate(Racc~trialNreversal,dat,mean)
aggChoicepp <- aggregate(Racc~trialNreversal, aggregate(Racc~trialNreversal*postn,pp,mean), quantile, c(0.025, .5, .975))

par(mfrow=c(2,1))
plot(aggChoice$trialNreversal, aggChoice$Racc, type='b', lwd=2, ylab='Choice = accurate prerev', xlab='Trial N (relative to reversal)')
abline(v=0, lty=2)
polygon(c(aggChoicepp$trialNreversal, rev(aggChoicepp$trialNreversal)),
        c(aggChoicepp$Racc[,1], rev(aggChoicepp$Racc[,3])), col=adjustcolor(2, alpha.f=.4))

plot(aggRT$trialNreversal, aggRT$rt[,2], type='b', lwd=2, ylab='RT (s)', xlab='Trial N (relative to reversal)')
abline(v=0, lty=2)
polygon(c(aggRTpp$trialNreversal, rev(aggRTpp$trialNreversal)),
        c(aggRTpp$`50%`[,1], rev(aggRTpp$`50%`[,3])), col=adjustcolor(2, alpha.f=.4))




# RL-DDM ------------------------------------------------------------------




# Debugging stuff ---------------------------------------------------------


## Debugging
p_vector <- sampled_pars(samplers)
p_vector['B'] <- 1
p_vector['v'] <- 0
p_vector['t0'] <- log(.15)
p_vector['v.B0'] <- 1
p_vector['v.alpha'] <- .1

#undebug(EMC2:::run_trend)
debug(EMC2:::get_pars_matrix)
debug(EMC2:::prep_trend)
EMC2:::get_pars_matrix(p_vector, dadm=samplers[[1]]$data[[1]], model=samplers[[1]]$model())

d_tmp <- d
d_tmp$D <- NA

d1 <- samplers[[1]]$data[[1]]
d2 <- samplers[[1]]$data[[2]]
dc <- rbind(d1, d2)

#debug(EMC2:::run_trend)
trend_pars_c <- t(replicate(nrow(dc), c('v.B0'=1, 'v.q0'=0.5, 'v.alpha'=.1)))
out_R <- EMC2:::run_trend(dadm=dc, trend = trend_RL$v, param=rep(1, nrow(dc)), trend_pars = trend_pars)

trend_pars1 <- t(replicate(nrow(d1), c('v.B0'=1, 'v.q0'=0.5, 'v.alpha'=.1)))
trend_pars2 <- t(replicate(nrow(d2), c('v.B0'=1, 'v.q0'=0.5, 'v.alpha'=.1)))
out_c <- c(EMC2:::run_trend_rcpp(data=d1, trend = trend_RL$v,
                                 param=rep(1, nrow(d1)), trend_pars = trend_pars1),
           EMC2:::run_trend_rcpp(data=d2, trend = trend_RL$v,
                                 param=rep(1, nrow(d2)), trend_pars = trend_pars2))
all(round(out_R[[1]],4) == round(out_c,4))
head(out_R[[1]])
head(out_c)

#out_R <- out_R[[1]]
tail(out_R[[1]][dc$subjects==8],30)
tail(out_c[dc$subjects==8],30)

out_R[[2]][dc$subjects==8,]
# second subject?



## posterior predictives -- generate feedback?
debug(EMC2:::make_data)
debug(EMC2:::make_data_unconditional)

debug(EMC2:::run_trend)

