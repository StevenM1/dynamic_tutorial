wd <- './'
source(file.path(wd, 'plotting_utils.R'))
library(EMC2)

# testing: covariance between alpha and timing estimation? in MVM exp 1
print(load('~/Projects/dynamicEAMs.bak/datasets/timingDM.Rdata'))

colnames(choiceDat) <- c('subjects', 'gender', 'age', 'handedness', 'trials', 'rt', 'response_key', 'accuracy', 'maxtrial')
choiceDat$R <- ifelse(choiceDat$response_key == 'z', 'left', 'right')
choiceDat$S <- NA
choiceDat[choiceDat$R=='left'&choiceDat$accuracy==1,'S'] <- 'left'
choiceDat[choiceDat$R=='left'&choiceDat$accuracy==0,'S'] <- 'right'
choiceDat[choiceDat$R=='right'&choiceDat$accuracy==1,'S'] <- 'right'
choiceDat[choiceDat$R=='right'&choiceDat$accuracy==0,'S'] <- 'left'
choiceDat <- choiceDat[,c('subjects', 'S', 'R', 'rt','trials')]
choiceDat$task <- 'flash'
choiceDat$S <- as.factor(choiceDat$S)
choiceDat$R <- factor(choiceDat$R, levels=levels(choiceDat$S))
choiceDat$task <- as.factor(choiceDat$task)

## timing
colnames(timingDat) <- c('subjects', 'gender', 'age', 'handedness', 'trials', 'rt')
timingDat$S <- timingDat$R <- 'left'
timingDat$S <- factor(timingDat$S, levels=c('left', 'right'))
timingDat$R <- factor(timingDat$R, levels=c('left', 'right'))
timingDat$task <- as.factor('timing')
timingDat <- timingDat[,c('subjects', 'S', 'R', 'rt', 'trials','task')]


# maybe a joint model? -----------------------------------------------------
trend_FM=make_trend(kernel='deltab', base='lin',
                    cov_names ='rt',
                    par_names='B', premap = TRUE, pretransform = FALSE, filter_lR=TRUE)

design_choice <- design(model=RDM,
                 data=choiceDat,
                 contrast=list(lM=ADmat),
                 covariates=c('rt'),
                 matchfun=function(d) d$S==d$lR,
                 formula=list(B ~ 1, v ~ lM, t0 ~ 1),
                 transform=list(func=c(B='identity'),          # should also update prior in this case
                                lower=c(B.alpha=0.05)),
                 constants=c('v_lMd'=0, 'B.q0'=3),
                 trend=trend_FM)

design_timing <- design(model=RDM,
                        data=timingDat[timingDat$trials>50,],
                        matchfun=function(d) d$S==d$lR,
                        formula=list(B ~ 1, v ~ 1, t0 ~ 1, s~1),
                        constants=c(B=1))

samplers <- make_emc(list(choiceDat, timingDat[timingDat$trials>50,]), design=list(design_choice, design_timing), compress=FALSE)
samplers <- fit(samplers, iter=1000, cores_per_chain=8, cores_for_chains=3, fileName=file.path(wd, 'samples_tmp/mvm1_joint2.RData'))
samplers <- EMC2:::loadRData(file.path(wd, 'samples_tmp/mvm1_joint2.RData'))

plot_pars(samplers)

par(mfrow=c(1,1))
plot_relations(samplers, only_cred = T, plot_cred = TRUE)

pars <- parameters(samplers, selection='alpha')
pars$m <- exp(pars$`2|s`)/sqrt(pars$`2|v`)
pars$m
pairs(data.frame(aggregate(.~subjects,pars,mean))[,-1])

#cor(data.frame(aggregate(.~subjects,pars,mean))[,-1])



# separate models below ---------------------------------------------------


# 
# # with FM -----------------------------------------------------------------
# trend_FM=make_trend(kernel='deltab', base='lin',
#                     cov_names ='rt',
#                     par_names='B_taskflash', premap = TRUE, pretransform = FALSE, filter_lR=TRUE)
# 
# design <- design(model=RDM,
#                  data=dat,
#                  contrast=list(lM=ADmat),
#                  covariates=c('rt'),
#                  matchfun=function(d) d$S==d$lR,
#                  formula=list(B ~ task, v ~ lM*task, t0 ~ 1),
#                  transform=list(func=c(B='identity'),          # should also update prior in this case
#                                 lower=c(B_taskflash.alpha=0.05)),
#                  constants=c('v_lMd'=0, 'B_taskflash.q0'=3),
#                  trend=trend_FM)
# samplers <- make_emc(dat, design=design, compress=FALSE)
# 
# samplers <- fit(samplers, iter=1000, cores_per_chain=8, cores_for_chains=3, 
#                 fileName=file.path(wd, 'samples_tmp/mvm1_choice_and_timing_FM.RData'))
# 
# samplers <- EMC2:::loadRData(file.path(wd, 'samples_tmp/mvm1_choice_and_timing_FM.RData'))
# 
# #samplers
# plot_relations(samplers, plot_cred, only_cred = F)
# 
# dat$task <- factor(dat$task, levels=c('flash', 'timing'))
# 
# design <- design(model=RDM,
#                  data=dat,
#                  contrast=list(lM=ADmat),
#                  covariates=c('rt'),
#                  matchfun=function(d) d$S==d$lR,
#                  formula=list(B ~ task, v ~ lM*task, t0 ~ 1, s~task),
#                  transform=list(func=c(B='identity'),          # should also update prior in this case
#                                 lower=c(B_taskflash.alpha=0.05)),
#                  constants=c('v_lMd'=0, 'B_taskflash.q0'=3),
#                  trend=trend_FM)
# samplers <- make_emc(dat, design=design, compress=FALSE)
# 
# samplers <- fit(samplers, iter=1000, cores_per_chain=8, cores_for_chains=3, 
#                 fileName=file.path(wd, 'samples_tmp/mvm1_choice_and_timing_FM2.RData'))
# 
# samplers <- EMC2:::loadRData(file.path(wd, 'samples_tmp/mvm1_choice_and_timing_FM2.RData'))
# 
# 
# 
# design_timing <- design(model=RDM,
#                         data=timingDat, #droplevels(timingDat[as.numeric(as.character(timingDat$subjects))<20,]),
#                         # contrast=list(lM=ADmat),
#                         matchfun=function(d) d$S==d$lR,
#                         formula=list(B ~ 1, v ~ 1, t0 ~ 1))
# 
# samplers <- make_emc(timingDat, design=design, compress=FALSE)
# samplers <- fit(samplers, iter=1000, cores_per_chain=5, cores_for_chains=3, fileName=file.path(wd, 'samples_tmp/mvm1_timing.RData'))
# samplers <- EMC2:::loadRData(file.path(wd, 'samples_tmp/mvm1_timing.RData'))
# 
# pp <- predict(samplers)
# #plot_cdf(timingDM, pp)
# 
# dat <- rbind(timingDat, choiceDat)
