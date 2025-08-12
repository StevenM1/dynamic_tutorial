## RL functions
plot_exp1 <- function(dat, pp) {
  dat$accuracy <- dat$S==dat$R
  pp$accuracy <- pp$S==pp$R

  ## Aggregations for Exp 1
  if(!'trials' %in% colnames(dat)) dat <- EMC2:::add_trials(dat)

  dat$bin <- dat$trialBin #as.numeric(cut(dat$trials, breaks=10))
  pp$bin <- pp$trialBin #as.numeric(cut(pp$trials, breaks=10))

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
}


plot_revl <- function(dat, pp, plot_all_RT_quantiles=TRUE) {
  dat$RS <- RS(dat)
  dat$Racc <- dat$chosen_symbol_was_correct_prereversal <- Smatch_prereversal(dat)
  pp$RS <- RS(pp)
  pp$Racc <- pp$chosen_symbol_was_correct_prereversal <- Smatch_prereversal(pp)
  
  ##
  nReversals <- sum(grepl('trialNrelativetoreversal', colnames(dat)))
  par(mfrow=c(2,nReversals))
  for(i in 1:nReversals) {
    dat$trialNreversal <- dat[,paste0('trialNrelativetoreversal',i)]
    pp$trialNreversal <- pp[,paste0('trialNrelativetoreversal',i)]
  
    aggRT <- aggregate(rt~trialNreversal,dat,quantile, c(.1, .5,.9))
    aggRTpp <- aggregate(rt~trialNreversal, aggregate(rt~trialNreversal*postn, pp, quantile, c(.1, .5, .9)), quantile, c(.025, .5,.975))
  
    aggChoice <- aggregate(Racc~trialNreversal,dat,mean)
    aggChoicepp <- aggregate(Racc~trialNreversal, aggregate(Racc~trialNreversal*postn,pp,mean), quantile, c(0.025, .5, .975))
  
    plot(aggChoice$trialNreversal, aggChoice$Racc, type='b', lwd=2, ylab='Choice = accurate prerev', xlab='Trial N (relative to reversal)')
    polygon(c(aggChoicepp$trialNreversal, rev(aggChoicepp$trialNreversal)),
            c(aggChoicepp$Racc[,1], rev(aggChoicepp$Racc[,3])), col=adjustcolor(2, alpha.f=.4))
    abline(v=0, lty=2)
    
    
    ## RTs: median
    if(plot_all_RT_quantiles) {
      ylim <- range(c(aggRT$rt[,2], quantile(as.matrix(aggRTpp[,-1]), c(0.025, .975))))
    } else {
      ylim <- range(c(aggRT$rt[,2], aggRTpp$`50%`))
    }
    plot(aggRT$trialNreversal, aggRT$rt[,2], type='b', lwd=2, ylab='RT (s)', xlab='Trial N (relative to reversal)', ylim=ylim)
    polygon(c(aggRTpp$trialNreversal, rev(aggRTpp$trialNreversal)),
            c(aggRTpp$`50%`[,1], rev(aggRTpp$`50%`[,3])), col=adjustcolor(2, alpha.f=.4))
    abline(v=0, lty=2)
    
    if(plot_all_RT_quantiles) {
      # and 10th, 90th quantile
      points(aggRT$trialNreversal, aggRT$rt[,1], type='b', lwd=2)
      points(aggRT$trialNreversal, aggRT$rt[,3], type='b', lwd=2)
      polygon(c(aggRTpp$trialNreversal, rev(aggRTpp$trialNreversal)),
              c(aggRTpp$`10%`[,1], rev(aggRTpp$`10%`[,3])), col=adjustcolor(2, alpha.f=.4))
      polygon(c(aggRTpp$trialNreversal, rev(aggRTpp$trialNreversal)),
              c(aggRTpp$`90%`[,1], rev(aggRTpp$`90%`[,3])), col=adjustcolor(2, alpha.f=.4))
    }
  }
}
