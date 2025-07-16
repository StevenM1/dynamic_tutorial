# SM effects
add_prev <- function(dat, n_hist, levels=c('l', 'r')) {
  for(hist_type in c('S', 'R')) {
    for(trial in 1:n_hist) {
      dat[,paste0(hist_type,'minus',trial)] <- factor(substr(Hmisc::Lag(dat[,hist_type], trial), 1, 1), levels=levels)
    }
  }
  return(dat)
}

plot_history_effects <- function(dat, pp, n_hist=3) {
  dat <- add_prev(dat, n_hist=n_hist, levels=sapply(levels(dat$S), substr, 1, 1))
  pp <- add_prev(pp, n_hist=n_hist, levels=sapply(levels(dat$S), substr, 1, 1))
  
  ## plot
  is_Rone <- levels(dat$R)[1]
  is_Sone <- levels(dat$S)[1]
  
  pp$choice <- pp$R==is_Rone
  dat$choice <- dat$R==is_Rone
  dat$stim <- dat$S==is_Sone
  pp$stim <- pp$S==is_Sone
  
  ## remove
  dat <- dat[!is.na(dat$Sminus3),]
  pp <- pp[!is.na(pp$Sminus3),]
  
  par(mfcol=c(3,2))
  for(hist_type in c('stimuli', 'responses')) {
    colN <- toupper(substr(hist_type,1,1))
    prevD <- prevPP <- NULL
    for(trial in 1:n_hist) {
      prevD <- paste0(prevD, dat[,paste0(colN,'minus',trial)])
      prevPP <- paste0(prevPP, pp[,paste0(colN,'minus',trial)])
    }
    dat$prev <- prevD
    pp$prev <- prevPP
    
    # RTs ~ history. Aggregate: First by subject, then across subjects
    agg1 <- aggregate(rt~subjects*R*prev, dat, mean)
    agg2 <- aggregate(rt~R*prev, agg1, mean)
    ## pp. Aggregate: First by subject, then across subjects, then CI across posteriors
    agg1pp <- aggregate(rt~postn*subjects*R*prev, pp, mean)
    agg2pp <- aggregate(rt~postn*R*prev, agg1pp, mean)
    agg3pp <- aggregate(rt~R*prev, agg2pp, quantile, c(0.025, 0.975))
    
    # 1. RT ~ history
    x <- 1:length(agg2$prev[agg2$R==is_Rone])
    plot(x,agg2$rt[agg2$R==is_Rone], type='b', pch=4, lwd=2, xaxt='n', ylab='RT', xlab=paste0('Previous ', hist_type), ylim=range(c(agg2$rt,agg3pp$rt)))
    points(x,agg2$rt[agg2$R!=is_Rone], type='b', pch=4, col=2, lwd=2)
    axis(side=1, at=x, labels=agg2$prev[agg2$R!=is_Rone])
    legend('topright', c(paste0('choice=', is_Rone), paste0('choice!=', is_Rone)), col=c(1,2), lwd=c(1,1), lty=c(1,1))
    
    arrows(x0=x, y0=agg3pp$rt[agg3pp$R==is_Rone,1], y1=agg3pp$rt[agg3pp$R==is_Rone,2], angle=90, code=3, length=0.025, col=1)
    arrows(x0=x, y0=agg3pp$rt[agg3pp$R!=is_Rone,1], y1=agg3pp$rt[agg3pp$R!=is_Rone,2], angle=90, code=3, length=0.025, col=2)
    
    
    # 2. Choice ~ history
    agg3 <- aggregate(choice~subjects*prev, dat, mean)
    agg4 <- aggregate(choice~prev, agg3, mean)
    
    agg3pp <- aggregate(choice~postn*subjects*prev, pp, mean)
    agg4pp <- aggregate(choice~prev*postn, agg3pp, mean)
    agg5pp <- aggregate(choice~prev, agg4pp, quantile, c(0.025, 0.975))
    
    ylim <- range(c(range(agg4$choice), range(agg5pp$choice)))
    x <- 1:length(agg4$prev)
    plot(x,agg4$choice, type='b', pch=4, lwd=2, xaxt='n', ylab=paste0('Choice == ', is_Rone), xlab=paste0('Previous ', hist_type), ylim=ylim)
    axis(side=1, at=x, labels=agg4$prev)
    abline(h=0.5)
    # PP
    arrows(x0=x, y0=agg5pp[,2][,1], y1=agg5pp[,2][,2], angle=90, code=3, length=0.025, col=1)
    
    # 3. Stimulus ~ history
    agg5 <- aggregate(stim~subjects*prev, dat, mean)
    agg6 <- aggregate(stim~prev, agg5, mean)
    x <- 1:length(agg6$prev)
    plot(x,agg6$stim, type='b', pch=4, lwd=2, xaxt='n', ylab=paste0('Stim == ', is_Rone), xlab=paste0('Previous ', hist_type), ylim=c(0.40, 0.60))
    axis(side=1, at=x, labels=agg6$prev)
    abline(h=0.5)
  }
}

## Spectrum
load_save_spectrum <- function(dat, dat_fn, by.postn=FALSE, detrend, demean, overwrite=FALSE) {
  powerSpectraData <- NULL
  if(!is.null(dat_fn) & !overwrite) {
    fourier_fn <- gsub('.RData', '_fourier.RData', dat_fn)
    if(file.exists(fourier_fn)) {
      powerSpectraData <- EMC2:::loadRData(fourier_fn)
    }
  }
  
  if(is.null(powerSpectraData)) {
    powerSpectraData <- getPowerSpectra(dat, by.postn = by.postn, detrend=detrend, demean=demean)
    if(!by.postn) powerSpectraData <- powerSpectraData[order(powerSpectraData$freq),]  # only order data(?)
    if(!is.null(dat_fn)) {
      fourier_fn <- gsub('.RData', '_fourier.RData', dat_fn)
      save(powerSpectraData, file=fourier_fn)
    }
  }
  return(powerSpectraData)
}

## Fourier plots
getPowerSpectra <- function(data, mean.pp=FALSE, by.postn=FALSE,
                            spans=c(3,5), detrend=TRUE, demean=FALSE, get.log=FALSE) {
  library(parallel)
  if(mean.pp) {
    pp <- data
    pp <- pp[order(pp$subjects, pp$postn, pp$trials),]
    powerBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$spec))
    freqBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$freq))
    
    if(is.list(powerBySubByPostN$rt)) {
      powerBySubByPostN$rt <- do.call(rbind, powerBySubByPostN$rt)
    }
    powerBySub <- lapply(unique(powerBySubByPostN$subjects), function(x) apply(powerBySubByPostN[powerBySubByPostN$subjects==x,'rt'],2,mean))
    meanPower <- apply(do.call(rbind, powerBySub), 2, mean)
    
    freqBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$freq))
    if(is.list(freqBySubByPostN$rt)) {
      freqBySubByPostN$rt <- do.call(rbind, freqBySubByPostN$rt)
    }
    freqBySub <- lapply(unique(freqBySubByPostN$subjects), function(x) apply(freqBySubByPostN[freqBySubByPostN$subjects==x,'rt'],2,mean))
    meanFreqs <- apply(do.call(rbind, freqBySub), 2, mean)
    return(data.frame(freq=meanFreqs, power=meanPower))
  } else if(by.postn) {
    pp <- data   # the 'data' that were passed are actually posterior predictives
    pp <- pp[order(pp$subjects, pp$postn, pp$trials),]  # ensure correct ordering
    
    # get spectra by subject by postN
    spectraByPostN <- mclapply(unique(pp$subjects), function(subject) lapply(unique(pp[pp$subjects==subject,'postn']), function(postn) {
      spectrum(pp[pp$subjects==subject&pp$postn==postn,'rt'], plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)
    }), mc.cores=10)
    
    freqBySubByPostN <- lapply(unique(pp$subjects), function(subject) spectraByPostN[[subject]][[1]]$freq)  # these should all be the same within each subject
    powerBySubByPostN <- lapply(spectraByPostN, function(x) lapply(x, function(y) y$spec))
    npostn <- length(powerBySubByPostN[[1]])
    nsubs <- length(freqBySubByPostN)
    
    # Are are frequencies the same?
    nUniqueLengths <- length(unique(sapply(freqBySubByPostN, length)))
    if(nUniqueLengths > 1) {
      # We need to interpolate because varying subjects have varying trial numbers. Find the subject with fewest trials
      nFrequencies <- min(sapply(freqBySubByPostN, length))
      frequencies <- freqBySubByPostN[[which.min(sapply(freqBySubByPostN, length))]]
      
      # interpolation: ugly for loop, sorry
      interpolatedPowers <- vector(mode='list', length=nsubs)
      for(subject in 1:nsubs) {
        interpolatedPowers[[subject]] <- matrix(NA, nrow=npostn, ncol=nFrequencies)
        for(postn in 1:npostn) {
          f <- approxfun(x=freqBySubByPostN[[subject]], y=powerBySubByPostN[[subject]][[postn]])
          interpolatedPowers[[subject]][postn,] <- f(frequencies)
        }
      }
      
      # now get across-subject mean per posterior predictive
      meanpowerByPostN <- do.call(rbind, lapply(1:npostn, function(postn) apply(do.call(rbind, lapply(interpolatedPowers, function(x) x[postn,])), 2, mean)))
    } else {
      frequencies <- freqBySubByPostN[[1]]
      tmp <- lapply(powerBySubByPostN, function(x) do.call(rbind, x))
      meanpowerByPostN <- do.call(rbind, lapply(1:npostn, function(postn) apply(do.call(rbind, lapply(tmp, function(x) x[postn,])), 2, mean)))
    }
    return(list(freq=frequencies, power=meanpowerByPostN))
  } else {
    # estimate power spectrum by subject
    spectra <- lapply(unique(data$subjects), function(x) spectrum(data[data$subjects==x, 'rt'], plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log))
    freqsBySub <- lapply(spectra, function(x) x$freq)
    powerBySub <- lapply(spectra, function(x) x$spec)
    
    # not all subjects have the same trial numbers.
    # As a result, the actual frequencies at which we have a power estimate differ per subject, so we need to interpolate
    approximators <- mapply(function(x,y) approxfun(x=x,y=y), freqsBySub, powerBySub)
    frequencies <- freqsBySub[[which.min(sapply(freqsBySub, length))]]    # get the frequencies corresponding to the shortest range(?)
    interpolatedPowers <- do.call(rbind, lapply(approximators, function(x) x(frequencies))) # this has shape [nSubjects x nFrequencies]
    
    meanPower <- apply(interpolatedPowers, 2, mean)
    return(data.frame(freq=frequencies, power=meanPower))
  }
}

plotSpectrum <- function(dat, pp, pp2=NULL,
                         dat_fn=NULL, pp_fn=NULL, pp2_fn=NULL,
                         xlab='', ylab='', main='', add.legend=FALSE, ylim=NULL, plot.log=TRUE, xlim=NULL, detrend=FALSE,
                         demean=TRUE,
                         trial_duration=NULL, full_x=FALSE, plot_xticklabels=TRUE, overwrite=TRUE) {
  powerSpectraData <- load_save_spectrum(dat, dat_fn, detrend=detrend, demean=demean, overwrite=overwrite)
  
  if(plot.log) { f <- function(x) log(x) } else { f <- function(x) x }
  if(is.null(trial_duration)) {
    plot(x=f(powerSpectraData$freq), y=f(powerSpectraData$power), xlab=xlab, ylab=ylab, type='n', main=main, ylim=ylim, xlim=xlim)
  } else {
    # frequencies in Hz
    x_axis_ticks_seconds <- c(3600, 30*60, 15*60, 5*60, 2*60, 60, 30, 5, 1)
    x_axis_ticks_hz <- 1/x_axis_ticks_seconds
    x_axis_ticks_1_over_trial <- x_axis_ticks_hz*trial_duration    ## NB
    # ylim <- range(f(powerSpectraData$power))
    
    # which()
    powers <- f(powerSpectraData$power)
    freqs <- f(powerSpectraData$freq)
    
    if(!full_x) {
      idx <- freqs>=(log(x_axis_ticks_1_over_trial)[4]-0.2)
    } else {
      idx <- rep(TRUE, length(powerSpectraData$freq))
    }
    #    xlim <- c(log(x_axis_ticks_1_over_trial)[4], max(freqs))
    #    ylim <- range(powers[freqs>=xlim[1]], na.rm=TRUE)
    
    plot(x=f(powerSpectraData$freq)[idx], y=f(powerSpectraData$power)[idx], xlab=xlab, ylab=ylab, type='n', main=main, ylim=ylim, xlim=xlim, xaxt='n')
    if(plot_xticklabels) {
      axis(side=1, at=log(x_axis_ticks_1_over_trial), labels=c('1 hr', '30 m', '15 m', '5 m', '2 m', '1 m', '30 s', '5 s', '1 s'), las=2)
    } else {
      axis(side=1, at=log(x_axis_ticks_1_over_trial), labels=rep('', length(x_axis_ticks_1_over_trial)), las=2)
    }
  }
  
  
  if(!is.null(pp)) {
    powerSpectraModel <- load_save_spectrum(pp, pp_fn, by.post=TRUE, detrend=detrend, demean=demean, overwrite=overwrite)
    #    powerSpectraModel <- powerSpectraModel(pp, by.postn=TRUE, detrend=detrend, demean=demean)
    for(i in 1:nrow(powerSpectraModel$power)) {
      lines(x=f(powerSpectraModel[[1]]),
            y=f(powerSpectraModel[[2]][i,]),
            col=adjustcolor(3, alpha.f=.1))
    }
    
    if(!is.null(pp2)) {
      powerSpectraModel2 <- load_save_spectrum(pp2, pp2_fn, by.post=TRUE, detrend=detrend, demean=demean)
      #      powerSpectraModel2 <- getPowerSpectra(pp2, by.postn=TRUE, detrend=detrend, demean=demean)
      for(i in 1:nrow(powerSpectraModel2$power)) {
        lines(x=f(powerSpectraModel2[[1]]),
              y=f(powerSpectraModel2[[2]][i,]),
              col=adjustcolor(3, alpha.f=.1))
      }
    }
    
    ## get mean of pp
    freqs <- powerSpectraModel[[1]]
    powers <- apply(powerSpectraModel[[2]],2,mean)
    lines(x=f(freqs), y=f(powers), col='dark red')
    
    if(!is.null(pp2)) {
      freqs2 <- powerSpectraModel2[[1]]
      powers2 <- apply(powerSpectraModel2[[2]],2,mean)
      lines(x=f(freqs2), y=f(powers2), col='dark red')
    }
  }
  lines(x=f(powerSpectraData$freq), y=f(powerSpectraData$power), col=par()$col)
  
  if(!is.null(pp)) lines(x=f(freqs), y=f(powers), col='dark green')
  if(!is.null(pp2)) lines(x=f(freqs2), y=f(powers2), col='dark red')
}

## Post-error effects
plotPES <- function(data_PES=NULL,
                    pp_PES_CI=NULL,
                    mean_rt=NULL,
                    xlab='Trial relative to error',
                    rtmeasure='rtabs',
                    main='Error-related effects', ylim=NULL,
                    ylab='RT (s)') {
  
  ## Plot
  if(is.null(ylim)) ylim <- range(c(pp_PES_CI[pp_PES_CI$measure==rtmeasure,'value'],
                                    data_PES[data_PES$measure==rtmeasure,'value']))
  ylim <- ylim*c(.95, 1.05)
  xs <- unique(data_PES[data_PES$measure==rtmeasure, 'trialNposterror'])
  
  plot(xs, xs, type='n', xlab=xlab, ylab=ylab,  xaxt='n', ylim=ylim, main=main, xlim=range(xs)+c(-.25, .25))
  
  axis(side=1, at=xs)
  if(!is.null(mean_rt)) abline(h=mean_rt, col=par()$col, lty=2)
  abline(v=0, lty=2, col=par()$col)
  
  ## data
  points(data_PES[data_PES$measure==rtmeasure,'trialNposterror'],
         data_PES[data_PES$measure==rtmeasure,'value'], pch=4, lwd=2)
  
  ## posterior predictives
  arrows(x0=pp_PES_CI[pp_PES_CI$measure==rtmeasure,'trialNposterror'],
         y0=pp_PES_CI[pp_PES_CI$measure==rtmeasure,'value'][,1],
         y1=pp_PES_CI[pp_PES_CI$measure==rtmeasure,'value'][,3],
         angle=90, code=3, length=.02, lwd=2, col=2)
  xs <- pp_PES_CI[pp_PES_CI$measure==rtmeasure,'trialNposterror']
  lines(xs, pp_PES_CI[pp_PES_CI$measure==rtmeasure,'value'][,2], col=2)
}

plotPESBySubject <- function(data_PES=NULL,
                             pp_PES_CI=NULL, 
                             orderSubjects=TRUE,
                             xlab='Trial relative to error',
                             rtmeasure='rt',
                             main='', 
                             ylim=NULL,
                             plotAxisTickLabels=TRUE,
                             ylab='RT difference (post-pre error)') {
  
  pp_PES_CI <- pp_PES_CI[pp_PES_CI$probs %in% c('mean', '0.5'),]
  data_PES <- data_PES[data_PES$probs %in% c('mean', '0.5'),]
  
  ## plot all subjects
  idx <- 1:nrow(data_PES)
  if(is.null(ylim)) ylim <- range(c(pp_PES_CI[pp_PES_CI $measure=='rt','value'],
                                    data_PES[data_PES $measure=='rt','value']))
  
  if(orderSubjects) {
    data_PES$subjects <- as.numeric(data_PES$subjects)
    pp_PES_CI$subjects <- as.numeric(pp_PES_CI $subjects)
    tmp <- data_PES[data_PES$probs=='mean'& data_PES$measure=='rt'&data_PES$trialNposterror==1,]
    mapping <- data.frame(plotOrder=1:nrow(tmp), subjects=as.numeric(tmp[order(tmp$value),'subjects']))
    data_PES <- merge(data_PES, mapping)
    pp_PES_CI <- merge(pp_PES_CI, mapping)
    data_PES$subjects <- data_PES$plotOrder
    pp_PES_CI$subjects <- pp_PES_CI$plotOrder
  }
  xs <- as.numeric(unique(data_PES$subjects))
  plot(xs, xs, type='n', xlab=xlab, ylab=ylab,  xaxt='n', ylim=ylim, main=main, xlim=range(xs)+c(-.25, +.25))
  abline(h=0, lty=2, col='grey')
  if(plotAxisTickLabels) {
    axis(side=1, at=xs, labels=unique(data_PES$subjects))
  }
  
  ## data
  points(as.numeric(data_PES[data_PES$measure=='rt'&data_PES$trialNposterror==1,'subjects']),
         data_PES[data_PES$measure=='rt'&data_PES$trialNposterror==1,'value'], pch=4, lwd=2)
  
  ## PP
  arrows(x0=as.numeric(pp_PES_CI[pp_PES_CI$measure=='rt'& pp_PES_CI$trialNposterror==1,'subjects']),
         y0= pp_PES_CI[pp_PES_CI$measure=='rt'&pp_PES_CI$trialNposterror==1,'value'][,1],
         y1= pp_PES_CI[pp_PES_CI$measure=='rt'&pp_PES_CI$trialNposterror==1,'value'][,3],
         angle=90, code=3, length=.02, lwd=2, col=2)
  xs <- as.numeric(pp_PES_CI[pp_PES_CI $measure=='rt'&pp_PES_CI$trialNposterror==1,'subjects'])
  ys <- pp_PES_CI[pp_PES_CI $measure=='rt'& pp_PES_CI$trialNposterror==1,'value'][,2]
  lines(xs[order(xs)], ys[order(xs)], lty=1,col=2)
}


getPostPreDifference <- function(dat, rtcolname='rt', trialRange=c(-3,7)) {
  errors <- dat$accuracy == 0
  ## remove first trial (per participant) and last trial (per participant)
  trialRangeBySub <- aggregate(trials~subjects, dat, range)
  
  dat$exclude_trial <- FALSE
  for(subject in unique(trialRangeBySub$subjects)) {
    ## don't include errors on the first and last trial; can't calculate post-error slowing because of lack of a pre / post error trial
    dat[dat$subjects==subject,'exclude_trial'] = dat[dat$subjects==subject,'trials'] %in% trialRangeBySub[trialRangeBySub$subjects==subject,'trials']
  }
  errors <- errors & !dat$exclude_trial
  errors <- which(errors)
  
  # RT difference post-pre
  rtdf <- expand.grid(trialNposterror=1, probs='mean', measure='rt', value=NA)
  rtdf[rtdf$trialNposterror==1,'value'] <- mean(dat[errors+1,rtcolname] - dat[errors-1,rtcolname])
  
  # Absolute RTs from trials ~-3,7
  rtabsdf <- expand.grid(trialNposterror=seq(trialRange[1],trialRange[2]), probs='mean', measure='rtabs', value=NA)
  for(trialNposterror in seq(trialRange[1],trialRange[2])) {
    trialIdx <- errors+trialNposterror
    trialIdx <- trialIdx[trialIdx>0]  ## remove potential trial numbers below 0
    if(trialNposterror != 0) trialIdx <- trialIdx[!trialIdx %in% errors]  ## only inspect accurate responses, dont include errors in pre/post error trials
    rtabsdf[rtabsdf$trialNposterror==trialNposterror,'value'] <- mean(dat[trialIdx,rtcolname], na.rm=TRUE)
  }
  
  return(rbind(rtabsdf, rtdf))
}

getErrorEffects <- function(dat, mc.cores=1) {
  is_data <- FALSE
  if(!'postn' %in% colnames(dat)) {
    dat$postn <- 1
    is_data <- TRUE
  }
  postpre <-  data.frame(expand.grid(postn=unique(dat$postn), subjects=unique(dat$subjects)))
  diffs = parallel::mclapply(1:nrow(postpre),
                             function(x) {
                               subject <- postpre[x,'subjects']
                               postn <- postpre[x,'postn']
                               out <- getPostPreDifference(dat[dat$subjects==subject&dat$postn==postn,])
                               out$subjects <- subject
                               out$postn <- postn
                               out
                             }, mc.cores=mc.cores)
  diffs <- do.call(rbind, diffs)
  
  if(is_data) {
    # no postn, so only 1 value -- no CI over posterior predictives possible
    avg <- aggregate(value~measure*probs*trialNposterror, diffs, mean)
  } else {
    avg <- aggregate(value~measure*probs*trialNposterror,
                     aggregate(value~measure*probs*trialNposterror*postn, diffs, mean), 
                     quantile, c(0.025, 0.5, 0.975))
  }
  
  return(list(average=avg,diffs=diffs))
}
