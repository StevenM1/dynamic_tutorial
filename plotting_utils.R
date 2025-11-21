# plot_trend <- function(input_data, emc, par_name, subject=1,
#                        lM_filter=NULL, lR_filter=NULL, on_x_axis='trials',
#                        pp_shaded=TRUE,
#                        ...) {
#   
#   
#   if(!is.list(input_data)) {
#     # user supplied p_vector
#     dadm <- emc[[1]]$data[[subject]]
#   } else {
#     dadm <- do.call(rbind, emc[[1]]$data)
#   }
#   filter <- rep(TRUE, nrow(dadm)) & dadm$subjects==subject
#   if(!is.null(lM_filter)) filter <- filter & dadm$lM==lM_filter
#   if(!is.null(lR_filter)) filter <- filter & dadm$lR==lR_filter
#   if(!is.list(input_data)) {
#     updated <- get_pars_matrix(p_vector=input_data,
#                                dadm=dadm,
#                                model=emc[[1]]$model())
#     trend <- updated[filter, par_name]
#     credible_interval <- NULL
#     ylim <- range(trend)
#     x <- dadm[filter, on_x_axis]
#     trend <- trend[order(x)]
#     x <- x[order(x)]
#   } else {
#     # user supplied posterior predictives
#     updated <- lapply(seq_along(input_data), function(i) data.frame(input_data[[i]][filter,par_name,drop=FALSE], x=dadm[filter,on_x_axis]))
#     credible_interval <- aggregate(.~x, do.call(rbind, updated), quantile, c(0.025, .5, .9))
#     
#     ## order
#     credible_interval <- credible_interval[order(credible_interval[[1]]),]
#     trend <- credible_interval[[2]][,2]  # median
#     x <- credible_interval[[1]]
#     if(pp_shaded) {
#       ylim <- range(credible_interval[[2]])
#     } else {
#       ylim <- range(sapply(updated, function(x) range(x[,par_name])))
#     }
#   }
#   
#   # parse dots for plotting arguments
#   dots <- list(...)
#   if(!'xlab' %in% names(dots)) dots$xlab=on_x_axis
#   if(!'ylim' %in% names(dots)) dots$ylim=ylim
#   
#   full_args <- c(list(x=x,
#                       y=trend,
#                       type='l', ylab=par_name), dots)
#   do.call(plot, full_args)
#   
#   if(!is.null(credible_interval)) {
#     if(pp_shaded) {
#       polygon(x=c(x, rev(x)),
#               y=c(credible_interval[[2]][,1], rev(credible_interval[[2]][,3])), col=adjustcolor(3, alpha.f = .3), border=adjustcolor(1, alpha.f=.4))
#     } else {
#       for(i in updated) lines(i$x, i[,par_name], lwd=.2)
#     }
#   }
# }


# SM effects
add_prev <- function(dat, n_hist, levels=c('l', 'r')) {
  for(hist_type in c('S', 'R')) {
    for(trial in 0:n_hist) {
      dat[,paste0(hist_type,'minus',trial)] <- factor(substr(Hmisc::Lag(dat[,hist_type], trial), 1, 1), levels=levels)
    }
  }
  return(dat)
}

legend_outside <- function(position = c("topright", "topleft"),
                           legend_text, fill = NULL, col = NULL, lty = NULL, lwd = NULL, pch = NULL,
                           border = NA, bty = "n",
                           inset_frac = 0.005,
                           ...) {
  position <- match.arg(position)
  
  usr <- par("usr")
  x_range <- usr[2] - usr[1]
  y_range <- usr[4] - usr[3]
  
  # Build args list for legend measurement
  legend_args <- list(x = usr[1], y = usr[4], legend = legend_text,
                      border = border, bty = bty, plot = FALSE, ...)
  
  if (!is.null(fill)) legend_args$fill <- fill
  if (!is.null(col)) legend_args$col <- col
  if (!is.null(lty)) legend_args$lty <- lty
  if (!is.null(lwd)) legend_args$lwd <- lwd
  if (!is.null(pch)) legend_args$pch <- pch
  
  # Measure legend size without plotting
  leg_info <- do.call(graphics::legend, legend_args)
  
  legend_height <- leg_info$rect$h
  legend_width <- leg_info$rect$w
  
  offset_y <- legend_height + inset_frac * y_range
  inset_x <- inset_frac * x_range
  
  par(xpd = TRUE)
  
  # Prepare args for actual legend drawing
  legend_args$x <- NULL
  legend_args$y <- NULL
  legend_args$plot <- NULL
  
  if (position == "topleft") {
    legend_args$x <- usr[1] + inset_x
    legend_args$y <- usr[4] + offset_y
    # legend_args$adj <- c(0, 0.5)
  } else if (position == "topright") {
    legend_args$x <- usr[2] - legend_width - inset_x
    legend_args$y <- usr[4] + offset_y
  # legend_args$adj <- c(0, 0.5)  # left align, vertical center
  }
  
  do.call(graphics::legend, legend_args)
  
  par(xpd = FALSE)
}

rev_str <- function(x) paste(rev(strsplit(x, "")[[1]]), collapse = "")

recode_RA <- function(dat) {
  for(hist_type in c('S', 'R')) {
    cname <- paste0(hist_type, 'rep')
    # Recode S and R into repetition/alternation
    tmp1 <- ifelse(as.numeric(dat[,hist_type])==1,1,-1)
    tmp2 <- Hmisc::Lag(tmp1, 1)
    dat[,cname] <- tmp1*tmp2
    dat[dat$trials==1,cname] <- 0
    dat[is.na(dat[,cname]),cname] <- 0
  }
  
  dat[,'S'] <- factor(dat$Srep, levels=c(1,-1), labels=c('R', 'A'))
  dat[,'R'] <- factor(dat$Rrep, levels=c(1,-1), labels=c('R', 'A'))
  dat[,!colnames(dat) %in% c('Srep', 'Rrep')]
}

plot_history_effects <- function(dat, pp=NULL, hist_type='S', n_hist=3, hist_from=1, p_random=0.5, degree=1, ...) {
  opts <- list(...)
  
  if(length(degree) > 1) {
    par(mfrow=c(2,2))
    for(degree_ in degree) {
      plot_history_effects(dat=dat, pp=pp, hist_type=hist_type, n_hist=n_hist, p_random=p_random, degree=degree_, nopar=TRUE, ...)
    }
    return(invisible(NULL))
  } else {
    if(!'nopar' %in% names(opts)) par(mfcol=c(1,2))
    #plot_history_effects(dat=dat, pp=pp, hist_type=hist_type, n_hist=n_hist, p_random=p_random, degree=degree, ...)
  }
  
  if(degree == 2) {
    dat <- recode_RA(dat)
    dat <- dat[!is.na(dat$S),]
    
    if(!is.null(pp)) {
      pp <- recode_RA(pp)
      pp <- pp[!is.na(pp$S),]
    }
  }
  dat <- add_prev(dat, n_hist=n_hist, levels=sapply(levels(dat$S), substr, 1, 1))
  if(!is.null(pp)) {
    pp <- add_prev(pp, n_hist=n_hist, levels=sapply(levels(dat$S), substr, 1, 1))
  }
  ## plot
  is_Rone_label <- is_Rone <- levels(dat$R)[1]
  is_Rtwo_label <- is_Rtwo <- levels(dat$R)[2]
  is_Sone_label <- is_Sone <- levels(dat$S)[1]
  
  if(!is.null(pp)) {
    pp$choice <- pp$R==is_Rone
  }
  dat$choice <- dat$R==is_Rone
  dat$stim <- dat$S==is_Sone
  if(!is.null(pp)) {
    pp$stim <- pp$S==is_Sone
  }
  
  if(degree == 2) {
    is_Rone_label <- ifelse(is_Rone=='R', 'Rep.', 'Alt.') 
    is_Rtwo_label <- ifelse(is_Rtwo=='R', 'Rep.', 'Alt.') 
    is_Sone_label <- ifelse(is_Sone=='R', 'Rep.', 'Alt.') 
  }
  
  ## remove
  dat <- dat[!is.na(dat[,paste0('Sminus', n_hist)]),]
  if(!is.null(pp)) {
    pp <- pp[!is.na(pp[,paste0('Sminus', n_hist)]),]
  }
  
  colN <- toupper(substr(hist_type,1,1))
  prevD <- prevPP <- NULL
  for(trial in hist_from:n_hist) {
    prevD <- paste0(prevD, dat[,paste0(colN,'minus',trial)])
    if(!is.null(pp)) {
      prevPP <- paste0(prevPP, pp[,paste0(colN,'minus',trial)])
    }
  }
  dat$prev <- prevD
  if(!is.null(pp)) {
    pp$prev <- prevPP
  }
  
  # RTs ~ history. Aggregate: First by subject, then across subjects
  agg1 <- aggregate(rt~subjects*R*prev, dat, mean)
  agg2 <- aggregate(rt~R*prev, agg1, mean)

  ylim1 <- range(c(agg2$rt))
  if(!is.null(pp)) {
    
    ## pp. Aggregate: First by subject, then across subjects, then CI across posteriors
    agg1pp <- aggregate(rt~postn*subjects*R*prev, pp, mean)
    agg2pp <- aggregate(rt~postn*R*prev, agg1pp, mean)
    agg3pp <- aggregate(rt~R*prev, agg2pp, quantile, c(0.025, .5, 0.975))
    #agg3pp$rt[is.infinite(agg3pp$rt)] <- max(agg3pp$rt[!is.infinite(agg3pp$rt)])
    ylim1 <- range(c(agg2$rt,agg3pp$rt))
  }
  # 1. RT ~ history
  x <- 1:length(agg2$prev[agg2$R==is_Rone])
  plot(x,agg2$rt[agg2$R==is_Rone], pch=4, lwd=1, xaxt='n', ylab='RT', xlab=paste0('Previous ', hist_type), ylim=ylim1)
  points(x,agg2$rt[agg2$R!=is_Rone], pch=4, col=2, lwd=1)
  
  labels_rev <- sapply(agg2$prev[agg2$R!=is_Rone], rev_str)
  axis(side=1, at=x, labels=labels_rev, las=3)
  
  if(!is.null(pp)) {
    legend_outside('topleft', legend_text=c('data', 'model'), col=c(1,1), lwd=c(NA,1), lty=c(NA,1), pch=c(4,NA), bty = "n")
  } else {
    legend_outside('topleft', legend_text=c('data'), col=c(1), lwd=c(1), lty=c(1), pch=c(4), bty = "n")
  }
  legend_outside('topright', legend_text=c(paste0('R=', is_Rone_label), paste0('R=', is_Rtwo_label)), fill=c(1,2), border=NA, bty='n')
  
  if(!is.null(pp)) {
    arrows(x0=x, y0=agg3pp$rt[agg3pp$R==is_Rone,1], y1=agg3pp$rt[agg3pp$R==is_Rone,3], angle=90, code=3, length=0.025, col=1)
    arrows(x0=x, y0=agg3pp$rt[agg3pp$R!=is_Rone,1], y1=agg3pp$rt[agg3pp$R!=is_Rone,3], angle=90, code=3, length=0.025, col=2)
    lines(x, agg3pp$rt[agg3pp$R==is_Rone,2], col=1)
    lines(x, agg3pp$rt[agg3pp$R!=is_Rone,2], col=2)
  } else {
    lines(x,agg2$rt[agg2$R==is_Rone])
    lines(x,agg2$rt[agg2$R!=is_Rone], col=2)
  }
  # 2. Choice ~ history
  agg3 <- aggregate(choice~subjects*prev, dat, mean)
  agg4 <- aggregate(choice~prev, agg3, mean)
  
  if(!is.null(pp)) {
    agg3pp <- aggregate(choice~postn*subjects*prev, pp, mean)
    agg4pp <- aggregate(choice~prev*postn, agg3pp, mean)
    agg5pp <- aggregate(choice~prev, agg4pp, quantile, c(0.025,.5, 0.975))
  }
  # stimulus~history
  agg5 <- aggregate(stim~subjects*prev, dat, mean)
  agg6 <- aggregate(stim~prev, agg5, mean)
  
  if(!is.null(pp)) {
    ylim2 <- range(c(range(agg4$choice), range(agg5pp$choice), range(agg6$stim)))
  } else {
    ylim2 <- range(c(range(agg4$choice), range(agg6$stim)))
  }
  x <- 1:length(agg4$prev)
  plot(x,agg4$choice, type='n', pch=4, lwd=1, xaxt='n', ylab=paste0('P (', is_Rone_label, ')'), xlab=paste0('Previous ', hist_type), ylim=ylim2)
  abline(h=p_random, lty=2)
  # PP
  if(!is.null(pp)) {
    arrows(x0=x, y0=agg5pp[,2][,1], y1=agg5pp[,2][,3], angle=90, code=3, length=0.025, col=1)
    lines(x=x, y=agg5pp[,2][,2], col=1)
  } else {
    lines(x,agg4$choice, col=1)
  }
  # 3. Stimulus ~ history
  x <- 1:length(agg6$prev)
  points(x, agg6$stim, pch=20, lwd=1, col=2) # 
  points(x, agg4$choice, pch=4, lwd=1)  # plot again to overwrite filled circle
  
  labels_rev <- sapply(agg6$prev, rev_str)
  axis(side=1, at=x, labels=labels_rev, las=3)
  if(!is.null(pp)) {
    legend_outside('topleft', legend_text=c('data', 'model'), col=c(1,1), lwd=c(NA,1), lty=c(NA,1), pch=c(4,NA), bty = "n")
  } else {
    legend_outside('topleft', legend_text=c('data'), col=c(1), lwd=c(1), lty=c(1), pch=c(4), bty = "n")
  }
  legend_outside('topright', legend_text=c('stimulus'), pch=c(20), col=2, border=NA, bty='n')
}

# plot_history_effects(dat, pp)

# ## Spectrum
# load_save_spectrum <- function(dat, dat_fn, by.postn=FALSE, detrend, demean, overwrite=FALSE) {
#   powerSpectraData <- NULL
#   if(!is.null(dat_fn) & !overwrite) {
#     fourier_fn <- gsub('.RData', '_fourier.RData', dat_fn)
#     if(file.exists(fourier_fn)) {
#       powerSpectraData <- EMC2:::loadRData(fourier_fn)
#     }
#   }
#   
#   if(is.null(powerSpectraData)) {
#     powerSpectraData <- getPowerSpectra(dat, by.postn = by.postn, detrend=detrend, demean=demean)
#     if(!by.postn) powerSpectraData <- powerSpectraData[order(powerSpectraData$freq),]  # only order data(?)
#     if(!is.null(dat_fn)) {
#       fourier_fn <- gsub('.RData', '_fourier.RData', dat_fn)
#       save(powerSpectraData, file=fourier_fn)
#     }
#   }
#   return(powerSpectraData)
# }
# 
# ## Fourier plots
# getPowerSpectra <- function(data, mean.pp=FALSE, by.postn=FALSE,
#                             spans=c(3,5), detrend=TRUE, demean=FALSE, get.log=FALSE) {
#   library(parallel)
#   if(mean.pp) {
#     pp <- data
#     pp <- pp[order(pp$subjects, pp$postn, pp$trials),]
#     powerBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$spec))
#     freqBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$freq))
#     
#     if(is.list(powerBySubByPostN$rt)) {
#       powerBySubByPostN$rt <- do.call(rbind, powerBySubByPostN$rt)
#     }
#     powerBySub <- lapply(unique(powerBySubByPostN$subjects), function(x) apply(powerBySubByPostN[powerBySubByPostN$subjects==x,'rt'],2,mean))
#     meanPower <- apply(do.call(rbind, powerBySub), 2, mean)
#     
#     freqBySubByPostN <- aggregate(rt~postn*subjects, pp, function(x) as.numeric(spectrum(x, plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)$freq))
#     if(is.list(freqBySubByPostN$rt)) {
#       freqBySubByPostN$rt <- do.call(rbind, freqBySubByPostN$rt)
#     }
#     freqBySub <- lapply(unique(freqBySubByPostN$subjects), function(x) apply(freqBySubByPostN[freqBySubByPostN$subjects==x,'rt'],2,mean))
#     meanFreqs <- apply(do.call(rbind, freqBySub), 2, mean)
#     return(data.frame(freq=meanFreqs, power=meanPower))
#   } else if(by.postn) {
#     pp <- data   # the 'data' that were passed are actually posterior predictives
#     pp <- pp[order(pp$subjects, pp$postn, pp$trials),]  # ensure correct ordering
#     
#     # get spectra by subject by postN
#     spectraByPostN <- mclapply(unique(pp$subjects), function(subject) lapply(unique(pp[pp$subjects==subject,'postn']), function(postn) {
#       spectrum(pp[pp$subjects==subject&pp$postn==postn,'rt'], plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log)
#     }), mc.cores=10)
#     
#     freqBySubByPostN <- lapply(unique(pp$subjects), function(subject) spectraByPostN[[subject]][[1]]$freq)  # these should all be the same within each subject
#     powerBySubByPostN <- lapply(spectraByPostN, function(x) lapply(x, function(y) y$spec))
#     npostn <- length(powerBySubByPostN[[1]])
#     nsubs <- length(freqBySubByPostN)
#     
#     # Are are frequencies the same?
#     nUniqueLengths <- length(unique(sapply(freqBySubByPostN, length)))
#     if(nUniqueLengths > 1) {
#       # We need to interpolate because varying subjects have varying trial numbers. Find the subject with fewest trials
#       nFrequencies <- min(sapply(freqBySubByPostN, length))
#       frequencies <- freqBySubByPostN[[which.min(sapply(freqBySubByPostN, length))]]
#       
#       # interpolation: ugly for loop, sorry
#       interpolatedPowers <- vector(mode='list', length=nsubs)
#       for(subject in 1:nsubs) {
#         interpolatedPowers[[subject]] <- matrix(NA, nrow=npostn, ncol=nFrequencies)
#         for(postn in 1:npostn) {
#           f <- approxfun(x=freqBySubByPostN[[subject]], y=powerBySubByPostN[[subject]][[postn]])
#           interpolatedPowers[[subject]][postn,] <- f(frequencies)
#         }
#       }
#       
#       # now get across-subject mean per posterior predictive
#       meanpowerByPostN <- do.call(rbind, lapply(1:npostn, function(postn) apply(do.call(rbind, lapply(interpolatedPowers, function(x) x[postn,])), 2, mean)))
#     } else {
#       frequencies <- freqBySubByPostN[[1]]
#       tmp <- lapply(powerBySubByPostN, function(x) do.call(rbind, x))
#       meanpowerByPostN <- do.call(rbind, lapply(1:npostn, function(postn) apply(do.call(rbind, lapply(tmp, function(x) x[postn,])), 2, mean)))
#     }
#     return(list(freq=frequencies, power=meanpowerByPostN))
#   } else {
#     # estimate power spectrum by subject
#     spectra <- lapply(unique(data$subjects), function(x) spectrum(data[data$subjects==x, 'rt'], plot=FALSE, spans=spans, detrend=detrend, demean=demean, log=get.log))
#     freqsBySub <- lapply(spectra, function(x) x$freq)
#     powerBySub <- lapply(spectra, function(x) x$spec)
#     
#     # not all subjects have the same trial numbers.
#     # As a result, the actual frequencies at which we have a power estimate differ per subject, so we need to interpolate
#     approximators <- mapply(function(x,y) approxfun(x=x,y=y), freqsBySub, powerBySub)
#     frequencies <- freqsBySub[[which.min(sapply(freqsBySub, length))]]    # get the frequencies corresponding to the shortest range(?)
#     interpolatedPowers <- do.call(rbind, lapply(approximators, function(x) x(frequencies))) # this has shape [nSubjects x nFrequencies]
#     
#     meanPower <- apply(interpolatedPowers, 2, mean)
#     return(data.frame(freq=frequencies, power=meanPower))
#   }
# }
# 
# plotSpectrum <- function(dat, pp, pp2=NULL,
#                          dat_fn=NULL, pp_fn=NULL, pp2_fn=NULL,
#                          xlab='', ylab='', main='', add.legend=FALSE, ylim=NULL, plot.log=TRUE, xlim=NULL, detrend=FALSE,
#                          demean=TRUE,
#                          trial_duration=NULL, full_x=TRUE, plot_xticklabels=TRUE, overwrite=TRUE) {
#   powerSpectraData <- load_save_spectrum(dat, dat_fn, detrend=detrend, demean=demean, overwrite=overwrite)
#   
#   if(plot.log) { f <- function(x) log(x) } else { f <- function(x) x }
#   if(is.null(trial_duration)) {
#     plot(x=f(powerSpectraData$freq), y=f(powerSpectraData$power), xlab=xlab, ylab=ylab, type='n', main=main, ylim=ylim, xlim=xlim)
#   } else {
#     # frequencies in Hz
#     x_axis_ticks_seconds <- c(3600, 30*60, 15*60, 5*60, 2*60, 60, 30, 5, 1)
#     x_axis_ticks_hz <- 1/x_axis_ticks_seconds
#     x_axis_ticks_1_over_trial <- x_axis_ticks_hz*trial_duration    ## NB
#     # ylim <- range(f(powerSpectraData$power))
#     
#     # which()
#     powers <- f(powerSpectraData$power)
#     freqs <- f(powerSpectraData$freq)
#     
#     if(!full_x) {
#       idx <- freqs>=(log(x_axis_ticks_1_over_trial)[4]-0.2)
#     } else {
#       idx <- rep(TRUE, length(powerSpectraData$freq))
#     }
#     #    xlim <- c(log(x_axis_ticks_1_over_trial)[4], max(freqs))
#     #    ylim <- range(powers[freqs>=xlim[1]], na.rm=TRUE)
#     
#     plot(x=f(powerSpectraData$freq)[idx], y=f(powerSpectraData$power)[idx], xlab=xlab, ylab=ylab, type='n', main=main, ylim=ylim, xlim=xlim, xaxt='n')
#     if(plot_xticklabels) {
#       axis(side=1, at=log(x_axis_ticks_1_over_trial), labels=c('1 hr', '30 m', '15 m', '5 m', '2 m', '1 m', '30 s', '5 s', '1 s'), las=2)
#     } else {
#       axis(side=1, at=log(x_axis_ticks_1_over_trial), labels=rep('', length(x_axis_ticks_1_over_trial)), las=2)
#     }
#   }
#   
#   
#   if(!is.null(pp)) {
#     powerSpectraModel <- load_save_spectrum(pp, pp_fn, by.post=TRUE, detrend=detrend, demean=demean, overwrite=overwrite)
#     #    powerSpectraModel <- powerSpectraModel(pp, by.postn=TRUE, detrend=detrend, demean=demean)
#     for(i in 1:nrow(powerSpectraModel$power)) {
#       lines(x=f(powerSpectraModel[[1]]),
#             y=f(powerSpectraModel[[2]][i,]),
#             col=adjustcolor(3, alpha.f=.1))
#     }
#     
#     if(!is.null(pp2)) {
#       powerSpectraModel2 <- load_save_spectrum(pp2, pp2_fn, by.post=TRUE, detrend=detrend, demean=demean)
#       #      powerSpectraModel2 <- getPowerSpectra(pp2, by.postn=TRUE, detrend=detrend, demean=demean)
#       for(i in 1:nrow(powerSpectraModel2$power)) {
#         lines(x=f(powerSpectraModel2[[1]]),
#               y=f(powerSpectraModel2[[2]][i,]),
#               col=adjustcolor(3, alpha.f=.1))
#       }
#     }
#     
#     ## get mean of pp
#     freqs <- powerSpectraModel[[1]]
#     powers <- apply(powerSpectraModel[[2]],2,mean)
#     lines(x=f(freqs), y=f(powers), col='dark red')
#     
#     if(!is.null(pp2)) {
#       freqs2 <- powerSpectraModel2[[1]]
#       powers2 <- apply(powerSpectraModel2[[2]],2,mean)
#       lines(x=f(freqs2), y=f(powers2), col='dark red')
#     }
#   }
#   lines(x=f(powerSpectraData$freq), y=f(powerSpectraData$power), col=par()$col)
#   
#   if(!is.null(pp)) lines(x=f(freqs), y=f(powers), col='dark green')
#   if(!is.null(pp2)) lines(x=f(freqs2), y=f(powers2), col='dark red')
# }

## Post-error effects
plotPES <- function(data_PES=NULL,
                    pp_PES_CI=NULL,
                    mean_rt=NULL,
                    xlab='Trial relative to error',
                    rtmeasure='rtabs',
                    main='', ylim=NULL,
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
         data_PES[data_PES$measure==rtmeasure,'value'], pch=4, lwd=1)
  
  ## posterior predictives
  arrows(x0=pp_PES_CI[pp_PES_CI$measure==rtmeasure,'trialNposterror'],
         y0=pp_PES_CI[pp_PES_CI$measure==rtmeasure,'value'][,1],
         y1=pp_PES_CI[pp_PES_CI$measure==rtmeasure,'value'][,3],
         angle=90, code=3, length=.02, lwd=1, col=1)
  xs <- pp_PES_CI[pp_PES_CI$measure==rtmeasure,'trialNposterror']
  lines(xs, pp_PES_CI[pp_PES_CI$measure==rtmeasure,'value'][,2], col=1)
  legend('topleft', c('data', 'model'), lwd=c(NA,1), lty=c(NA,1), col=c(1,1), pch=c(4,NA), bty='n')
}

plot_error_effects <- function(dat, pp, cores=10) {
  if(!'accuracy' %in% colnames(dat)) dat$accuracy <- dat$S==dat$R
  if(!'accuracy' %in% colnames(pp)) pp$accuracy <- pp$S==pp$R
  
  data_PES <- getErrorEffects(dat)
  pp_PES <- getErrorEffects(pp, mc.cores=cores)
  
  par(mfrow=c(1,1))
  plotPES(data_PES=data_PES$average, pp_PES=pp_PES$average,
          mean_rt=mean(aggregate(rt~subjects,dat,mean)[,2]),
          main='')
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
