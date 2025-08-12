RS <- function(d) {
  ## which stimulus *symbol* was chosen?
  ifelse(d$R=='left', d$s_left, d$s_right)
}

Smatch_prereversal <- function(d) {
  # this figures out which responses correspond to the symbol that had the high pay-off probability at the start of the experiment.
  d$out <- NA
  for(subject in unique(d$subjects)) {
    tmp <- d[d$subjects==subject,]
    tmp <- tmp[!duplicated(tmp[,c('s_left', 'p_left')]),]
    tmp$correct_symbols <- ifelse(tmp$p_left>tmp$p_right, tmp$s_left, tmp$s_right)
    correct_symbols_prereversal <- tmp[!duplicated(tmp$s_left),'correct_symbols']
    d[d$subjects==subject,'out'] <- d[d$subjects==subject,'RS'] %in% correct_symbols_prereversal
  }
  d$out
}

# For RL paradigms, we need to use a combination of response-coding of accumulators (left/right)
# and *stimulus/symbol-coding* of accumulators (mapping each accumulator onto a symbol that represents each choice option),
lS <- function(d) {
  ## For race models: In dadms, there's a column 'lR', which codes which latent response (left/right) each accumulator corresponds to.
  ## In RL paradigms we also need to know to which latent response *symbol* the accumulator corresponds. That's lRS
  factor(d[cbind(1:nrow(d), match(paste0('s_', d$lR), colnames(d)))])
}

Scorrect <- function(d) {
  ## which stimulus *symbol* was correct?
  ifelse(d$p_left > d$p_right, d$s_left, d$s_right)
}

# this is a function generator that generates one column per symbol, with the appropriate structure in the dadm for updating.
covariate_column_generator <- function(col_name) {
  function(d) {
    d[,col_name] <- NA
    d[d$RS == col_name, col_name] <- d[d$RS == col_name, 'reward']
    #d[d$lR!=levels(d$lR)[1],col_name] <- NA
    d[d$RS != d$lS, col_name] <- NA
    d[[col_name]]
  }
}
