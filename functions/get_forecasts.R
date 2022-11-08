################################################################################
## get forecasts

# function to map forecasts to a discrete number of possible options
convert_to_discrete <- function(f, disc_f){
  n_bins <- length(disc_f)
  p <- sapply(seq_along(f), function(i) disc_f[which.min(abs(f[i] - disc_f))])
  return(p)
}


# function to extract probability forecast from an ensemble
get_ens_probability <- function(ens, lead, threshold){
  n_ens <- dim(ens)[3]
  p <- (rowSums(ens[lead, , ] > threshold) + 1)/(n_ens + 2)
  return(p)
}


# function to extract climatological probability forecasts
get_clim_probability <- function(obs, n_test, threshold, disc_fcsts, states = NULL){
  if(is.null(states)){
    p <- mean(obs > threshold, na.rm = T)
    p <- rep(p, n_test)
  }else{
    n_states <- length(unique(states$tr))
    p_states <- sapply(1:n_states, function(i) mean(obs[states$tr == i] > threshold, na.rm = T))
    p <- numeric(n_test)
    for(i in 1:n_states) p[states$ts == i] <- p_states[i]
  }
  p <- convert_to_discrete(p, disc_fcsts)
  return(p)
}
  

# function to extract post processed probability forecasts
get_pp_probability <- function(tr_y, ts_y, tr_ens, ts_ens, lead, threshold, disc_fcsts, states = NULL){
  tr_p <- get_ens_probability(tr_ens, lead, threshold)
  ts_p <- get_ens_probability(ts_ens, lead, threshold)
  if(is.null(states)){
    train_df <- data.frame(obs = tr_y, ens = tr_p)
    PPmodel <- glm(obs ~ ens, data = train_df, family = binomial)
    test_df <- data.frame(obs = ts_y, ens = ts_p)
    p <- predict(PPmodel, test_df, type = "response")
  }else{
    train_df <- data.frame(obs = tr_y, ens = tr_p, st = states$tr)
    PPmodel <- glm(obs ~ ens + st, data = train_df, family = binomial)
    test_df <- data.frame(obs = ts_y, ens = ts_p, st = states$ts)
    p <- predict(PPmodel, test_df, type = "response")
  }
  p <- convert_to_discrete(p, disc_fcsts)
  return(p)
}

