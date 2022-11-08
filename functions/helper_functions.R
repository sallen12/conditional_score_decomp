################################################################################
# helper functions


# function to obtain unconditional Brier score decomposition terms from isotonic regression
bs_decomp_iso <- function(p, o, states = NULL){
  na_ind <- is.na(o)
  
  diag_obj <- reliabilitydiag::reliabilitydiag(x = p[!na_ind], y = o[!na_ind])
  terms <- c(summary(diag_obj)$uncertainty, 
             summary(diag_obj)$discrimination,
             summary(diag_obj)$miscalibration,
             summary(diag_obj)$mean_score)
  names(terms) <- c("UNC", "RES", "REL", "TOT")
  return(terms)
}


# function to obtain unconditional Brier score decomposition terms
bs_decomp <- function(p, o, method = "bias_corrected"){

  p <- p[!is.na(o)] 
  o <- o[!is.na(o)]

  if(method == "isotonic"){
    
    terms <- bs_decomp_iso(p, o)
    
  }else if(method %in% c("classical", "bias_corrected")){
    
    forecasts <- unique(p)
    N <- length(o) 
    obar <- mean(o)
    
    n_k <- sapply(seq_along(forecasts), function(k) sum(p == forecasts[k]))
    o_k <- sapply(seq_along(forecasts), function(k) sum(o[p == forecasts[k]]))
    obar_k <- sapply(seq_along(forecasts), function(k) mean(o[p == forecasts[k]]))
    
    unc <- obar*(1 - obar)
    res <- sum((n_k/N)*((obar_k - obar)^2))
    rel <- sum((n_k/N)*((obar_k - forecasts)^2))
    
    if(method == "bias_corrected"){
      unc <- unc + unc/(N - 1)
      res <- res + unc/(N - 1) - sum(n_k*obar_k*(1 - obar_k)/(n_k - 1))/N
      rel <- rel - sum(n_k*obar_k*(1 - obar_k)/(n_k - 1))/N
    }
    
    tot <- unc - res + rel
    terms <- c(unc, res, rel, tot)
    names(terms) <- c("UNC", "RES", "REL", "TOT")
    
  }
  
  return(terms)

}


# function to obtain conditional Brier score decomposition terms  
bs_decomp_cond <- function(p, o, states, method = "bias_corrected"){
  
  states <- states[!is.na(o)]
  p <- p[!is.na(o)] 
  o <- o[!is.na(o)]

  groups <- unique(states)
  n_j <- sapply(seq_along(groups), function(j) sum(states == groups[j]))
  N <- length(o) 

  if(method == "isotonic"){

    terms_un <- bs_decomp_iso(p, o) # unconditional terms
    terms_cnd <- sapply(seq_along(groups), function(j) bs_decomp_iso(p[states == j], o[states == j])) # conditional terms
    
    unc_A <- sum((n_j/N)*terms_cnd['UNC', ])
    res_fA <- sum((n_j/N)*terms_cnd['RES', ])
    rel_fA <- sum((n_j/N)*terms_cnd['REL', ])
    
    res_A <- terms_un['UNC'] - unc_A
    res_Af <- rel_fA - terms_un['REL']
    
  }else if(method %in% c("classical", "bias_corrected")){
    
    forecasts <- unique(p)
    obar <- mean(o)
    
    n_k <- sapply(seq_along(forecasts), function(k) sum(p == forecasts[k]))
    o_k <- sapply(seq_along(forecasts), function(k) sum(o[p == forecasts[k]]))
    obar_k <- sapply(seq_along(forecasts), function(k) mean(o[p == forecasts[k]]))
    
    o_j <- sapply(seq_along(groups), function(j) sum(o[states == groups[j]]))
    obar_j <- sapply(seq_along(groups), function(j) mean(o[states == groups[j]]))
    
    n_kj <- sapply(seq_along(forecasts), function(k) 
      sapply(seq_along(groups), function(j) sum(p == forecasts[k] & states == groups[j])))
    o_kj <- sapply(seq_along(forecasts), function(k) 
      sapply(seq_along(groups), function(j) sum(o[p == forecasts[k] & states == groups[j]])))
    obar_kj <- sapply(seq_along(forecasts), function(k) 
      sapply(seq_along(groups), function(j) mean(o[p == forecasts[k] & states == groups[j]])))
    
    unc_A <- sum((n_j/N)*obar_j*(1 - obar_j))
    res_A <- sum((n_j/N)*((obar_j - obar)^2))
    res_fA <- sum((n_kj/N)*((obar_j - obar_kj)^2), na.rm = T)
    res_Af <- sum(t(n_kj/N)*((obar_k - t(obar_kj))^2), na.rm = T)
    rel_fA <- sum(t(n_kj/N)*((forecasts - t(obar_kj))^2), na.rm = T)
    tot <- unc_A - res_fA + rel_fA
    
    if(method == "bias_corrected"){
      bias1 <- sum(n_j*obar_j*(1 - obar_j)/(n_j - 1))/N
      bias2 <- obar*(1 - obar)/(N - 1)
      bias3 <- sum(n_kj*obar_kj*(1 - obar_kj)/(n_kj - 1), na.rm = T)/N
      bias4 <- sum(n_k*obar_k*(1 - obar_k)/(n_k - 1))/N
      
      unc_A <- unc_A + bias1
      res_A <- res_A - bias1 + bias2
      res_fA <- res_fA + bias1 - bias3
      res_Af <- res_Af - bias3 + bias4
      rel_fA <- rel_fA - bias3
    }
    
  }
  
  tot <- unc_A - res_fA + rel_fA
  terms <- c(unc_A, res_A, res_fA, res_Af, rel_fA, tot)
  names(terms) <- c("UNC_A", "RES_A", "RES_F|A", "RES_A|F", "REL_F|A", "TOT")
  
  return(terms)
  
}


# function to plot reliability diagrams in the presence of NAs
reldiag <- function(p, o, level = 0.95, region.position = "estimate"){
  na_ind <- is.na(o)
  reliabilitydiag(p = p[!na_ind], y = o[!na_ind], region.position = region.position, region.level = level)
}


# function to create a waterfall plot
waterfall_plot <- function(p, o, states, title = "", scale = 10000){
  
  # plot parameters may need adjusting for scales different from 10000
  if(is.null(scale)){
    scale <- 1
  }
  
  # unconditional
  mth_cl <- bs_decomp(p, o)*scale
  mth_cl[2] <- -mth_cl[2]
  df <- data.frame(begin = c(0, cumsum(mth_cl[1:2]), 0), end = mth_cl + c(0, cumsum(mth_cl[1:2]), 0), id = c(1.5, 4, 6.5, 8),
                   names = c(" UNC_Y", "-RES_F", " REL_F", " Total"),
                   cols = as.factor(c(1, 2, 3, 4)))
  plot_cl <- ggplot(df) + geom_rect(aes(ymin = c(1, 3, 6, 8) - 0.45, ymax = c(2, 5, 7, 8) + 0.45, 
                                        xmin = 0, xmax = end - begin, fill = cols)) +
    geom_text(aes(x = 0.35*scale, y = id, label = sprintf("%.1f", mth_cl), fontface = 2), 
              colour = c(scales::hue_pal()(3), "black"), hjust = 1) +
    scale_fill_manual(values = c(scales::hue_pal()(3), "black")) + 
    scale_x_continuous(name = NULL, limits = c(-0.15, 0.35)*scale, breaks = seq(-0.1, 0.2, 0.1)*scale) +
    scale_y_continuous(name = NULL, breaks = df$id, position = "right",
                       labels = c(expression(paste("   ", UNC[Y])), 
                                  expression(-RES["F"]), 
                                  expression(paste("   ", REL["F"])),
                                  paste(" ", "Total"))) +
    geom_vline(aes(xintercept = 0)) +
    ggtitle("") + theme_classic() + 
    theme(legend.position = "none", axis.text = element_text(size = 11), axis.ticks.y = element_blank())
  
  # conditional
  mth <- bs_decomp_cond(p, o, states)*scale
  mth <- c(mth[1], mth[2], -mth[2], -mth[3], mth[4], -mth[4], mth[5], mth[6])
  df <- data.frame(begin = c(0, cumsum(mth[1:6]), 0), end = mth + c(0, cumsum(mth[1:6]), 0), id = 1:8,
                   names = c("UNC_Y|A", "RES_A", "-RES_A", "-RES_F|A", "RES_A|F", "-RES_A|F", "REL_F|A", "Total"),
                   cols = as.factor(c(1, 1, 2, 2, 2, 3, 3, 4)))
  plot_new <- ggplot(df) + geom_rect(aes(ymin = id - 0.45, ymax = id + 0.45, xmin = 0, xmax = end - begin, fill = cols)) +
    geom_text(aes(x = -0.15*scale, y = id, label = sprintf("%.1f", mth), fontface = 2), 
              colour = c(scales::hue_pal()(3), "black")[df$cols], hjust = 1) +
    scale_fill_manual(values = c(scales::hue_pal()(3), "black")) + 
    scale_x_continuous(name = NULL, limits = c(-0.25, 0.25)*scale, breaks = seq(-0.1, 0.2, 0.1)*scale, expand = c(0, 0)) +
    scale_y_continuous(name = NULL, breaks = df$id, 
                       labels = c(expression(paste("   ", UNC["Y|A"])),
                                  expression(paste("   ", RES[A])),
                                  expression(-RES[A]), 
                                  expression(-RES["F|A"]), 
                                  expression(paste("   ", RES["A|F"])),
                                  expression(-RES["A|F"]),
                                  expression(paste("   ", REL["F|A"])), 
                                  paste("  ", "Total"))) +
    geom_vline(aes(xintercept = 0)) +
    ggtitle(title) + theme_classic() + 
    theme(legend.position = "none", axis.text = element_text(size = 11), 
          axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0))
  
  gridExtra::grid.arrange(plot_new, plot_cl, nrow = 1)
}


# function to create an alternative waterfall plot
waterfall_plot_alt <- function(p, o, states, title = "", scale = 10000){

  # plot parameters may need adjusting for scales different from 10000
  if(is.null(scale)){
    scale <- 1
  }
  
  # unconditional
  mth_cl <- bs_decomp(p, o)*scale
  mth_cl[2] <- -mth_cl[2]
  df <- data.frame(begin = c(0, cumsum(mth_cl[1:2]), 0), end = mth_cl + c(0, cumsum(mth_cl[1:2]), 0), 
                   id = c(1.5, 4, 6.5, 8),
                   names = c("UNC_Y", "-RES_F", "REL_F", "Total"),
                   cols = as.factor(c(1, 2, 3, 4)))
  plot_cl <- ggplot(df) + geom_rect(aes(ymin = c(1, 3, 6, 8) - 0.45, ymax = c(2, 5, 7, 8) + 0.45, 
                                        xmin = begin, xmax = end, fill = cols)) +
    geom_text(aes(x = (begin + end)/2, y = id, label = sprintf("%.1f", mth_cl)), colour = c(rep("black", 3), "white")) +
    geom_vline(aes(xintercept = 0)) + 
    geom_linerange(aes(x = end, ymin = c(2, 5, 7, 8) + 0.45, ymax = c(2, 5, 7, 8) + 1 - 0.45), 
                   size = 1, colour = c(rep("black", 3), "white")) +
    scale_fill_manual(values = c(scales::hue_pal()(3), "black")) + 
    scale_x_continuous(name = NULL, limits = c(0, 2600), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, breaks = df$id, position = "right",
                       labels = c(expression(UNC[Y]), expression(-RES["F"]), expression(REL["F"]), "Total")) +
    ggtitle("") + theme_classic() + 
    theme(legend.position = "none", axis.text = element_text(size = 11))
  
  # conditional
  mth <- bs_decomp_cond(p, o, states)*scale
  mth <- c(mth[1], mth[2], -mth[2], -mth[3], mth[4], -mth[4], mth[5], mth[6])
  df <- data.frame(begin = c(0, cumsum(mth[1:6]), 0), end = mth + c(0, cumsum(mth[1:6]), 0), id = 1:8,
                   names = c("UNC_Y|A", "RES_A", "-RES_A", "-RES_F|A", "RES_A|F", "-RES_A|F", "REL_F|A", "Total"),
                   cols = as.factor(c(1, 1, 2, 2, 2, 3, 3, 4)))
  plot_new <- ggplot(df) + geom_rect(aes(ymin = id - 0.45, ymax = id + 0.45, xmin = begin, xmax = end, fill = cols)) +
    geom_text(aes(x = (begin + end)/2, y = id, label = sprintf("%.1f", mth)), colour = c(rep("black", 7), "white")) +
    geom_vline(aes(xintercept = 0)) +
    geom_linerange(aes(x = end, ymin = id + 0.45, ymax = id + 1 - 0.45), 
                   size = 1, colour = c(rep("black", 7), "white")) +
    scale_fill_manual(values = c(scales::hue_pal()(3), "black")) + 
    scale_x_continuous(name = NULL, limits = c(0, 2600), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, breaks = df$id, 
                       labels = c(expression(UNC["Y|A"]), expression(RES[A]), expression(-RES[A]), 
                                  expression(-RES["F|A"]), expression(RES["A|F"]), expression(-RES["A|F"]),
                                  expression(REL["F|A"]), "Total")) +
    ggtitle(title) + theme_classic() + 
    theme(legend.position = "none", axis.text = element_text(size = 11))
  
  gridExtra::grid.arrange(plot_new, plot_cl, nrow = 1)
}


# function to plot decomposition terms against lead time
plot_vs_lt <- function(threshold, method = "COSMO", states, unconditional = T){
  
  if(unconditional){
    if(method == "COSMO"){
      p_decomp <- sapply(1:n_lt, function(lt){
        ts_y <- as.numeric(ts_obs[lt, ] > threshold)
        pEns <- get_ens_probability(ts_fcst, lt, threshold)
        bs_decomp(pEns, ts_y)*sc 
      }) # raw ensemble
    }else if(method == "clim"){
      p_decomp <- sapply(1:n_lt, function(lt){
        ts_y <- as.numeric(ts_obs[lt, ] > threshold)
        pClim <- get_clim_probability(tr_obs[1, ], n_test, threshold, disc_fcsts = discrete_fcsts)
        bs_decomp(pClim, ts_y)*sc 
      }) # climatology
    }else if(method == "cond_clim"){
      p_decomp <- sapply(1:n_lt, function(lt){
        ts_y <- as.numeric(ts_obs[lt, ] > threshold)
        pClim_cond <- get_clim_probability(tr_obs[1, ], n_test, threshold, disc_fcsts = discrete_fcsts, states = states)
        bs_decomp(pClim_cond, ts_y)*sc 
      }) # conditional climatology
    }else if(method == "post-proc"){
      p_decomp <- sapply(1:n_lt, function(lt){
        tr_y <- as.numeric(tr_obs[lt, ] > threshold)
        ts_y <- as.numeric(ts_obs[lt, ] > threshold)
        pPP <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lt, threshold, disc_fcsts = discrete_fcsts)
        bs_decomp(pPP, ts_y)*sc 
      }) # post-processed
    }else if(method == "cond_post-proc"){
      p_decomp <- sapply(1:n_lt, function(lt){
        tr_y <- as.numeric(tr_obs[lt, ] > threshold)
        ts_y <- as.numeric(ts_obs[lt, ] > threshold)
        pPP_cond <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lt, threshold, disc_fcsts = discrete_fcsts, states = states)
        bs_decomp(pPP_cond, ts_y)*sc 
      }) # conditional post-processed
    }
    
    df <- data.frame(lead = rep(1:n_lt, each = 4), term = rep(c("UNC", "RES", "REL", " Tot"), n_lt),
                     s = as.vector(p_decomp))
    ggplot(df) + geom_line(aes(x = lead, y = s, col = term), size = 1) + 
      geom_hline(aes(yintercept = 0), lty = "dotted") + ylim(c(0, 2500)) +    
      scale_y_continuous(name = "") +
      scale_x_continuous(name = "Lead time (days)") +
      scale_color_manual(values = c("black", "#619CFF", "#00BA38", "#F8766D"), 
                         labels = c("Tot", expression(REL[F]), expression(RES[F]), expression(UNC[Y]))) +
      theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom")
  }else{
  if(method == "COSMO"){
    p_decomp <- sapply(1:n_lt, function(lt){
      ts_y <- as.numeric(ts_obs[lt, ] > threshold)
      pEns <- get_ens_probability(ts_fcst, lt, threshold)
      bs_decomp(pEns, ts_y, states$ts)*sc 
    }) # raw ensemble
  }else if(method == "clim"){
    p_decomp <- sapply(1:n_lt, function(lt){
      ts_y <- as.numeric(ts_obs[lt, ] > threshold)
      pClim <- get_clim_probability(tr_obs[1, ], n_test, threshold, disc_fcsts = discrete_fcsts)
      bs_decomp(pClim, ts_y, states$ts)*sc
    }) # climatology
  }else if(method == "cond_clim"){
    p_decomp <- sapply(1:n_lt, function(lt){
      ts_y <- as.numeric(ts_obs[lt, ] > threshold)
      pClim_cond <- get_clim_probability(tr_obs[1, ], n_test, threshold, disc_fcsts = discrete_fcsts, states = states)
      bs_decomp(pClim_cond, ts_y, states$ts)*sc
    }) # conditional climatology
  }else if(method == "post-proc"){
    p_decomp <- sapply(1:n_lt, function(lt){
      tr_y <- as.numeric(tr_obs[lt, ] > threshold)
      ts_y <- as.numeric(ts_obs[lt, ] > threshold)
      pPP <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lt, threshold, disc_fcsts = discrete_fcsts)
      bs_decomp(pPP, ts_y, states$ts)*sc
    }) # post-processed
  }else if(method == "cond_post-proc"){
    p_decomp <- sapply(1:n_lt, function(lt){
      tr_y <- as.numeric(tr_obs[lt, ] > threshold)
      ts_y <- as.numeric(ts_obs[lt, ] > threshold)
      pPP_cond <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lt, threshold, disc_fcsts = discrete_fcsts, states = states)
      bs_decomp(pPP_cond, ts_y, states$ts)*sc
    }) # conditional post-processed
  }

  df <- data.frame(lead = rep(1:n_lt, each = 6), s = as.vector(p_decomp),
                   term = rep(c("UNC_Y|A", "RES_A", "RES_F|A", "RES_A|F", "REL_F|A", " Tot"), n_lt))
  ggplot(df) + geom_line(aes(x = lead, y = s, col = term), size = 1) +
    geom_hline(aes(yintercept = 0), lty = "dotted") +
    scale_y_continuous(name = "", limits = c(0, 2500)) +
    scale_x_continuous(name = "Lead time (days)") + 
    scale_color_manual(values = c("black", "#D39200", "#00C19F", "#DB72FB", "#93AA00", "#00B9E3"), 
                       labels = c("Tot      ", expression(REL["F|A"]), expression(RES[A  ]),
                                  expression(RES["A|F"]),
                                  expression(RES["F|A"]), expression(UNC["Y|A"]))) +
    theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom")
  }
}


# function to plot decomposition terms against threshold
plot_vs_t <- function(threshold_vec, lead, method = "COSMO", states, unconditional = T){
  
  if(unconditional){
    if(method == "COSMO"){
      p_decomp <- sapply(seq_along(threshold_vec), function(i){
        ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
        pEns <- get_ens_probability(ts_fcst, lead, threshold_vec[i])
        bs_decomp(pEns, ts_y)*sc
      }) # raw ensemble
    }else if(method == "clim"){
      p_decomp <- sapply(seq_along(threshold_vec), function(i){
        ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
        pClim <- get_clim_probability(tr_obs[1, ], n_test, threshold_vec[i], disc_fcsts = discrete_fcsts)
        bs_decomp(pClim, ts_y)*sc
      }) # climatology
    }else if(method == "cond_clim"){
      p_decomp <- sapply(seq_along(threshold_vec), function(i){
        ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
        pClim_cond <- get_clim_probability(tr_obs[1, ], n_test, threshold_vec[i], disc_fcsts = discrete_fcsts, states = states)
        bs_decomp(pClim_cond, ts_y)*sc
      }) # conditional climatology
    }else if(method == "post-proc"){
      p_decomp <- sapply(seq_along(threshold_vec), function(i){
        tr_y <- as.numeric(tr_obs[lead, ] > threshold_vec[i])
        ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
        pPP <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lead, threshold_vec[i], disc_fcsts = discrete_fcsts)
        bs_decomp(pPP, ts_y)*sc
      }) # post-processed
    }else if(method == "cond_post-proc"){
      p_decomp <- sapply(seq_along(threshold_vec), function(i){
        tr_y <- as.numeric(tr_obs[lead, ] > threshold_vec[i])
        ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
        pPP_cond <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lead, threshold_vec[i], disc_fcsts = discrete_fcsts, states = states)
        bs_decomp(pPP_cond, ts_y)*sc 
      }) # conditional post-processed
    }
    
    
    df <- data.frame(x = as.vector(p_decomp),  t = rep(threshold_vec, each = 4),
                     y = rep(c("UNC", "RES", "REL", " Tot"), length(threshold_vec)))
    ggplot(df) + geom_line(aes(x = t, y = x, col = y), size = 1) +
      geom_hline(aes(yintercept = 0), lty = "dotted") +
      scale_x_continuous(name = "Threshold (\u00B0C)") +
      scale_y_continuous(name = "", limits = c(0, 2500)) + 
      scale_color_manual(values = c("black", "#619CFF", "#00BA38", "#F8766D"), 
                         labels = c("Tot", expression(REL[F]), expression(RES[F]), expression(UNC[Y]))) +
      theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom") 
      
  }else{
    if(method == "COSMO"){
      p_decomp <- sapply(seq_along(threshold_vec), function(i){
        ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
        pEns <- get_ens_probability(ts_fcst, lead, threshold_vec[i])
        bs_decomp(pEns, ts_y, states$ts)*sc
      }) # raw ensemble
    }else if(method == "clim"){
      p_decomp <- sapply(seq_along(threshold_vec), function(i){
        ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
        pClim <- get_clim_probability(tr_obs[1, ], n_test, threshold_vec[i], disc_fcsts = discrete_fcsts)
        bs_decomp(pClim, ts_y, states$ts)*sc
      }) # climatology
    }else if(method == "cond_clim"){
      p_decomp <- sapply(seq_along(threshold_vec), function(i){
        ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
        pClim_cond <- get_clim_probability(tr_obs[1, ], n_test, threshold_vec[i], disc_fcsts = discrete_fcsts, states = states)
        bs_decomp(pClim_cond, ts_y, states$ts)*sc
      }) # conditional climatology
    }else if(method == "post-proc"){
      p_decomp <- sapply(seq_along(threshold_vec), function(i){
        tr_y <- as.numeric(tr_obs[lead, ] > threshold_vec[i])
        ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
        pPP <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lead, threshold_vec[i], disc_fcsts = discrete_fcsts)
        bs_decomp(pPP, ts_y, states$ts)*sc
      }) # post-processed
    }else if(method == "cond_post-proc"){
      p_decomp <- sapply(seq_along(threshold_vec), function(i){
        tr_y <- as.numeric(tr_obs[lead, ] > threshold_vec[i])
        ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
        pPP_cond <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lead, threshold_vec[i], disc_fcsts = discrete_fcsts, states = states)
        bs_decomp(pPP_cond, ts_y, states$ts)*sc 
      }) # conditional post-processed
    }
    
    df <- data.frame(x = as.vector(p_decomp),  t = rep(threshold_vec, each = 6),
                     y = rep(c("UNC_Y|A", "RES_A", "RES_F|A", "RES_A|F", "REL_F|A", " Tot"),
                             length(threshold_vec)))
    ggplot(df) + geom_line(aes(x = t, y = x, col = y), size = 1) + 
      geom_hline(aes(yintercept = 0), lty = "dotted") +
      scale_x_continuous(name = "Threshold (\u00B0C)") +
      scale_y_continuous(name = "", limits = c(0, 2500)) + 
      scale_color_manual(values = c("black", "#D39200", "#00C19F", "#DB72FB", "#93AA00", "#00B9E3"), 
                         labels = c("Tot      ", expression(REL["F|A"]), expression(RES[A  ]),
                                    expression(RES["A|F"]), expression(RES["F|A"]), expression(UNC["Y|A"]))) +
      theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom")
      
  }
}


# function to estimate confidence regions for skill scores using nonparametric bootstrapping
bss_unc <- function(s, s_ref, N = 1000, level = 0.95){
  n <- length(s)
  ss <- sapply(1:N, function(i){
    ind <- sample(1:n, n, replace = TRUE)
    skill <- 1 - mean(s[ind], na.rm = T)/mean(s_ref[ind], na.rm=T)
    return(skill)
  })
  ss_CI <- quantile(ss, c((1 - level)/2, (1 + level)/2))
  return(unname(ss_CI))
}


# function to estimate confidence regions for decomposition terms using nonparametric bootstrapping
rel_unc <- function(p, o, N = 1000, level = 0.95, states = NULL){
  n <- length(o)
  if(is.null(states)){
    rel <- sapply(1:N, function(i){
      ind <- sample(1:n, n, replace = TRUE)
      terms <- bs_decomp(p[ind], o[ind])
      impr <- terms["REL"]/terms["TOT"]
      return(unname(impr))
    })
  }else{
    rel <- sapply(1:N, function(i){
      ind <- sample(1:n, n, replace = TRUE)
      terms <- bs_decomp_cond(p[ind], o[ind], states = states[ind])
      impr <- terms["RES_A|F"]/terms["TOT"]
      return(unname(impr))
    })
  }
  rel_CI <- quantile(rel, c((1 - level)/2, (1 + level)/2))
  return(unname(rel_CI))
}


# plot improvement gained by the ensemble forecast as a function of threshold
plot_pp_gain <- function(threshold_vec, lead){
  
  ss_CI <- sapply(seq_along(threshold_vec), function(i){
    print(i)
    tr_y <- as.numeric(tr_obs[lead, ] > threshold_vec[i])
    ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
    pEns <- get_ens_probability(ts_fcst, lead, threshold_vec[i])
    pPP <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lead, threshold_vec[i], disc_fcsts = discrete_fcsts)
    
    # ss
    s_ens <- (pEns - ts_y)^2
    s_pp <- (pPP - ts_y)^2
    ss <- 1 - mean(s_pp, na.rm = T)/mean(s_ens, na.rm = T)
    ss_unc <- bss_unc(s_pp, s_ens)
    
    # rel
    impr <- bs_decomp(pEns, ts_y)
    impr <- unname(impr["REL"]/impr["TOT"])
    impr_unc <- rel_unc(pEns, ts_y)
    
    return(c(ss, ss_unc, impr, impr_unc))
  })
  
  df <- data.frame(t = threshold_vec, 
                   x = c(ss_CI[1, ], ss_CI[4, ]),
                   se_low = c(ss_CI[2, ], ss_CI[5, ]),
                   se_high = c(ss_CI[3, ], ss_CI[6, ]),
                   mth = rep(c("Actual", "Theoretical"), each = length(threshold_vec)))
  ggplot(df) + geom_line(aes(x = t, y = x, col = mth), size = 1) + 
    geom_hline(aes(yintercept = 0), lty = "dotted") +
    geom_ribbon(aes(x = t, ymin = se_low, ymax = se_high, fill = mth), alpha = 0.2) +
    scale_x_continuous(name = "Threshold (\u00B0C)") +
    scale_y_continuous(name = "Brier skill score", limits = c(-0.05, 0.35)) + 
    theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom") +
    guides(fill = "none")
}


# plot improvement gained by conditional post-processing as a function of threshold
plot_cond_pp_gain <- function(threshold_vec, lead, states){
  
  ss_CI <- sapply(seq_along(threshold_vec), function(i){
    print(i)
    tr_y <- as.numeric(tr_obs[lead, ] > threshold_vec[i])
    ts_y <- as.numeric(ts_obs[lead, ] > threshold_vec[i])
    pPP <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lead, threshold_vec[i], disc_fcsts = discrete_fcsts)
    pPP_cond <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lead, threshold_vec[i], disc_fcsts = discrete_fcsts, states = states)
    
    # ss
    s_pp <- (pPP - ts_y)^2
    s_pp_cond <- (pPP_cond - ts_y)^2
    ss <- 1 - mean(s_pp_cond, na.rm = T)/mean(s_pp, na.rm = T)
    ss_unc <- bss_unc(s_pp_cond, s_pp)
    
    # rel
    impr <- bs_decomp(pPP, ts_y, states = states$ts)
    impr <- unname(impr["RES_A|F"]/impr["TOT"])
    impr_unc <- rel_unc(pPP, ts_y, states = states$ts)
    
    return(c(ss, ss_unc, impr, impr_unc))
  })
  
  df <- data.frame(t = threshold_vec, 
                   x = c(ss_CI[1, ], ss_CI[4, ]),
                   se_low = c(ss_CI[2, ], ss_CI[5, ]),
                   se_high = c(ss_CI[3, ], ss_CI[6, ]),
                   mth = rep(c("Actual", "Theoretical"), each = length(threshold_vec)))
  ggplot(df) + geom_line(aes(x = t, y = x, col = mth), size = 1) +
    geom_hline(aes(yintercept = 0), lty = "dotted") +
    geom_ribbon(aes(x = t, ymin = se_low, ymax = se_high, fill = mth), alpha = 0.2) +
    scale_x_continuous(name = "Threshold (\u00B0C)") +
    scale_y_continuous(name = "Brier skill score", limits = c(-0.05, 0.35)) +
    theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom") +
    guides(fill = "none")
}
