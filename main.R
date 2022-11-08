
source('./functions/ConfigureEnv.R')
ConfigureEnv()


################################################################################
## load data

source("./functions/load_data.R")


################################################################################
## get forecasts at given lead time and threshold

source("./functions/get_forecasts.R")

lead <- 3
threshold <- 20
states <- alt_groups # choose states

## observations
tr_y <- as.numeric(tr_obs[lead, ] > threshold)
ts_y <- as.numeric(ts_obs[lead, ] > threshold)

## raw
pEns <- get_ens_probability(ts_fcst, lead, threshold)

## climatology
pClim <- get_clim_probability(tr_obs[1, ], n_test, threshold, disc_fcsts = discrete_fcsts)

## conditional climatology
pClim_cond <- get_clim_probability(tr_obs[1, ], n_test, threshold, disc_fcsts = discrete_fcsts, states = states)

## post-processed (logistic regression)
pPP <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lead, threshold, disc_fcsts = discrete_fcsts)

## conditional post-processed (logistic regression)
pPP_cond <- get_pp_probability(tr_y, ts_y, tr_fcst, ts_fcst, lead, threshold, disc_fcsts = discrete_fcsts, states = states)


################################################################################
## evaluate at given lead time and threshold

source("./functions/helper_functions.R")

##### brier score decompositions

## unconditional decomposition

bs_decomp(pClim, ts_y) # climatology
bs_decomp(pClim_cond, ts_y) # conditional climatology
bs_decomp(pEns, ts_y) # raw ensemble
bs_decomp(pPP, ts_y) # post-processed
bs_decomp(pPP_cond, ts_y) # conditional post-processed


## conditional decomposition 

bs_decomp_cond(pClim, ts_y, states$ts) # climatology
bs_decomp_cond(pClim_cond, ts_y, states$ts) # conditional climatology
bs_decomp_cond(pEns, ts_y, states$ts) # raw ensemble
bs_decomp_cond(pPP, ts_y, states$ts) # post-processed
bs_decomp_cond(pPP_cond, ts_y, states$ts) # conditional post-processed


##### figures

## reliability diagrams

reldiag(pClim, ts_y, region.position = "diagonal") # climatology
reldiag(pClim_cond, ts_y) # conditional climatology
reldiag(pEns, ts_y) # raw ensemble
reldiag(pPP, ts_y) # post-processed
reldiag(pPP_cond, ts_y) # conditional post-processed


## waterfall plots

sc <- 10000 # scale all score decomposition terms by 10000 for interpretability

waterfall_plot(pClim, ts_y, states$ts, scale = sc) # climatology
waterfall_plot(pClim_cond, ts_y, states$ts, scale = sc) # conditional climatology
waterfall_plot(pEns, ts_y, states$ts, scale = sc) # raw ensemble
waterfall_plot(pPP, ts_y, states$ts, scale = sc) # post-processed
waterfall_plot(pPP_cond, ts_y, states$ts, scale = sc) # conditional post-processed



################################################################################
## evaluate at all lead times

threshold <- 20

## classical decomposition

plot_vs_lt(threshold, method = "COSMO", states = states)


## conditional decomposition

plot_vs_lt(threshold, method = "COSMO", states = states, classical = F)


################################################################################
## evaluate at all thresholds

threshold_vec <- seq(8, 32, 1)
lead <- 3

## classical decomposition

plot_vs_t(threshold_vec, lead, method = "COSMO", states = states)


## conditional decomposition

plot_vs_t(threshold_vec, lead, method = "COSMO", states = states, classical = F)


## compare theoretical and actual improvements of post-processing

plot_pp_gain(threshold_vec, lead)


## compare theoretical and actual improvements of conditional post-processing

plot_cond_pp_gain(threshold_vec, lead, states)


################################################################################
## analysis of information

terms_unc <- bs_decomp(pEns, ts_y)
terms_cnd <- bs_decomp_cond(pEns, ts_y, states$ts)

# proportion of information explained by the forecasts (generalised R^2 value)
(terms_unc['RES'] - terms_unc['REL'])/terms_unc['UNC']

# potential proportion of information explained by the forecasts
terms_unc['RES']/terms_unc['UNC']

# proportion of information lost due to miscalibration
terms_unc['REL']/terms_unc['UNC']

# proportion of information explained by the states
terms_cnd['RES_A']/terms_unc['UNC']

# proportion of information explained by the states that is not captured by the forecast (generalised partial R^2 value)
terms_cnd['RES_A|F']/(terms_unc['UNC'] - terms_unc['RES'])


