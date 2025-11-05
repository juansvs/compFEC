library(ggplot2)
# Pooled Fecal Egg Count simulation
# Juan S. Vargas
# November 2025
#

#' Code to simulate collection and analysis of fecal samples from goats. The
#' simulation focuses on understanding the differences between analyzing samples
#' individually and as a pooled sample
library(tidyverse)
library(ggpubr)

# set.seed(345)
rm(list = ls())
gen_cv <- function(n) VGAM::rtriangle(n, lower = 0.2, theta = 0.3, upper = 0.73)# coefficient of variation in weight samples, random value from a triangular dist. (0.2, 0.3, 0.73)

niter <- 100            # number of iterations per scenario
Na <- 30                 # number of animals
ms <- c(100, 500, 2000)  # population mean EPG values
ks <- c(0.1, 0.5, 2, 10) # population EPG aggregation values
targ_wt <- 5             # sample (farm level) target weight
subsamp_wts <- c(0.5,1,2) # subsample weights for composite (at the lab level)
subsamp_sds <- 1/1.96*c(0.1,0.2,0.5,1) # subsample weight sds: tolerance in g divided by 1.96 (95%CI of normal dist)
flot_wt <- 2      # subsample weight for indiv. FEC. 2 g for McMaster technique
dls <- c(5, 10, 20, 50)  # test detection limit (epg)
mix_effs <- c(0,0.5,1)   # fecal sample mixing efficiency
resamp <- c(T,F)         # with resampling Y/N   

# scenarios, each one has a different combination of m, the mean parasite
# burden, and k, the inverse aggregation parameter
scenarios <- expand.grid(m = ms, k = ks, resamp = resamp, 
                         ss_wt = subsamp_wts, ss_sd = subsamp_sds, 
                         mef = mix_effs, dl = dls
                         )

out <- lapply(1:nrow(scenarios), \(i) {
  with(scenarios[i,],{
    db <- replicate(niter,{
      # outdb <- data.frame(m = numeric(nrow(scenarios)), k = 0, sw = 0, cv = 0, dl = 0, n = 0, epg_obs = 0, epg_tru = 0)
      # random coef of var from triangular dist, random target weight between 1
      # and 6
      n_samp <- sample(seq(10,Na,2),1)
      wgt_cv <- gen_cv(1)
      #wgt_sd <- s_targ_wt*wgt_cv
      # generate EPG values for the population, one per individual. This
      # represents the mean FEC per individual
      feces_epg <- rnbinom(n = Na, size = k, mu = m)
      
      # simulate samples of these following a Poisson dist, depending on random n
      # samples per individual
      inds_sampled <- sample(1:Na, n_samp, replace = resamp)
      n_inds_sampled <- length(unique(inds_sampled))
      
      # simulate variable weight taken from the animals
      # using a normal distribution in some cases generates negative weights, so
      # it might be better to use a log-normal distribution instead
      sdlog <- sqrt(log(wgt_cv^2+1))
      sample_weights <- rlnorm(n_samp, meanlog = log(targ_wt), sdlog = sdlog)
      
      # total eggs in samples
      sample_eggs <- rpois(n_samp, sample_weights*feces_epg[inds_sampled])
      
      # sample epg
      sample_epg <- sample_eggs/sample_weights
      
      # subsample in lab
      # 1. individual
      # take subsample with the weight for flotation, with some error
      ind_subsamp_wts <- rnorm(n_samp,flot_wt,ss_sd)
      # calculate number of eggs in subsample, and corresponding epg
      ind_subsamp_eggs <- rpois(n_samp, flot_wt*sample_epg)
      ind_epg_tru <- ind_subsamp_eggs/ind_subsamp_wts
      # draw random number of how many are observed, this depends on the
      # detection limit
      ind_eggs_obs <- rpois(n_samp, ind_epg_tru/dl)
      # calculate epg
      ind_epg <- ind_eggs_obs*dl
      
      # 2. pooled 
      # subsample weights, with some tolerance
      comp_subsamp_wts <- rnorm(n_samp,ss_wt,ss_sd)
      # number of eggs in each subsample, drawn from Poisson dist, and
      # corresonding epg
      comp_subsamp_eggs <- rpois(n_samp, comp_subsamp_wts*sample_epg)
      comp_subsamp_epg <- comp_subsamp_eggs/comp_subsamp_wts
      # mixing with some efficiency mef redistributes the eggs. Poor mixing
      # (mef~0) would result in taking a sample from only one individual. We
      # simulate this using a weighted mean.
      samples_needed <- flot_wt/ss_wt # how many individuals needed to fill the weight
      wm_wts <- c(rep(1,samples_needed),numeric(n_samp-samples_needed)) #
      # redistribute the weights given some mixing efficiency
      wm_wts_mxd <- wm_wts*(1-mef)+mean(wm_wts)*mef
      comp_flot_eggs <- sum(floor(wm_wts_mxd*comp_subsamp_eggs))
      comp_flot_epg_tru <- sum(comp_flot_eggs)/flot_wt
      comp_egg_obs <- rpois(1, comp_flot_epg_tru/dl)
      comp_epg <- comp_egg_obs*dl
      
      # output
      c(i,n_samp, resamp, n_inds_sampled, mef, comp_epg, mean(ind_epg), var(ind_epg), mean(feces_epg), var(feces_epg))
    }
    )
    return(t(db))
  }
  )
}
)
#### Post processing ####
# join resulting dbs
outdb <- do.call(rbind, out)
outdb <- as.data.frame(outdb)
names(outdb) <- c("scenario","n_samp", "resamp", "n_inds", "mef", "epg_comp", "epg_ind_avg", "epg_ind_var", "epg_tru_avg", "epg_tru_sd")
outdb <- left_join(outdb, mutate(scenarios,scenario = row_number()))
head(outdb)

# export
saveRDS(outdb, "sim_output.rds")

