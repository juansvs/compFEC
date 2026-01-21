# Pooled Fecal Egg Count simulation
# Juan S. Vargas
# November 2025
#

#' Code to simulate collection and analysis of fecal samples from goats. The
#' simulation focuses on understanding the differences between analyzing samples
#' individually and as a pooled sample
library(tidyverse)

set.seed(345)
rm(list = ls())
gen_cv <- function(n) rlnorm(n, log(0.3431), 0.3184)

niter <- 10              # number of iterations per scenario
Na <- c(10, 20, 30)     # number of animals sampled
ms <- c(100, 500, 2000)  # population mean EPG values
ks <- c(0.1, 0.5, 2, 10) # population EPG aggregation values
targ_wt <- 5             # sample (farm level) target weight
wgt_cv <- c(0.1, 0.4, 0.7) # coefficient of variation in sample weights
subsamp_wts <- c(0.5, 1, 2)# subsample weights for composite (at the lab level)
subsamp_sds <- 1/1.96 * c(0.1, 0.2, 0.5) # subsample weight sds: tolerance in g divided by 1.96 (95%CI of normal dist)
flot_wt <- 2             # subsample weight for FEC. e.g. 2 g for McMaster method
# flot_vol <- 30           # flotation solution volume (mL)
# slide_vol <- 0.3         # volume in counting slide (mL)
dilution_factors <- c(1/50, 1/20, 1/5) # egg counting test dilution factor
mix_effs <- c(0, 0.5, 1)   # fecal sample mixing efficiency
resamp <- c(TRUE, FALSE)         # with resampling Y/N   
count_duplicate <- c(TRUE, FALSE)

# scenarios, each one has a different combination of parameter values
scenarios <- expand.grid(n_samp = Na, m = ms, k = ks, targ_wt = targ_wt,
                         resamp = resamp, ss_wt = subsamp_wts, ss_sd = subsamp_sds,
                         mef = mix_effs, dl = dilution_factors, wgt_cv = wgt_cv,
                         dup = count_duplicate)

# main simulation function. Do the sampling and FEC for every parameter
# combination 10 times over and put out a database with the mean and median
# result for every one.
sampled_scenarios <- slice_sample(scenarios, n = 5000)
out <- lapply(1:nrow(sampled_scenarios), \(i) {
  with(sampled_scenarios[i,], {
    # db <- replicate(niter,{
      # generate EPG values for the population, one per individual. This
      # represents the true mean FEC per individual
      feces_epg <- rnbinom(n = 30, size = k, mu = m)
      
      # simulate samples of these following a Poisson dist, depending on random n
      # samples per individual. First get the set of individuals sampled.
      inds_sampled <- sample(1:30, n_samp, replace = resamp)
      n_inds_sampled <- length(unique(inds_sampled))
      
      # simulate variable weight taken from the animals
      # using a normal distribution in some cases generates negative weights, so
      # it might be better to use a log-normal distribution instead
      sdlog <- sqrt(log(wgt_cv^2 + 1))
      sample_weights <- rlnorm(n_samp, meanlog = log(targ_wt), sdlog = sdlog)
      
      # total eggs in samples
      sample_eggs <- rpois(n_samp, sample_weights * feces_epg[inds_sampled])
      
      # sample epg
      sample_epg <- sample_eggs / sample_weights
      
      # subsample in lab
      # 1. individual
      # take subsample with the weight for flotation, with some error
      ind_subsamp_wts <- rnorm(n_samp, flot_wt, ss_sd)
      # calculate number of eggs in subsample, and corresponding epg
      ind_subsamp_eggs <- rpois(n_samp, ind_subsamp_wts * sample_epg)
      ind_epg_tru <- ind_subsamp_eggs / ind_subsamp_wts
      # draw random number of how many are observed, this depends on the test
      # and its detection limit
      ind_eggs_obs <- rpois(n_samp, ind_epg_tru * flot_wt * dl)
      # calculate epg
      ind_epg <- ind_eggs_obs / (dl * flot_wt)
      
      # 2. pooled 
      # subsample weights, with some tolerance
      comp_subsamp_wts <- rnorm(n_samp, ss_wt, ss_sd)
      # replace negative weights with the target weight - sd
      comp_subsamp_wts[comp_subsamp_wts < 0] <- ss_wt - ss_sd
      # number of eggs in each subsample, drawn from Poisson dist
      comp_subsamp_eggs <- rpois(n_samp, comp_subsamp_wts * sample_epg)
      # mixing with some efficiency mef redistributes the eggs. Poor mixing
      # (mef~0) would result in taking a sample from only one individual. We
      # simulate this using a weighted mean.
      samples_needed <- flot_wt / ss_wt # number of individuals needed to fill the weight
      # proportion of each subsample effectively taken under the worst possible
      # mixing.
      wm_wts <- c(rep(1, samples_needed),
                  numeric(n_samp-samples_needed))
      # redistribute the weights given some mixing efficiency
      wm_wts_mxd <- wm_wts * (1 - mef) + mean(wm_wts) * mef
      # number of eggs drawn in flotation subsample. The weights multiplied by
      # the number of eggs from each subsample. These egg numbers are reshuffled
      # to simulate duplicate counting
      comp_subsamp_eggs_reordered <-
        matrix(comp_subsamp_eggs[replicate(3, sample(n_samp))], ncol = 3)
      comp_flot_eggs <- apply(wm_wts_mxd * comp_subsamp_eggs_reordered, 2, \(x) sum(floor(x)))
      # number of eggs observed
      comp_egg_obs <- rpois(3, comp_flot_eggs * dl)
      comp_epg <- comp_egg_obs / (flot_wt * dl)
      # count multiple slides
      if(dup) {
        # check that slides 1 and 2 are within 10% of each other,
        # if not read a third slide
        if(all(comp_epg[1:2] > 0)) {
          comp_count_dif <- abs(comp_epg[1] - comp_epg[2]) / comp_epg[2]
          if(comp_count_dif <= 0.1) {
            comp_epg <- mean(comp_epg[1:2])
          } else {
            comp_epg <- mean(comp_epg)
          }
        } else {
          comp_epg <- mean(comp_epg)
        }
      } else {
        comp_epg <- comp_epg[1]
      }
      # output
      c(i, n_inds_sampled, comp_epg,
        mean(ind_epg), sd(ind_epg))
    }
    )
}
)

### Post processing ###
# join resulting dbs
outdb <- as.data.frame(do.call(rbind, out))
outdb <- outdb[complete.cases(outdb),]
names(outdb) <- c("scenario", "n_inds", "epg_comp",
                  "epg_ind_avg", "epg_ind_sd")
outdb <- left_join(outdb, mutate(sampled_scenarios, scenario = row_number()))
head(outdb)

# export
saveRDS(outdb, "sim_output.rds")

