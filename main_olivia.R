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

gen_cv <- function(n) VGAM::rtriangle(n, lower = 0.2, theta = 0.3, upper = 0.73)# coefficient of variation in weight samples, random value from a triangular dist. (0.2, 0.3, 0.73)

niter <- 2000 # number of iterations per scenario
ms <- c(100, 500, 2000) # population mean EPG values
ks <- c(0.1, 0.5, 2, 10) # population aggregation values
Na <- 30 # number of animals

# scenarios, each one has a different combination of m, the mean parasite
# burden, and k, the inverse aggregation parameter
scenarios <- expand.grid(m = ms, k = ks)

out <- lapply(1:nrow(scenarios), \(i) {
  m <- scenarios$m[i]
  k <- scenarios$k[i]
  db <- replicate(niter,{
    # outdb <- data.frame(m = numeric(nrow(scenarios)), k = 0, sw = 0, cv = 0, dl = 0, n = 0, epg_obs = 0, epg_tru = 0)
    # random coef of var from triangular dist, random target weight between 1
    # and 6
    s_targ_wt <- sample(1:6,1)
    dl <- sample(1:50, 1)
    ind_samp <- sample(1:Na,1)
    wgt_cv <- gen_cv(1)
    wgt_sd <- s_targ_wt*wgt_cv
    # generate EPG values for the population, one per individual
    feces_epg <- rnbinom(n = ind_samp, size = k, mu = m)
    
    # simulate samples of these following a Poisson dist, depending on random n
    # samples per individual
    # sample_epg <- matrix(rpois(Na*s_targ_wt,lambda = feces_epg), nrow = Na) |> apply(1, mean)

    # simulate variable weight taken from the sample
    # sd <- sqrt(log(wgt_cv^2+1))
    sample_weights <- rnorm(ind_samp, s_targ_wt, wgt_sd)
    
    # total eggs in sample
    # sample_eggs <- sample_epg*sample_weights
    sample_eggs <- rpois(ind_samp, sample_weights*feces_epg)
    
    # cum total weight
    cum_samp_weight <- sum(sample_weights)
    
    # cum egg number
    cum_egg_number <- sum(sample_eggs)
    
    # composite EPG
    composite_epg <- cum_egg_number/cum_samp_weight
    
    # eggs counted
    exp_egg_count <- floor(composite_epg/dl)
    eggs_count_pois <- rpois(1, lambda = eggs_counted)
    
    # observed FEC
    epg_obs <- (eggs_count_pois*dl)[ind_samp]
    
    # output
    c(m, k, s_targ_wt, wgt_sd, dl, ind_samp, epg_obs, mean(sample_epg[1:ind_samp]), mean(feces_epg))
  }
  )
  return(t(db))
}
)
#### Post processing ####
# join resulting dbs
outdb <- do.call(rbind, out)
outdb <- as.data.frame(outdb)
names(outdb) <- c("m","k", "tgt_wt", "wt_sd", "dlim", "inds", "epg_comp", "epg_ind", "epg_tru")
head(outdb)

# export
saveRDS(outdb, "sim_output.rds")

