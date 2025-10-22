library(ggplot2)
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
    feces_epg <- rnbinom(n = Na, size = k, mu = m)
    
    # simulate samples of these following a Poisson dist, depending on random n
    # samples per individual
    sample_epg <- matrix(rpois(Na*s_targ_wt,lambda = feces_epg), nrow = Na) |> apply(1, mean)
    # sample_epg <- rpois(Na, feces_epg*sample_weights)
    # simulate variable weight taken from the sample
    # sd <- sqrt(log(wgt_cv^2+1))
    sample_weights <- rnorm(Na, s_targ_wt, wgt_sd)
    
    # total eggs in sample
    sample_eggs <- sample_epg*sample_weights
    
    # cum total weight
    cum_samp_weight <- cumsum(sample_weights)
    
    # cum egg number
    cum_egg_number <- cumsum(sample_eggs)
    
    # composite EPG
    composite_epg <- cum_egg_number/cum_samp_weight
    
    # eggs counted
    eggs_counted <- floor(composite_epg/dl)
    eggs_count_pois <- rpois(Na, lambda = eggs_counted)
    
    # observed FEC
    epg_obs <- (eggs_count_pois*dl)[ind_samp]
    
    # output
    c(m, k, s_targ_wt, wgt_sd, dl, ind_samp, epg_obs, mean(sample_epg[1:ind_samp]), mean(feces_epg))
  }
  )
  return(t(db))
}
)

# join resulting dbs
outdb <- do.call(rbind, out)
outdb <- as.data.frame(outdb)
names(outdb) <- c("m","k", "tgt_wt", "wt_sd", "dlim", "inds", "epg_comp", "epg_ind", "epg_tru")
head(outdb)
# fig 2 shows the absolute difference of the composite FEC with respect to the population m
p_absdif <- ggplot(outdb, aes(inds,abs(epg_comp-m)+0.1))+
  # geom_point()+
  facet_wrap(~m, labeller = labeller(m = \(x) paste("m =",x)))+
  geom_smooth(method = 'lm', aes(lty = as.factor(k)), color = 'black')+
  labs(x = "Individuals sampled", y = expression(paste(Log[10], " absolute difference")), lty = "k")+
  scale_y_log10()+
  theme_pubr(legend = 'right')

# relative difference plot
p_reldif <- ggplot(outdb, aes(inds,(epg_comp-m)/m))+
  geom_hline(yintercept = 0, lty = 2)+
  # geom_point()+
  facet_wrap(~m, labeller = labeller(m = \(x) paste("m =",x)))+
  geom_smooth(method = 'lm', aes(lty = as.factor(k)), color = 'black')+
  labs(x = "Individuals sampled", y = "Relative difference", lty = "k")+
  theme_pubr(legend = 'right')

# combine the figures
ggarrange(p_absdif,p_reldif, nrow = 2, common.legend = T, legend = 'right')

