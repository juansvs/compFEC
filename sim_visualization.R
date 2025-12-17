library(tidyverse)
library(ggplot2)
library(ggpubr)

# import data
outdb <- readRDS("sim_output.rds")

# density plot of relative error wrt m, k, and method
p_difdens <- mutate(outdb,
                    comp = ((epg_comp - m) / m),
                    ind = ((epg_ind_avg - m) / m)) %>%
  filter(dl == 1/50) %>%
  pivot_longer(cols = c(comp, ind),
               names_to = "method",
               values_to = "rdif") %>%
  ggplot(aes(rdif)) +
  geom_density(aes(y = after_stat(scaled),
                   fill = method),
               show.legend = F, lwd = 0.5) +
  geom_vline(xintercept = 0, lty = 2) +
  facet_grid(m ~ k, scales = "free",
             labeller = labeller(m = \(x) paste("m =", x), 
                                 k = \(x) paste("k =", x))) +
  labs(x = "Relative error", y = "Density", fill = "Method")

# Export figure
p_difdens + theme_classic() + scale_fill_manual(values = c("gray80", "gray20"))
png("figures/relerrdensplot.png", 6, 4, units = 'in', res = 300)
p_difdens +
  theme_classic() +
  coord_cartesian(xlim = c(NA, 3)) +
  scale_fill_discrete(palette = hcl.colors(2, alpha = 0.5))
dev.off()

# dens plot of rel. error wrt mixing efficiency for diff. mean burdens, high
# aggregation
p_mef_dens <- filter(outdb, resamp == 0, k <= 0.5, dl == 1/50) %>%
  ggplot(aes((epg_comp - m) / m)) +
  geom_density(aes(fill = factor(mef)), show.legend = F) +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap( ~ m,
              labeller = labeller(m = \(x) paste("m =", x)),
              scales = 'free') +
  labs(x = "Relative error", y = "Density") +
  coord_cartesian(xlim = c(NA, 3)) +
  scale_fill_discrete(palette = hcl.colors(3, "PuRd", alpha = 0.5))
# export
png("figures/mefdensplot.png", 5, 2, 'in', res = 300)
p_mef_dens + theme_classic()
dev.off()


# fig 2 shows the absolute difference of the composite FEC with respect to the
# population m
p_absdif <- ggplot(outdb, aes(n_samp,abs(epg_comp-m)))+
  # geom_point()+
  facet_wrap(~m, labeller = labeller(m = \(x) paste("m =",x)))+
  geom_point()+
  labs(x = "Number of samples", y = "Absolute difference")
p_absdif + theme_pubr(legend = 'right')

p_absdif_lm <- ggplot(outdb, aes(n_samp,abs(epg_comp-m)))+
  # geom_point()+
  facet_wrap(~m, labeller = labeller(m = \(x) paste("m =",x)))+
  geom_smooth(method = 'lm', aes(lty = as.factor(k)), color = 'black')+
  labs(x = "Number of samples", y = "Absolute difference", lty = "k")+
  scale_y_continuous()
p_absdif_lm + theme_pubr(legend = 'right')

# relative difference dot plot
p_reldif_pts <- filter(outdb, resamp == 0, mef == 1) %>% 
  ggplot(aes(n_samp,(epg_comp-m)/m))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_jitter(alpha = 0.1)+
  facet_grid(m~k, labeller = labeller(m = \(x) paste("m =",x),
                                      k = \(x) paste("k =",x)))+
  labs(x = "Individuals sampled", y = "Relative difference")
p_reldif_pts

# relative difference lm plot
p_reldif1 <- filter(outdb, mef == 1) %>% 
  ggplot(aes(n_samp,(epg_comp-m)/m))+
  geom_hline(yintercept = 0, lwd = 0.5)+
  facet_wrap(~m, labeller = labeller(m = \(x) paste("m =",x)))+
  geom_smooth(method = 'lm',aes(color = as.factor(k), lty = factor(resamp)))+
  labs(x = "Individuals sampled", y = "Relative error", lty = "Resampling", color = "k")
  theme_pubr(legend = 'right')
  scale_y_continuous(transform = 'log1p')
p_reldif1
# same for indiv
p_reldif2 <-  ggplot(outdb, aes(n_samp,(epg_ind_avg-m)/m))+
  geom_hline(yintercept = 0, lwd = 0.5)+
  facet_wrap(~m, labeller = labeller(m = \(x) paste("m =",x)))+
  geom_smooth(method = 'lm',aes(color = as.factor(k), lty = factor(resamp)))+
  labs(x = "Individuals sampled", y = "Relative error", lty = "Resampling", color = "k")
  theme_pubr(legend = 'right')
  scale_y_continuous(transform = 'log1p')
p_reldif2
ggarrange(p_reldif1, p_reldif2, labels = 'auto', legend.grob = get_legend(p_reldif1), nrow = 2, legend = 'right')

# pointrange vs number of samples
p_var <-filter(outdb, resamp == 0, mef ==1) %>% 
  ggplot(aes(n_samp,(epg_comp-m)/m))+
  geom_hline(yintercept = 0, lty = 2)+
  stat_summary(aes(color = factor(k)), 
               fun.data = \(x) data.frame(ymin = quantile(x, p = 0.25), 
                                          y = median(x), 
                                          ymax = quantile(x, p = 0.75)), 
               geom = "pointrange", position = position_dodge(width = 1)) +
  facet_grid(~m, labeller = labeller(m = \(x) paste("m =", x))) +
  labs(x = "Individuals sampled", y = "Relative difference", color = "k")
p_var

# Variability (IQR) vs number of samples
ggplot(outdb, aes(inds, abs(epg_comp - m) / m))+
  stat_summary(fun = \(x) sd(x) / mean(x)) +
  geom_point(alpha = 0.2) +
  facet_grid(m ~ k)

# Effect of subsample precision
filter(outdb, resamp == 0, mef ==1) %>%
  ggplot(aes(ss_sd * 1.96, group = factor(ss_sd), log1p((epg_comp - m) / m))) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot() +
  # scale_y_continuous(transform = 'log10')+
  facet_grid( ~ k) +
  labs(x = "Subsample tolerance (g)", y = "Relative error")
  

# effect of variance at farm levele
filter(outdb, resamp == 0, mef == 1) %>% 
  ggplot(aes(wt_cv, log1p(abs(epg_comp - m) / m))) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() +
  facet_grid(m ~ k) +
  labs(x = "Sample weight CV", y = "Relative error")

# comp vs ind
pcomp <- lapply(c(100, 500, 2000), \(x) {
  ggplot(filter(outdb, m == x), aes(epg_comp, epg_ind_avg)) +
  geom_abline(slope = 1, linetype = 2) +
  geom_hline(aes(yintercept = m), lty = 3) +
  geom_vline(aes(xintercept = m), lty = 3) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = 'lm', formula = y ~ x - 1) +
  facet_grid(~k, scales = "free",  labeller = labeller(m = \(x) paste("m =",x), 
                                                        k = \(x) paste("k =",x)))+
  # scale_y_log10()+scale_x_log10()+
  labs(x = "Composite FEC", y = "Individual FEC") +
  theme_classic()})
png("figures/compscatter.png", 6, 5,'in', res = 300)
ggarrange(pcomp[[1]], pcomp[[2]], pcomp[[3]],
          labels = 'auto', nrow = 3)
dev.off()
# comp vs ind (ratios)
p_comp <- ggplot(outdb, aes(epg_comp / m, epg_ind_avg / m)) +
  geom_abline(slope = 1, linetype = 2)+
  geom_point(alpha = 0.2) +
  geom_smooth(method = 'lm', se = F) +
  facet_grid(m~k, scales = "free")+
  # scale_y_log10()+scale_x_log10()+
  labs(x = "Composite FEC ratio", y = "Individual FEC ratio")
p_comp
# compare errors
p_err_comp <- mutate(outdb, comp = (epg_comp - m) / m,
                     ind = (epg_ind_avg - m) / m) %>%
  pivot_longer(c(comp, ind), 
                           names_to = "method",values_to = "epg") %>%
  ggplot(aes(method,epg)) +
  geom_violin(aes(fill = method), show.legend = F) +
  stat_summary(geom = 'pointrange',
               fun.data = \(x) data.frame(ymin = quantile(x, p = 0.25),
                                          y = median(x),
                                          ymax = quantile(x, p = 0.75))) +
  geom_hline(yintercept = 0,lty = 2)+
  facet_grid(m~k, labeller = labeller(m = \(x) paste("m =", x), 
                                      k = \(x) paste("k =", x))) +
  theme_classic(base_size = 12)+
  labs(x = "Method", y = "Relative error")+
  scale_x_discrete(labels = c("Pooled", "Indiv."))
p_err_comp

## density plot of estimated values
pivot_longer(outdb, cols = c(epg_comp, epg_ind_avg), names_to = "method",
             values_to = "est") %>%
  ggplot(aes(est)) +
  geom_vline(aes(xintercept = m), lty = 2)+
  geom_density(aes(y = after_stat(scaled),fill = method),
               alpha = 0.4, show.legend = F) +
  facet_grid(m~k, scales = "free",
             labeller = labeller(m = \(x) paste("m =", x), 
                                 k = \(x) paste("k =", x))) +
  labs(x = "FEC estimate (epg)", y = "Density", fill = "Method")





# scatter plot comp vs ind
filter(outdb, resamp == 0) %>% 
  ggplot(aes(epg_comp,y = epg_ind_avg)) +
  geom_abline(slope = 1, lty = 2) +
  geom_point(alpha = 0.3) +
  facet_grid(m ~ k, scales = 'free')

## Bland-Altman plot
mutate(outdb, FECmean = (epg_comp + epg_ind_avg) / 2,
       FECdif = epg_comp - epg_ind_avg) %>%
  ggplot(aes(FECmean, FECdif)) +
  geom_vline(aes(xintercept = m)) +
  geom_point()+
  facet_grid(m ~ k, scales = "free")

## relative difference between methods
ggplot(outdb, aes((epg_comp - epg_ind_avg) / epg_ind_avg)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_density(aes(fill = factor(m)), alpha = 0.5, show.legend = F) +
  scale_x_continuous(limits = c(NA, 5)) +
  scale_fill_discrete(palette = hcl.colors(3, "PuRd")) +
  labs(x = "Relative difference (comp. - ind.)", y = "Density") +
  facet_wrap( ~ k, labeller = \(x) label_both(x, sep = " = ")) +
  theme_classic()
  
# effect of mixing efficiency
filter(outdb, resamp == 0, dl == 1/50) %>%
  ggplot(aes((epg_comp - m) / m)) +
  geom_density(aes(y = after_stat(scaled),
                   fill = factor(mef))) +
  geom_vline(xintercept = 0, lty = 2) +
  # scale_x_continuous(transform = 'log10')+
  facet_grid(m~k,
             labeller = labeller(m = \(x) paste("m =", x),
                                      k = \(x) paste("k =", x))) +
  labs(x = "Relative error", y = "Density", fill = "Mixing efficiency") +
  coord_cartesian(xlim = c(NA, 4)) +
  scale_fill_discrete(palette = hcl.colors(3, alpha = 0.5))
  
## mixing efficiency violin plot
filter(outdb, resamp == 0, k <= 0.5) %>%
  ggplot(aes(factor(mef), (epg_comp - m) / m)) +
  geom_violin(aes(fill = factor(mef)), show.legend = F) +
  geom_hline(yintercept = 0, lty = 2) +
  # scale_x_continuous(transform = 'log10')+
  facet_wrap(~m,
             labeller = labeller(m = \(x) paste("m =", x)),
             scales = 'free') +
  labs(x = "Mixing efficiency f", y = "Density")+
  theme_classic()



filter(outdb, resamp == 0) %>% 
  ggplot(aes(factor(mef), (epg_comp - m) / m)) +
  geom_boxplot() +
  facet_grid(m ~ k, labeller = labeller(m = \(x) paste("m =", x), 
                                      k = \(x) paste("k =", x))) +
  labs(x = "Mixing efficiency", y = "Relative error")

## Absolute error ##
pivot_longer(outdb, cols = c("epg_comp", "epg_ind_avg"), 
             names_to = "method", values_to = "epg") %>%
  mutate(abserr = abs(epg - m)) %>%
  ggplot(aes(abserr)) +
  geom_density(aes(fill = method), alpha = 0.5, show.legend = F) +
  # scale_x_continuous(transform = 'log1p')+
  labs(x = "Absolute error") +
  facet_grid(m ~ k, scales='free',
             labeller = labeller(m = \(x) paste("m =", x),
                                 k = \(x) paste("k =", x))) +
  theme_pubr()
## Absolute relative error ##
pivot_longer(outdb, cols = c("epg_comp", "epg_ind_avg"), 
             names_to = "method", values_to = "epg") %>%
  mutate(absrelerr = (abs(epg - m) / m)) %>%
  ggplot(aes(absrelerr)) +
  geom_density(aes(fill = method), alpha = 0.5, show.legend = F) +
  # geom_vline(aes(xintercept = med, color = method, linetype = method), data = summarise(outdb_lg, med=median(absrelerr, na.rm=T), .by = c(m,k,method)), show.legend = F)+
  # scale_x_continuous(transform = 'log10')+
  labs(x = "Absolute relative error") +
  facet_grid(m ~ k, scales = 'free',
             labeller = labeller(m = \(x) paste("m =", x), 
                                 k = \(x) paste("k =", x))) +
  theme_pubr()

pivot_longer(outdb, cols = c("epg_comp", "epg_ind_avg"), 
             names_to = "method", values_to = "epg") %>%
  mutate(absrelerr = abs(epg - m) / m) %>%
  ggplot() +
  geom_density(aes(n_samp, y = absrelerr)) +
  labs(x = "Number of samples") +
  facet_grid(m ~ k, scales = 'free',
             labeller = labeller(m = \(x) paste("m =", x), 
                                 k = \(x) paste("k =", x))) +
  theme_pubr()

### False negative rate plot, method as color, wrt m/k
pivot_longer(outdb, cols = c("epg_comp", "epg_ind_avg"), 
             names_to = "method", values_to = "epg") %>% 
  mutate(fn = epg == 0) %>% 
  ggplot(aes(method, fn))
  
pdf(width = 6, height = 4)
p_absdif
p_reldif
p_reldif_pts
p_var
p_comp
p_err_comp
p_difdens
dev.off()
