library(mcr)
library(ggpubr)
library(MASS)
library(tidyverse)

#### Load data ####
## Summarised database with mean and CME k per farm
# farm_res <- read.csv('../farm_results.csv')
# head(farm_res)
# ggplot(farm_res, aes(fec_mean, fec_comp))+
#   geom_abline(slope = 1, lty = 2)+
#   geom_smooth(method = 'lm', color = 'black', lwd=0.5) +
#   geom_point(size = 2)+
#   labs(x = "Mean FEC (epg)", y = "Composite FEC (epg)") +
#   theme_pubr()

## Raw epg data #
epg_raw_data <- read.csv('raw_epg.csv')

#### Calculate CME ####
cme <- summarise(epg_raw_data, m = mean(epg, na.rm =T),
                 v = var(epg, na.rm = T), n = n(), .by = id) %>%
  mutate(cme = (m^2 - (v / n)) / (v - m))

####----- Fit distributions ----
epgdata <- split(epg_raw_data, ~id)|>
  lapply('[[', "epg") |>
  lapply(\(x) subset(x, subset = complete.cases(x)))
# Negative binomial
nbfits <- lapply(epgdata, 
  \(x) if(sum(x)>0) fitdistr(x, densfun = 'negative binomial') else NULL) |>
  lapply(coef)

nbfitsdb <- do.call(rbind, nbfits) %>%
  as.data.frame() %>%
  rownames_to_column("id")

# Poisson
poisfits <- lapply(epgdata,
  \(x) if(sum(x)>0) fitdistr(x, densfun = 'Poisson') else NULL) |>
  lapply(coef)
  
#### Goodness of fit ####

# Create a dataframe with the goodness of fit metric for the Poisson and 
# negative binomial distribution.
gof_df <- data.frame()
for(i in 1:length(nbfits)) {
  if(!is.null(nbfits[[i]])) {
    size <- nbfits[[i]][1]
    mu <- nbfits[[i]][2]
    lambda <- poisfits[[i]]
    ksnb <- ks.test(epgdata[[i]], pnbinom, size = size, mu = mu)
    kspoi <- ks.test(epgdata[[i]], ppois, lambda = lambda)
    gofi <- data.frame(id = names(nbfits)[i], 
                       nbstat = ksnb$statistic, nbp = ksnb$p.value,
                       poistat = kspoi$statistic, poip = kspoi$p.value)
    gof_df <- rbind(gof_df, gofi)
  }
}

gof_df
# See how many fit each distribution
sum(gof_df$nbp>0.05) # 18 of 26 adjust to the neg bin
sum(gof_df$poip>0.05) # none adjust to a Poisson

#### Composite epg ####
# read in data
comp_epg <- read.csv('pooled_epg.csv') %>% 
  pivot_longer(cols = starts_with('comp'), names_to = 'test',
               values_to = 'epg_comp')
# summarise mean epg and subsample weights
comp_epg_summ <- summarise(comp_epg, epg_comp = mean(epg_comp, na.rm=T),
                           subsamp_wt = mean(subsamp_wt, na.rm=T), .by = id)
comp_epg_summ

#### Join dataframes #####
nbfits_joint <- left_join(cme, nbfitsdb) %>% 
  left_join(comp_epg_summ) %>% 
  left_join(gof_df) %>% 
  mutate(mu1q = qnbinom(0.25, size, mu = mu), mu3q = qnbinom(0.75, size, mu = mu))# quartiles of estimated nb dist
nbfits_joint

#### Stats #####

# lmer of FEC wrt covariates
fec_glmm <- lmer(m ~ sample_month + (1 | farm), nbfits_joint)
summary(fec_glmm)

# we fit multiple models and there is no clear link between FEC and any of the 
# covariates.

k_glmm <- lmer(size ~ Size + (1 | farm), nbfits_joint)
summary(k_glmm)
k_lm <- glm(size ~ Size + System, family = gaussian(link = "log"), data = nbfits_joint)
summary(k_lm)
# correlation between individual and composite mean burden estimates
plot(nbfits_joint$m, nbfits_joint$epg_comp, pty='s')
abline(0,1)
cor.test(nbfits_joint$m, nbfits_joint$epg_comp)
# there is a significant correlation between the individual mean epg
# and the composite epg estimate
plot(size ~ m, nbfits_joint)
cor.test(nbfits_joint$size, nbfits_joint$m)
# r = 0.37, t = 2.00, p = 0.0566

# t test of difference between k values
t.test(nbfits_joint$cme, nbfits_joint$size, paired = T)
# t test of difference between FEC estimates (ind vs comp)
t.test(nbfits_joint$m, nbfits_joint$epg_comp, paired = T)
# neither are significant

# Passing-Bablok regression for method comparison
mc1 <- mcreg(nbfits_joint$m, nbfits_joint$epg_comp, 
             method.reg = "PaBa", na.rm = T)
range(nbfits_joint$epg_comp, na.rm = T) # 0 - 2650
mc1db <- tibble(epg_comp = 0:2650) %>% 
  mutate(est = mc1@para[1]+epg_comp*mc1@para[2],
         lc = mc1@para[5]+epg_comp*mc1@para[6],
         uc = mc1@para[7]+epg_comp*mc1@para[8])
summary(mc1)
# repeat with log-transformed data
mc1_log <- mcreg(log10(nbfits_joint$m+1), log10(nbfits_joint$epg_comp + 1),
                 method.reg = "PaBa", na.rm = TRUE)
mc1_log_db <- tibble(epg_comp = 0:2650) %>% 
  mutate(est = mc1_log@para[1]+epg_comp*mc1_log@para[2],
         lc = mc1_log@para[5]+epg_comp*mc1_log@para[6],
         uc = mc1_log@para[7]+epg_comp*mc1_log@para[8])

mc2 <- mcreg(nbfits_joint$cme, nbfits_joint$size,
             method.reg = "PaBa", na.rm = TRUE)
mc2@para
range(nbfits_joint$cme, na.rm = TRUE)
mc2db <- tibble(epg_comp = seq(0, 6, 0.5)) %>% 
  mutate(est = mc2@para[1]+epg_comp*mc2@para[2],
         lc = mc2@para[5]+epg_comp*mc2@para[6],
         uc = mc2@para[7]+epg_comp*mc2@para[8])

#### Visualizations ####
# k estimates, MLE vs CME
p1 <- ggplot(nbfits_joint, aes(size, cme)) +
  geom_abline(slope = 1) +
  # geom_smooth(method = 'lm', color = 'gray30')+
  # geom_ribbon(aes(epg_comp, est, ymin = lc, ymax = uc),
  #             data = mc2db, alpha = 0.5, fill = "dodgerblue") +
  # geom_line(aes(epg_comp, est),data = mc2db, color = 'dodgerblue', lty = 2) +
  geom_point(aes(color = nbp > 0.05, shape = nbp > 0.05), show.legend = FALSE) +
  labs(x = "MLE k", y = "CME k")
p1 + theme_pubr()

# epg, ind vs. comp
p2 <- mutate(nbfits_joint, ind_sem = sqrt(v/n)) %>% 
  ggplot(aes(epg_comp, m)) +
  geom_abline(slope = 1) +
  geom_ribbon(aes(epg_comp, est, ymin = lc, ymax = uc),
              data = mc1db, alpha = 0.5, fill = "dodgerblue") +
  geom_line(aes(epg_comp, est), data = mc1db, color = 'dodgerblue', lty = 2) +
  geom_pointrange(aes(ymin = m - ind_sem, ymax = m + ind_sem))+
  labs(x = "Composite FEC (epg)", y = "Individual FEC (epg)")
p2 + theme_pubr() + coord_flip()

# p2 on log scale, without PB regression
p2_log <- mutate(nbfits_joint, ind_sem = sqrt(v/n)) %>% 
  ggplot(aes(m + 1, epg_comp + 1, color = nbp > 0.05, shape = nbp > 0.05)) +
  geom_abline(slope = 1) +
  # geom_ribbon(aes(epg_comp, est, ymin = lc, ymax = uc),
  #             data = mc1_log_db, alpha = 0.5, fill = "dodgerblue") +
  # geom_line(aes(epg_comp, est), data = mc1_log_db, color = 'dodgerblue', lty = 2) +
  geom_linerange(aes(xmin = m - ind_sem + 1, xmax = m + ind_sem + 1), show.legend = FALSE) +
  geom_point(show.legend = FALSE) +
  labs(x = "Individual FEC (epg)", y = "Composite FEC (epg)")
p2_log + theme_pubr() + coord_transform(x = "log10", y = "log10")

# correlation between mean and k
p3 <- pivot_longer(nbfits_joint, c(m, epg_comp),
                   names_to = "method", values_to = "est") %>%
  ggplot() + geom_point(aes(est, size, shape = method, color = method)) +
    labs(x = "FEC (epg)", y = "MLE k")
p3

# arrange the plots together
ggarrange(p2_log +
            theme_pubr() +
            scale_x_continuous(breaks = c(20, 200, 2000)) +
            scale_y_continuous(breaks = c(20, 200, 2000)) +
            coord_transform(x = "log10", y = "log10"),
          p1 + theme_pubr(),
          # p3 + theme_classic(),
          labels = "auto", ncol = 2,
          legend = "none")
