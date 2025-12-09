library(tidyverse)
library(performance)
library(glmmTMB)
library(mcr)
outdb <- readRDS('sim_output.rds')

## Relative error
outdb_lg <- mutate(outdb, scenid = row_number()) %>% 
  pivot_longer(c(epg_comp, epg_ind_avg), names_to = "method",
               values_to = "epg") %>%
  mutate(err = epg-m, abserr = abs(err), relerr = (epg-m)/m,
         absrelerr = abs(relerr))
# quantiles for rel error with diff methods
split(outdb_lg, ~method) %>% lapply("[[", 'relerr') %>% lapply(quantile)
# quantiles of rel error separated additionally by m and k
split(outdb_lg, ~m+k+method) |> lapply("[[","relerr") |> lapply(summary)
# absolute relative error within 10% of real value
split(outdb_lg, ~m+k+method) |> lapply("[[","relerr") |>
  lapply(abs) |> lapply("<",0.1) |> sapply(mean)
filter(outdb_lg, method == 'epg_comp', k<=0.5)|> split(~mef+m) |>
  lapply("[[", "relerr") |> lapply(summary)
# proportion of false negatives by method
outdb_lg %>% mutate(falseneg = relerr == -1) %>%
  summarise(mean(falseneg), .by = c(method))
# 6.64% false negatives with compound method, 0.06% with ind. method.
outdb_lg %>% mutate(falseneg = relerr == -1) %>%
  summarise(fn = mean(falseneg), .by = c(method, m, k))
# with high aggregation false neg rate is even higher: 23.4% (m = 100, k = 0.1),
# compared to 0.3% with ind method. Also substantial at m = 500 and k = 0.1
# (13.1% vs 0.2%), and at m = 2000 (10.0% vs. 0%). 

## composite, influence of mixing efficiency
mutate(outdb, absrelerr = abs(epg_comp-m)/m) %>% split(~mef) %>%
  lapply("[[", "absrelerr") %>% lapply(summary)
filter(outdb, k<=0.5) %>% mutate(relerr = (epg_comp-m)/m) %>% split(~mef) %>%
  lapply("[[", "relerr") %>% sapply(quantile)

# misclassification
filter(outdb, m == 100, mef==1) %>% summarise(hic = mean(epg_comp>=500, na.rm=T), hii = mean(epg_ind_avg>=500, na.rm=T), .by = k)

filter(outdb, m == 500, mef ==1) %>% summarise(hic = mean(epg_comp>=2000, na.rm=T), hii = mean(epg_ind_avg>=2000, na.rm=T),
                                               loc = mean(epg_comp<=100, na.rm=T), loi = mean(epg_ind_avg<=100, na.rm=T),
                                               .by = k)
filter(outdb, m == 2000, mef ==1) %>% summarise(loc = mean(epg_comp<=500, na.rm=T), loi = mean(epg_ind_avg<=500, na.rm=T),
                                                .by = k)
#### Statistical tests ####
# false negatives
fn_db <- mutate(outdb_lg, fn = relerr == -1) %>% slice_sample(prop = 0.1)
fn_glmm <- glmmTMB(fn ~ n_samp + resamp + wt_cv + factor(mef) * method +
                     factor(m) + factor(k) + ss_sd + ss_wt + factor(dl) + (1|scenid),
                   family = 'binomial', data = fn_db)
summary(fn_glmm)

# same analysis but include only one set of m and k, and focus on factors that
# can be controlled
fn_db2 <- mutate(outdb_lg, fn = relerr == -1) %>%
  filter(m == 500, k == 0.5) %>%
  slice_sample(prop = 0.2)
fn_glmm2 <- glmmTMB(fn ~ n_samp + resamp + wt_cv + factor(mef) * method +
                     ss_sd + ss_wt + factor(dl) + (1|scenid),
                   family = 'binomial', data = fn_db,
                   control = glmmTMBControl(optimizer = optim,
                                            optArgs = list(method = "BFGS")))
summary(fn_glmm2)
# these models suggest that the most relevant factors for the false negative
# rate are the method, mixing efficacy, and aggregation level. The method has
# the greatest effect by far, using independent samples practically eliminates
# the probability of getting a false negative. Mixing well (mef = 1) reduces it
# by 96.2%, and even partial mixing reduces it by 94.2%. Among the factors that
# are beyond the samplers control, aggregation seems to be the most relevant. At
# k = 0.5 the probability is decreased by 88.7% compared to the highest
# aggregation, and at k = 2 the reduction is of 98.3%.

glmm_db <- slice_sample(outdb_lg, prop = 0.1) %>% mutate(resp = abs(0.1+epg-m)/m) 
absrelerr_glmm <- glmmTMB(resp~n_samp*resamp+wt_cv+factor(mef)*method+factor(m)+k+ss_sd+ss_wt+(1|scenid), 
                     family = Gamma(link = 'log'), data = glmm_db)
summary(absrelerr_glmm)
exp(fixef(absrelerr_glmm)$cond)-1

## COmparison of estimates ##
mutate(outdb, reldif = (epg_comp-epg_ind_avg)/epg_ind_avg) %>% 
  split(~m+k) %>% 
  lapply(pull, "reldif") %>% 
  sapply(quantile, na.rm=T)
  

#### Passing-Bablok regression for method comparison ####
mcsplit <- split(outdb, ~m+k) |> lapply(\(x) mcreg(x$epg_comp, x$epg_ind_avg, method.reg = 'PaBaLarge'))
mc1sim <- mcreg(outdb$epg_comp, outdb$epg_ind_avg, 
             method.reg = "PaBaLarge", na.rm = T)
range(nbfits_joint$epg_comp, na.rm = T) # 0 - 2650
mc1db <- tibble(epg_comp = 0:2650) %>% 
  mutate(est = mc1@para[1]+epg_comp*mc1@para[2],
         lc = mc1@para[5]+epg_comp*mc1@para[6],
         uc = mc1@para[7]+epg_comp*mc1@para[8])