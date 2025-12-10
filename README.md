---
editor_options: 
  markdown: 
    wrap: 72
---

# Composite Faecal Egg Counts in goats

This code simulates sampling feces in livestock farms, specifically from
goats, and processing them in laboratory to quantify the number of
helminth eggs. Fecal egg counts are a commonly used technique to
estimate the parasite burden in farms. Obtaining and analyzing
individual samples often is coslty and time-consuming, both for farmers
who must collect the samples and for laboratory analysts who must
perform many egg counts. Alternatively, collecting many samples
together, mixing them in the laboratory greatly reduces the number of
samples analyzed in the lab, while still providing an estimate of the
mean burden across all individuals in the farm. Parasites are commonly
aggregated, with few individuals harboring most of the adult parasites,
and this could create large errors in the estimation of the mean burden.
To assess the magnitude and relevance of the potential error, we can
simulate the process of sampling and analysis, assuming certain
distributions of parasites across the population based on values from
the literature. Analyses of empirical and simulated livestock fecal
egg counts, to assess error associated with composite FEC compared to
individual FEC. The composite method is much more practical to
implement at scale because it greatly reduces the processing time in
the laboratory. However, pooling samples greatly risks underestimating
the burden when parasites are aggregated, especially at low burdens.

The folder contains the following `R` scripts:

`lit_data_analysis.R`: summarizes aggregation patterns based on mean
and variance of FEC reported in the literature.

`sim_main.R`: simulates populations of goats with aggregated FEC values
and individual and composite FEC tests. Across scenarios there are
differences in the mean burden and aggregation patterns, in the
weight and variability of weights in the field and in the laboratory,
as well as in the test sensitivity.

`sim_visualization.R`: creates plots that summarize the simulations

`sim_analysis.R`: summary metrics and statistical analysis of the
influence of different factors on the relative error and the
probability of false negatives.

`empirical_analysis.R`

Fecal egg count
simulation Date: October 2025 Author: Juan S. Vargas Soto Institution:
Queen's University Belfast


