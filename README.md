# compFEC
Simulation of livestock fecal egg count to assess error associated with composite method with respect to individual method

Fecal egg count simulation
Date: October 2025
Author: Juan S. Vargas Soto
Institution: Queen's University Belfast


This code simulates sampling feces in livestock farms, specifically from goats, and processing them in laboratory to quantify the number of helminth eggs. Fecal egg counts are a commonly used technique to estimate the parasite burden in farms. 
Obtaining and analyzing individual samples often is coslty and time-consuming, both for farmers who must collect the samples and for laboratory analysts who must perform many egg counts. Alternatively, collecting many samples together, mixing them in the laboratory greatly reduces the number of samples analyzed in the lab, while still providing an estimate of the mean burden across all individuals in the farm.
Parasites are commonly aggregated, with few individuals harboring most of the adult parasites, and this could create large errors in the estimation of the mean burden. To assess the magnitude and relevance of the potential error, we can simulate the process of sampling and analysis, assuming certain distributions of parasites across the population based on values from the literature.