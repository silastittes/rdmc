# dmc
## Distinguishing among modes of convergent adaptation using population genomic data -- R package

NOTE: This package is very much a work in progress, especially the documentation. All effort will be made to keep updates backwards compatible, but there are no guarantees.

This is the source code page for an R package implementing methods presented in Lee and Coop (2017). [See this page for Kristin Lee's original code and exentions](https://github.com/kristinmlee/dmc/)

There are only minimal changes to the original code, and no changes were made to the mathematical underpinnings. Assume mistakes are caused by me and not Lee and Coop, and please do report issues to this page.  

The code below is a reimplementation of [Kristin Lee's original DMC example.](https://github.com/kristinmlee/dmc/blob/master/dmc_example.md).

If this package is used, please cite Lee and Coop (2017), and share the link this page (maybe?). 

```
devtools::install_github(repo = "RILAB/rDMC")
library(MASS)
library(tidyverse)
library(furrr)
library(dmc)

data(neutralAlleleFreqs_example)
data(selectedRegionPositions_example)
data(selectedRegionAlleleFreqs_example)

#specify parameters and input data.
barge <- 
  parameter_barge(
    Ne =  10000,
    rec = 0.005,
    allFreqs = allFreqs,  
    freq_notRand = freq_notRand, 
    selPops = c(1, 3, 5),
    positions = positions,
    sets = list(c(1, 3), 5),
    modes = c("sv", "ind"),
    n_sites = 10,
    sampleSizes = rep(10, 6),
    numPops = 6,
    numBins = 1000,
    sels = c(1e-4, 1e-3, 0.01, seq(0.02, 0.14, by = 0.01), seq(0.15, 0.3, by = 0.05), 
             seq(0.4, 0.6, by = 0.1)),
    times = c(0, 5, 25, 50, 100, 500, 1000, 1e4, 1e6),
    gs = c(1/(2*10000), 10^-(4:1)),
    migs = c(10^-(seq(5, 1, by = -2)), 0.5, 1),
    sources = selPops, 
    locus_name = "chow"
  )

#composite likelihood estimates for each model. Multiple cores will not be used if running code in Rstudio.
neut_cle <- cle_neutral(barge)
ind_cle <- cle_ind(barge, cores = 1)
mig_cle <- cle_mig(barge, cores = 1)
sv_cle <- cle_svsrc(barge, cores = 1)
multi_cle <- cle_multi(barge, cores = 1)

#combine all data sets
mergeby <- names(neut_cle)
all_mods <- 
  full_join(mig_cle, ind_cle, by = mergeby) %>%
  full_join(., multi_cle, by = mergeby) %>% 
  full_join(., sv_cle, by = mergeby)


#get max comp likelihood parameters estimates for all models
all_mods %>% 
  group_by(model) %>% 
  filter(cle == max(cle))

#visualize comp likelihoods by position relative to neutral sites
neut <- neut_cle$cle[1]

all_mods %>% 
  group_by(selSite, model) %>% 
  summarise(mcle = max(cle) - neut) %>% 
  ggplot(aes(selSite, mcle, colour = model)) +
  geom_line() +
  geom_point() +
  theme_bw()


#visualize likelihood surface wrt selection coefficients
all_mods %>% 
  group_by(sels, model) %>% 
  summarise(mcle = max(cle) - neut) %>% 
  ggplot(aes(sels, mcle, colour = model)) +
  geom_line() +
  geom_point() +
  theme_bw()

```
