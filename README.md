# rdmc
## R package for Distinguishing among modes of convergent adaptation using population genomic data.

NOTE: This package is a work in progress, especially the documentation. All effort will be made to keep updates backwards compatible, but there are no guarantees.

This is the source code page for an R package implementing methods presented in Lee and Coop (2017). [See this page for Kristin Lee's original code and exentions](https://github.com/kristinmlee/rdmc/)

While there is a good amount of additions, there are only minimal changes to the original code, and no changes were made to the mathematical underpinnings. Assume mistakes are caused by Silas Tittes and not by Kristin Lee or Graham Coop, and please do report issues to this page.

The code below is a reimplementation of [Kristin Lee's original DMC example.](https://github.com/kristinmlee/rdmc/blob/master/dmc_example.md), using the same example data.

If this package is used, *please* cite Lee and Coop (2017), and share a link this page. A preprint should be available to cite the r package in the near future. 

## Usage

```

devtools::install_github(repo = "RILAB/rdmc")
library(rdmc)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 15))


#load example data
data(neutralAlleleFreqs_example)
data(selectedRegionPositions_example)
data(selectedRegionAlleleFreqs_example)

#specify parameters and input data.
param_list <-
  parameter_barge(
    Ne =  10000,
    rec = 0.005,
    neutral_freqs = allFreqs,
    selected_freqs = freqs_notRand,
    selected_pops = c(1, 3, 5),
    positions = positions,
    n_sites = 10,
    sample_sizes = rep(10, 6),
    num_pops = 6,
    num_bins = 1000,
    sels = c(
      1e-4,
      1e-3,
      0.01,
      seq(0.02, 0.14, by = 0.01),
      seq(0.15, 0.3, by = 0.05),
      seq(0.4, 0.6, by = 0.1)
    ),
    times = c(0, 5, 25, 50, 100, 500, 1000, 1e4, 1e6),
    gs = c(1 / (2 * 10000), 10 ^ -(4:1)),
    migs = c(10 ^ -(seq(5, 1, by = -2)), 0.5, 1),
    sources = selected_pops,
    locus_name = "chow"
  )

#composite likelihood estimates for each model. note: multi core will no work in RStudio
neut_cle <- cle_neutral(param_list)
ind_cle <- cle_ind(param_list)
mig_cle <- cle_mig(param_list)
sv_cle <- cle_svsrc(param_list)

param_list <-
  update_mode(barge = param_list,
              sets = list(c(1, 3), 5),
              modes = c("sv", "ind"))
multi_svind <- cle_multi(param_list)

param_list <- update_mode(barge = param_list, sets = barge$sets, modes =  c("mig", "ind"))
multi_migind <- cle_multi(param_list)

#combine all data sets
mergeby <- names(neut_cle)
all_mods <-
  full_join(ind_cle, mig_cle, by = mergeby) %>%
  full_join(., sv_cle, by = mergeby) %>%
  full_join(., multi_svind, by = mergeby) %>%
  full_join(., multi_migind, by = mergeby)

readr::write_csv(all_mods, "all_mods.csv")

all_mods <- readr::read_csv("all_mods.csv")
#get max comp likelihood parameters estimates for all models
all_mods %>%
  group_by(model) %>%
  filter(cle == max(cle))

neut <- unique(neut_cle$cle)
all_mods %>%
  group_by(selected_sites, model) %>%
  summarise(mcle = max(cle) - neut) %>%
  ggplot(aes(selected_sites, mcle, colour = model)) +
  geom_line() +
  geom_point() +
  xlab("Position") +
  ylab("Composite likelihood") +
  theme(legend.position = "n") +
  scale_color_brewer(palette = "Set1")

#visualize likelihood surface wrt selection coefficients
all_mods %>%
  group_by(sels, model) %>%
  summarise(mcle = max(cle) - neut) %>%
  ggplot(aes(sels, mcle, colour = model)) +
  geom_line() +
  geom_point() +
  ylab("Composite likelihood") +
  xlab("Selection coefficient") +
  scale_color_brewer(palette = "Set1")

```
