# rdmc
## R package for Distinguishing among modes of convergent adaptation using population genomic data.


This package is a work in progress, especially the documentation. All effort will be made to keep updates backwards compatible, but there are no guarantees.

This is the source code page for an R package implementing methods presented in Lee and Coop (2017). [See this page for Kristin Lee's original code and exentions](https://github.com/kristinmlee/rdmc/)

While there is a good amount of additions, there are only minimal changes to the original code, and no changes were made to the mathematical underpinnings. Assume mistakes are caused by Silas Tittes and not by Kristin Lee or Graham Coop, and please do report issues to this page.

The code below is a reimplementation of [Kristin Lee's original DMC example.](https://github.com/kristinmlee/rdmc/blob/master/dmc_example.md), using the same example data.


A manuscript describing this package is in revision. An earlier version is available [here](https://www.biorxiv.org/content/10.1101/2020.04.22.056150v1). If the package is used, *please* cite Lee and Coop (2017). Full citation information is available from the R package (after installation) by running the command `citation("rdmc")`.


## Installation

```
devtools::install_github('silastittes/rdmc')
```

*NOTE*: I have made several updates to the package compared to the version currently available on BioXiv. If you need to use the older version for some reason, please modify the above installation to:


```
devtools::install_github('silastittes/rdmc', ref = 'ginv')
```

The output from the two versions should be the same, but the master version is faster and has much better function documentation.

## Usage

```
library(rdmc)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 15))


#load example data
data(neutral_freqs)
data(selected_freqs)
data(positions)


param_list <-
  parameter_barge(
    Ne =  10000,
    rec = 0.005,
    neutral_freqs = neutral_freqs,
    selected_freqs = selected_freqs,
    selected_pops = c(1, 3, 5),
    positions = positions,
    n_sites = 10,
    sample_sizes = rep(10, 6),
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
    locus_name = "test_locus",
    cholesky = TRUE
  )



#fit composite likelihood models
neut_cle <- mode_cle(param_list, mode = "neutral")
ind_cle <- mode_cle(param_list, mode = "independent")
mig_cle <- mode_cle(param_list, mode = "migration")
sv_cle <- mode_cle(param_list, mode = "standing_source")



#update barge to fit a  mixed-mode model
param_list <-
  update_mode(
    barge = param_list,
    sets = list(c(1, 3), 5),
    modes = c("standing_source", "independent"))

#fit mixed-mode model
multi_svind <- mode_cle(param_list, "multi")

#update to another mixed-mode
param_list <-
  update_mode(
    barge = param_list,
    sets = list(c(1, 3), 5),
    modes =  c("migration", "independent"))

#fit mixed-mode model
multi_migind <- mode_cle(param_list, "multi")


#merge data frame of all fit models
all_mods <-
  bind_rows(
    ind_cle,
    mig_cle,
    sv_cle,
    multi_svind,
    multi_migind
  )


```


## Summarize


```
#max composite likelihood estimate 
#of all params over all models  
all_mods %>%
  group_by(model) %>%
  filter(cle == max(cle))

```


```
# A tibble: 5 x 10
# Groups:   model [5]
  selected_sites  sels   cle locus            gs times     migs sources sel_pops model                                      
           <dbl> <dbl> <dbl> <chr>         <dbl> <dbl>    <dbl>   <dbl> <chr>    <chr>                                      
1         0.0017  0.03 3708. test_locus NA          NA NA            NA 1-3-5    independent                                
2         0.0017  0.03 3404. test_locus NA          NA  0.00001       3 1-3-5    migration                                  
3         0.0017  0.05 3733. test_locus  0.01     1000 NA             3 1-3-5    standing_source                            
4         0.0017  0.03 3745. test_locus  0.01    10000  0.00001       1 1_3-5    standing_source-standing_source-independent
5         0.0017  0.03 3605. test_locus  0.00005     0  0.00001       1 1_3-5    migration-migration-independent
```

## Visualize

```
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
  geom_vline(xintercept = 0.05, lty = 2) +
  ylab("Composite likelihood") +
  xlab("Selection coefficient") +
  scale_color_brewer(palette = "Set1")

```

![](https://github.com/silastittes/rdmc/blob/choleskyd/man/figures/cle.png?raw=true)


![](https://github.com/silastittes/rdmc/blob/choleskyd/man/figures/sel.png?raw=true)


The figures show results for the example data. The first shows the composite likelihood score at each of the 10 proposed sites of selection for each model. The true selected site was modeled at position 0.  These figures match the original example generated by Kristin Lee [here](https://github.com/kristinmlee/dmc/blob/master/dmc_example.md).
