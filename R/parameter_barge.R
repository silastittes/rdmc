

#' Slightly faster Determinant
#' @param Sigma Covariance matrix
#' @noRd
#' @export
chol_det <- function(Sigma){
  exp(2*sum(log(diag(chol(Sigma)))))
  #tryCatch(exp(2*sum(log(diag(chol(Sigma))))),
  #         error = function(c) "Cholesky factorization didn't work. Change 'parameter_barge(..,cholesky = FALSE)'",
  #         warning = function(c) "Cholesky factorization didn't work. Change 'parameter_barge(..,cholesky = FALSE)'"
  #)
}


#' Fast matrix inversion, requires Sigma is singular definite.
#' @param Sigma Covariance matrix
#' @noRd
#' @export
chol_inv <- function(Sigma){
  chol2inv(chol(Sigma))
  #tryCatch(chol2inv(chol(Sigma)),
  #         error = function(c) "Cholesky factorization didn't work. Try 'parameter_barge(..,cholesky = FALSE)'",
  #         warning = function(c) "Cholesky factorization didn't work. Try 'parameter_barge(..,cholesky = FALSE)'"
  #)
}

#' Generate and transfer parameters, quantities, and objects used for fitting downstream convergent adaptation models.
#'
#'
#' @param neutral_freqs Matrix of allele frequencies at putatively neutral sites with dimensions,  number of populations x number of sites.
#' @param selected_freqs Matrix of allele frequencies at putatively selected sites with dimensions, number of populations x number of sites.
#' @param selected_pops Vector of indices for populations that experienced selection.
#' @param positions Vector of genomic positions for the selected region.
#' @param n_sites Integer for the number of sites to propose as the selected site. Sites are uniformly placed along \code{positions} using \code{seq(min(positions), max(positions), length.out = n_sites)}. Must be less than or equal to \code{length(positions)}. Cannot be used with \code{sel_sites}.
#' @param sel_sites Optional vector of sites to propose as selected site. Useful if particular loci are suspected to be under selection. Cannot be used with \code{n_sites}.
#' @param sample_sizes Vector of sample sizes of length number of populations. (i.e. twice the number of diploid individuals sampled in each population).
#' @param num_bins The number of bins in which to bin alleles a given distance from the proposed selected sites.
#' @param sels Vector of proposed selection coefficients.
#' @param times Vector of proposed times in generations the variant is standing in populations before selection occurs and prior to migration from source population.
#' @param gs Vector of proposed frequencies of the standing variant.
#' @param migs Vector of proposed migration rates (proportion of individuals of migrat origin each generation). Cannot be 0.
#' @param sources Vector of population indices to propose as the source population of the beneficial allele. Used for both the migration and standing variant with source models. Note: the source must be one of the populations contained in selPops.
#' @param Ne Effective population size (assumed equal for all populations).
#' @param rec Per base recombination rate for the putatively selected region.
#' @param locus_name String to name the locus. Helpful if multiple loci will be combined in subsequent analyses. Defaults to "locus".
#' @param sets  A list of population indices, where each element in the list contains a vector of populations with a given mode of convergence. For example, if populations 2 and 6 share a mode and population 3 has another, sets = list(c(2,6), 3). Required for modeling multiple modes. Only required for fitting models with mixed modes. Must be used in conjunction with the "modes".
#' @param modes Character vector of length sets defining mode for each set of selected populations ("independent", "standing", and/or "migration"). Only required for fitting models with mixed modes. More details about the modes is available on help page for \code{\link{mode_cle}}
#' @param cholesky Logical to use cholesky factorization of covariance matrix. Used for both inverse and determinant. Faster, but not guaranteed to work for all data sets. TRUE by default. if FALSE, \code{\link[MASS]{ginv}} from MASS is used.
#' @export

parameter_barge <-
  function(neutral_freqs,
           selected_freqs,
           selected_pops,
           positions,
           n_sites = NULL,
           sel_sites = NULL,
           sample_sizes,
           num_bins,
           sets = NULL,
           modes = NULL,
           sels,
           migs,
           times,
           gs,
           sources,
           Ne,
           rec,
           locus_name = "locus",
           cholesky = TRUE) {

    if(dim(neutral_freqs)[1] != dim(selected_freqs)[1]){
      stop("Number of populations for neutral_freqs and selected_freqs do not match.")
    }

    #convert variable names
    allFreqs = neutral_freqs
    freqs_notRand = selected_freqs
    selPops = selected_pops
    sampleSizes = sample_sizes
    numPops = dim(selected_freqs)[1]
    numBins = num_bins


    #generated stuff
    sources = selPops


    if(missing(sel_sites) & missing(n_sites)){
      stop("sel_sites or n_sites argument must be used, but not both! sel_sites requires a vector of length > 0, n_sites requires an integer > 0.")
    } else if(missing(sel_sites) & !missing(n_sites)){
      selSite = seq(min(positions), max(positions), length.out = n_sites)
    } else if(!missing(sel_sites) & missing(n_sites)){
      selSite = sel_sites
    } else{
      stop("sel_sites or n_sites argument must be used, but not both! sel_sites requires a vector of length > 0, n_sites requires an integer > 0.")
    }

    allRunFreq = apply(allFreqs, 2, function(my.freqs) {
      if (runif(1) < 0.5) {
        my.freqs = 1 - my.freqs
      }
      my.freqs
    })


    #Neutral covariance matrix
    numLoci = ncol(allRunFreq)
    my.means.rand = (allRunFreq %*% t(allRunFreq)) / numLoci

    diag(my.means.rand) = diag(my.means.rand) * sampleSizes / (sampleSizes - 1) - rowMeans(allRunFreq) /
      (sampleSizes - 1)

    dist.ij = which(my.means.rand == min(my.means.rand), arr.ind = TRUE)[1, ]

    A.rand = mean(allRunFreq[dist.ij[1], ] * allRunFreq[dist.ij[2], ])
    C.rand = mean(allRunFreq[dist.ij[1], ] * (1 - allRunFreq[dist.ij[2], ]))

    F_estimate = (my.means.rand - A.rand) / C.rand

    M = numPops
    Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M),
                     nrow = M - 1,
                     ncol = M)
    diag(Tmatrix) = (M - 1) / M
    sampleErrorMatrix = diag(1 / sampleSizes, nrow = numPops, ncol = numPops)

    Sigma <- Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix)

    if(cholesky){
      inv_FOmegas_neutral = chol_inv(Sigma)
      det_FOmegas_neutral = chol_det(Sigma)
    } else {
      det_FOmegas_neutral = det(Sigma)
      inv_FOmegas_neutral = ginv(Sigma)
    }



    #grids of parameter combinations to search over for each of the three main models. more to come?
    full_par <- expand_grid(sels, gs, times, migs, sources)

    ind_par <-
      mutate(distinct(dplyr::select(full_par, sels)), idx = 1:n())
    neut_par <-
      mutate(distinct(dplyr::select(ind_par, -sels)), idx = 1:n())
    mig_par <-
      mutate(distinct(dplyr::select(full_par, -c(times, gs))), idx = 1:n())
    sv_par <-
      mutate(distinct(dplyr::select(full_par, -migs, -sources)), idx = 1:n())
    svsrc_par <-
      mutate(distinct(dplyr::select(full_par, -migs)), idx = 1:n())

    if (!missing(modes)) {
      modes_s <- unique(sort(modes))
      multi_par <-
        ifelse(
          identical(modes_s, c("independent", "standing_source")),
          tibble(expand_grid(sels, gs, times, migs = migs[1], sources)),
          ifelse(
            identical(modes_s, c("independent", "migration")),
            tibble(expand_grid(
              sels, gs = gs[1], times = times[1], migs, sources
            )),
            ifelse(
              identical(modes_s, c("migration", "standing_source")),
              tibble(expand_grid(sels, gs, times, migs, sources)),
              ifelse(identical(modes_s, c("independent", "migration", "standing_source")),
              tibble(
                expand_grid(sels, gs, times, migs, sources)
              ),
              NA)
            )
          )
        )[[1]]
      #multi_par <- expand_grid(sels, gs, times, migs, sources)
      multi_par <- mutate(multi_par, idx = 1:n())
    } else {
      multi_par <- NULL
    }


    #matrix goodness1
    nonSelPops = seq(1, numPops)[-selPops]
    distances = sapply(1:length(selSite), function(i)
      abs(positions - selSite[i]))

    ##get distance
    my.seq = seq(min(distances), max(distances), length.out = (numBins + 1))
    midDistances = sapply(1:numBins, function(i)
      mean(c(my.seq[i], my.seq[i + 1])))

    ##MVN parameters
    k = numPops - 1
    mu = as.matrix(rep(0, k))
    rank = numPops - 1

    ##mean centering
    M = numPops
    Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M),
                     nrow = M - 1,
                     ncol = M)
    diag(Tmatrix) = (M - 1) / M

    ##selected matrix goodness##
    randFreqs = apply(freqs_notRand, 2, function(my.freqs) {
      if (runif(1) < 0.5) {
        my.freqs = 1 - my.freqs
      }
      my.freqs
    })

    #get site-specific mean allele frequencies across populations and mean-centered population allele frequencies
    freqs <- randFreqs

    #calculate distances from proposed selected sites and bin
    distances = sapply(1:length(selSite), function(i)
      abs(positions - selSite[i]))
    #numBins = 1000
    my.seq = seq(min(distances) - 0.001,
                 max(distances) + 0.001,
                 length.out = (numBins + 1))
    distBins = apply(distances, 2, function(i)
      as.numeric(cut(i, my.seq)))

    #mean centering
    M = numPops
    Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M),
                     nrow = M - 1,
                     ncol = M)
    diag(Tmatrix) = (M - 1) / M

    #get site-specific mean allele frequencies across populations and mean-centered population allele frequencies
    freqs = t(freqs)
    epsilons = rowMeans(freqs)
    freqs_MC = sapply(1:nrow(freqs), function(i)
      Tmatrix %*% freqs[i,])

    #MVN parameters
    k = numPops - 1
    mu = as.matrix(rep(0, k))
    rank = numPops - 1

    barge_list <-
      list(
        locus_name = locus_name,
        allFreqs = allFreqs,
        freqs_notRand = freqs_notRand,
        positions = positions,
        sampleSizes = sampleSizes,
        selSite = selSite,
        numPops = numPops,
        numBins = numBins,
        n_sites = n_sites,
        selPops = selPops,
        sources = sources,
        sets = sets,
        modes = modes,
        Ne = Ne,
        rec = rec,
        sels = sels,
        times = times,
        gs = gs,
        migs = migs,
        neut_par = neut_par,
        ind_par = ind_par,
        mig_par = mig_par,
        sv_par = sv_par,
        svsrc_par = svsrc_par,
        multi_par = multi_par,
        F_estimate = F_estimate,
        det_FOmegas_neutral = det_FOmegas_neutral,
        inv_FOmegas_neutral = inv_FOmegas_neutral,
        k = k,
        mu = mu,
        rank = rank,
        M = M,
        sampleErrorMatrix = sampleErrorMatrix,
        distances = distances,
        midDistances = midDistances,
        Tmatrix = Tmatrix,
        epsilons = epsilons,
        freqs = freqs,
        freqs_MC = freqs_MC,
        my.seq = my.seq,
        distBins = distBins,
        cholesky = cholesky
      )

    return(barge_list)
  }




#' Update an existing parameter barge with a new mode details.
#'
#'
#' @param barge An existing list object made by calling the parameter_barge() function.
#' @param sets  List of length number of different modes of convergence to be specified vector "modes" where each element in list contains vector of populations with a given single mode of convergence i.e. if populations 2 and 6 share a mode and populations 3 has another, sets = list(c(2,6), 3).
#' @param modes Character vector of length sets defining a new set of mixed modes for each set of selected populations ("independent", "standing_source", and/or "migration"). Other variables will be updated accordingly.
#' @export

update_mode <-
  function(barge, sets, modes) {
    sels <- barge$sels
    gs <- barge$gs
    migs <- barge$migs
    times <- barge$times
    sources <- barge$sources
    #grids of parameter combinations to search over for each of the three main models. more to come?
    full_par <- expand_grid(sels, gs, times, migs, sources)

    ind_par <-
      mutate(distinct(dplyr::select(full_par, sels)), idx = 1:n())
    neut_par <-
      mutate(distinct(dplyr::select(ind_par, -sels)), idx = 1:n())
    mig_par <-
      mutate(distinct(dplyr::select(full_par, -c(times, gs))), idx = 1:n())
    sv_par <-
      mutate(distinct(dplyr::select(full_par, -migs, -sources)), idx = 1:n())
    svsrc_par <-
      mutate(distinct(dplyr::select(full_par, -migs)), idx = 1:n())

    modes_s <- unique(sort(modes))

    if(length(setdiff(modes_s, c("independent", "migration", "standing_source")))){
      stop("Only 'independent', 'migration', or 'standing_source' are allowed as vector elements for the modes argument")
    }


    if(identical(modes_s, c("independent", "standing_source"))){
      multi_par <- tibble(expand_grid(sels, gs, times, migs = migs[1], sources))
    }
    if(identical(modes_s, c("independent", "migration"))){
      multi_par <- tibble(expand_grid(sels, gs = gs[1], times = times[1], migs, sources))
    }
    if(identical(modes_s, c("migration", "standing_source"))){
      multi_par <- tibble(expand_grid(sels, gs, times, migs, sources))
    }
    if(identical(modes_s, c("independent", "migration", "standing_source"))){
      multi_par <- tibble(expand_grid(sels, gs, times, migs, sources))
    }

    multi_par <- mutate(multi_par, idx = 1:n())


    barge$sets <- sets
    barge$modes <- modes
    barge$multi_par <- multi_par

    return(barge)
  }
