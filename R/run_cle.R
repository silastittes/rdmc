#' Generate data frame for neutral model
#'
#'	@param barge List of parameters and data generated using parameter_barge()
#'	@export

cle_neutral <-
  function(barge){

    grid_df <- expand_grid(
      selSite = barge$selSite,
      barge$neut_par)
    grid_df <- mutate(grid_df, site_idx = group_indices(grid_df, selSite))

    neutral_cle <-
      pmap_dfr(grid_df, function(selSite, site_idx, idx){
        cle <- calcCompLikelihood_par(site_idx, barge$det_FOmegas_neutral,  barge$inv_FOmegas_neutral, barge, neutral = TRUE)
        tibble(selSite, cle, locus = barge$locus_name, model = "neutral")
      })

    neutral_cle <- mutate(neutral_cle, sels = as.numeric(NA), gs = as.numeric(NA), times = as.numeric(NA), migs = as.numeric(NA), sources = as.numeric(NA))
    return(neutral_cle)
  }


#' Generate data frame for independent mutations model
#'
#'	@param barge List of parameters and data generated using parameter_barge()
#'	@param cores Number of cores to use. Defaults to 1. if More than one are used, furrr::feautre_pmap is used.
#'	@export

cle_ind <-
  function(barge, cores = 1){
    ind_df <-
      do(group_by(barge$ind_par, idx), {
        FOmegas_ind <- calcFOmegas_indSweeps(.$sels, barge)
        det <- map(FOmegas_ind, ~det(.x))
        inv <- map(FOmegas_ind, ~ginv(.x))
        tribble(
          ~FOmegas_ind, ~det, ~inv,
          FOmegas_ind, det, inv
        )
      })

    grid_df <- expand_grid(
      selSite = barge$selSite,
      barge$ind_par)

    grid_df <- mutate(grid_df, site_idx = group_indices(grid_df, selSite))

    if(cores > 1){
      plan(multicore, workers = cores)
      ind_cle <-
        future_pmap_dfr(grid_df, function(selSite, site_idx, sels, idx){
          cle <- calcCompLikelihood_par(site_idx, ind_df$det[[idx]], ind_df$inv[[idx]], barge)
          tibble(selSite, sels, cle, locus = barge$locus_name, model = "independent")
        })
    } else{
      ind_cle <-
        pmap_dfr(grid_df, function(selSite, site_idx, sels, idx){
          cle <- calcCompLikelihood_par(site_idx, ind_df$det[[idx]], ind_df$inv[[idx]], barge)
          tibble(selSite, sels, cle, locus = barge$locus_name, model = "independent")
        })
    }

    ind_cle <- mutate(ind_cle, gs = as.numeric(NA), times = as.numeric(NA), migs = as.numeric(NA), sources = as.numeric(NA))
    return(ind_cle)
  }

#' Generate data frame for migration from a source population model
#'
#'	@param barge List of parameters and data generated using parameter_barge()
#'	@param cores Number of cores to use. Defaults to 1. if More than one are used, furrr::feautre_pmap is used.
#'	@export

cle_mig <-
  function(barge, cores = 1){

    mig_df <-
      do(group_by(barge$mig_par, idx), {
        FOmegas_mig <- calcFOmegas_mig(.$sels, .$migs, .$sources, barge)
        det <- map(FOmegas_mig, ~det(.x))
        inv <- map(FOmegas_mig, ~ginv(.x))
        tribble(
          ~FOmegas_mig, ~det, ~inv,
          FOmegas_mig, det, inv
        )
      })

    grid_df <-
      expand_grid(
        selSite = barge$selSite,
        barge$mig_par)
    grid_df <- mutate(grid_df , site_idx = group_indices(grid_df, selSite))

    if(cores > 1){
      plan(multicore, workers = cores)
      mig_cle <-
        future_pmap_dfr(grid_df, function(selSite, idx, site_idx, sels, migs, sources){
          cle <- calcCompLikelihood_par(site_idx, mig_df$det[[idx]], mig_df$inv[[idx]], barge)
          tibble(selSite, sels, migs, sources, cle, locus = barge$locus_name, model = "migration")
        })
    } else {
      mig_cle <-
        pmap_dfr(grid_df, function(selSite, idx, site_idx, sels, migs, sources){
          cle <- calcCompLikelihood_par(site_idx, mig_df$det[[idx]], mig_df$inv[[idx]], barge)
          tibble(selSite, sels, migs, sources, cle, locus = barge$locus_name, model = "migration")
        })
    }

    mig_cle <- mutate(mig_cle, gs = as.numeric(NA), times = as.numeric(NA))
    return(mig_cle)
  }


#' Generate data frame for standing variation model
#'
#'	@param barge List of parameters and data generated using parameter_barge()
#'	@param cores Number of cores to use. Defaults to 1. if More than one are used, furrr::feautre_pmap is used.
#'	@export

cle_sv <-
  function(barge, cores = 1){

    sv_df <-
      do(group_by(barge$sv_par, idx), {
        FOmegas_sv <- calcFOmegas_stdVar(.$sels, .$gs, .$times, barge)
        det <- map(FOmegas_sv, ~det(.x))
        inv <- map(FOmegas_sv, ~ginv(.x))
        tribble(
          ~FOmegas_sv, ~det, ~inv,
          FOmegas_sv, det, inv
        )
      })

    grid_df <-
      expand_grid(
        selSite = barge$selSite,
        barge$sv_par)
    grid_df <- mutate(grid_df , site_idx = group_indices(grid_df, selSite))

    if(cores > 1){
      plan(multicore, workers = cores)
      sv_cle <-
        future_pmap_dfr(grid_df, function(selSite, idx, site_idx, sels, gs, times){
          cle <- calcCompLikelihood_par(site_idx, sv_df$det[[idx]], sv_df$inv[[idx]], barge)
          tibble(selSite, sels, gs, times, cle, locus = barge$locus_name, model = "standing_variant")
        })
    } else {
      sv_cle <-
        pmap_dfr(grid_df, function(selSite, idx, site_idx, sels, gs, times){
          cle <- calcCompLikelihood_par(site_idx, sv_df$det[[idx]], sv_df$inv[[idx]], barge)
          tibble(selSite, sels, gs, times, cle, locus = barge$locus_name, model = "standing_variant")
        })
    }

    sv_cle <- mutate(sv_cle, migs = as.numeric(NA), sources = as.numeric(NA))
    return(sv_cle)
  }


#' Generate data frame for standing variation from a source population model.
#'
#'	@param barge List of parameters and data generated using parameter_barge()
#'	@param cores Number of cores to use. Defaults to 1. if More than one are used, furrr::feautre_pmap is used.
#'	@export

cle_svsrc <-
  function(barge, cores = 1){

    sv_df <-
      do(group_by(barge$svsrc_par, idx), {
        FOmegas_sv <- calcFOmegas_stdVar.source(.$sels, .$gs, .$times, .$sources, barge)
        det <- map(FOmegas_sv, ~det(.x))
        inv <- map(FOmegas_sv, ~ginv(.x))
        tribble(
          ~FOmegas_sv, ~det, ~inv,
          FOmegas_sv, det, inv
        )
      })

    grid_df <-
      expand_grid(
        selSite = barge$selSite,
        barge$svsrc_par)
    grid_df <- mutate(grid_df, site_idx = group_indices(grid_df, selSite))

    if(cores>1){
      plan(multicore, workers = cores)
      svsrc_cle <-
        future_pmap_dfr(grid_df , function(idx, site_idx,
                                           selSite, sels, gs, times, sources){
          cle <- calcCompLikelihood_par(site_idx, sv_df$det[[idx]], sv_df$inv[[idx]], barge)
          tibble(selSite, sels, gs, times, sources, cle, locus = barge$locus_name, model = "standing_source")
        })
    } else {
      svsrc_cle <-
        pmap_dfr(grid_df, function(idx, site_idx,
                                   selSite, sels, gs, times, sources){
          cle <- calcCompLikelihood_par(site_idx, sv_df$det[[idx]], sv_df$inv[[idx]], barge)
          tibble(selSite, sels, gs, times, sources, cle, locus = barge$locus_name, model = "standing_source")
        })
    }
    svsrc_cle <- mutate(svsrc_cle, migs = as.numeric(NA))
    return(svsrc_cle)
  }



#' Generate data frame for multiple modes
#'
#'	@param barge List of parameters and data generated using parameter_barge()
#'	@param cores Number of cores to use. Defaults to 1. if More than one are used, furrr::feautre_pmap is used.
#'	@export

cle_multi <-
  function(barge, cores = 1){

    multi_df <-
      do(group_by(barge$multi_par, idx), {
        FOmegas_multi <- calcFOmegas_mixed(.$sels, .$gs, .$times, .$migs, .$sources, barge)
        det <- map(FOmegas_multi, ~det(.x))
        inv <- map(FOmegas_multi, ~ginv(.x))
        tribble(
          ~FOmegas_multi, ~det, ~inv,
          FOmegas_multi, det, inv
        )
      })

    grid_df <-
      expand_grid(
        selSite = barge$selSite,
        barge$multi_par)
    grid_df <- mutate(grid_df, site_idx = group_indices(grid_df, selSite))

    if(cores>1){
      plan(multicore, workers = cores)
      multi_cle <-
        future_pmap_dfr(grid_df , function(idx, site_idx,
                                           selSite, sels, gs, times, migs, sources){
          cle <- calcCompLikelihood_par(site_idx, multi_df$det[[idx]], multi_df$inv[[idx]], barge)
          tibble(selSite, sels, gs, times, migs, sources, cle, locus = barge$locus_name, model = "mixed")
        })
    } else {
      multi_cle <-
        pmap_dfr(grid_df , function(idx, site_idx,
                                    selSite, sels, gs, times, migs, sources){
          cle <- calcCompLikelihood_par(site_idx, multi_df$det[[idx]], multi_df$inv[[idx]], barge)
          tibble(selSite, sels, gs, times, migs, sources, cle, locus = barge$locus_name, model = "mixed")
        })
    }
    return(multi_cle)
  }

