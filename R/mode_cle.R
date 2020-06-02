#' Fit composite likelihood model of convergent adaption.
#'
#'	@param barge List of parameters and data generated using \code{\link{parameter_barge}}
#'	@param mode Specify the mode of convergent adaption. Options are:
#'
#'	"neutral" -- Neutral model. Return same composite likelihood across all sites.
#'
#'	"independent" --  Selected mutation occurred independently in all populations specified in \code{selected_pops}
#'
#'	 "migration" -- Selected mutation occurs in a source population (post-divergence) and migrates to other populations specified in \code{selected_pops} during the course of the sweep.
#'
#'	 "standing" -- Selected mutation is present at low frequency in ancestral population and sweeps in \code{selected_pops} populations after divergence.
#'
#'	 "standing_source" -- Selected mutation is present at low frequency, experiences migration from a source population prior to selection, sweeps in \code{selected_pops} after migration.
#'
#'	 "multi" --  Some mix of the modes above, with the exception of "standing", which is not currently supported for mixed modes. See \code{sets} and \code{modes} arguments for \code{\link{parameter_barge}} for proper input format.
#'
#'	@return Regardless of the mode, the function returns a data frame with the same columns. When a column does not pertain to a mode, \code{NA} values are added to ensure the data frames of different modes can be easily joined. The columns are:
#'
#'	selected_sites -- The proposed sites of selection constructed with or provided to \code{\link{parameter_barge}}.
#'
#'	cle -- The composite likelihood for that row's combination of parameters.
#'
#'	locus -- The user specified name of the input locus. Helpful if multiple loci will be combined in later analysis.
#'
#'	sels -- Selection coefficient.
#'
#'	gs -- Allele frequence of the standing variant prior to sweep.
#'
#'	times -- Time in generations the variant is standing in populations before selection occurs and prior to migration from source population.
#'
#'	migs -- Migration rate.
#'
#'	sources -- The population where the new mutation originated.
#'
#'	sel_pops -- Hyphen delimited string of the selected populations.
#'
#'	model -- The mode of convergent adaptation proposed. When multi modes are used, the modes apply to which
#'
#'	@export

mode_cle <-
  function(barge, mode){

    if(mode == "neutral"){

      cmodes <- "neutral"

      grid_df <- expand_grid(
        selSite = barge$selSite,
        barge$neut_par)
      grid_df <- mutate(grid_df, site_idx = group_indices(grid_df, selSite))

      neutral_cle <-
        pmap_dfr(grid_df, function(selSite, site_idx, idx){
          cle <- calcCompLikelihood_par(site_idx, barge$det_FOmegas_neutral,  barge$inv_FOmegas_neutral, barge, neutral = TRUE)
          tibble(selected_sites = selSite, cle, locus = barge$locus_name)
        })

      neutral_cle <- mutate(
        neutral_cle,
        sels = as.numeric(NA),
        gs = as.numeric(NA),
        times = as.numeric(NA),
        migs = as.numeric(NA),
        sources = as.numeric(NA),
        sel_pops = as.character(NA),
        model = as.character(NA),
        model = cmodes
      )
      return(neutral_cle)

    } else if(mode == "independent"){

      cmodes <- "independent"
      cpops <- paste0(barge$selPops, collapse = "-")

      ind_df <-
        do(group_by(barge$ind_par, idx), {
          FOmegas_ind <- calcFOmegas_indSweeps(.$sels, barge)
          if(barge$cholesky){
            det <- map(FOmegas_ind, ~chol_det(.x))
            inv <- map(FOmegas_ind, ~chol_inv(.x))
          } else {
            det <- map(FOmegas_ind, ~det(.x))
            inv <- map(FOmegas_ind, ~ginv(.x))
          }
          tribble(
            ~FOmegas_ind, ~det, ~inv,
            FOmegas_ind, det, inv
          )
        })

      grid_df <- expand_grid(
        selSite = barge$selSite,
        barge$ind_par)

      grid_df <- mutate(grid_df, site_idx = group_indices(grid_df, selSite))

      ind_cle <-
        pmap_dfr(grid_df, function(selSite, site_idx, sels, idx){
          cle <- calcCompLikelihood_par(site_idx, ind_df$det[[idx]], ind_df$inv[[idx]], barge)
          tibble(selected_sites = selSite, sels, cle, locus = barge$locus_name)
        })

      ind_cle <- mutate(ind_cle, gs = as.numeric(NA), times = as.numeric(NA), migs = as.numeric(NA),
                        sources = as.numeric(NA), sel_pops = cpops, model = cmodes)
      return(ind_cle)

    } else if(mode == "migration"){

      cmodes <- "migration"
      cpops <- paste0(barge$selPops, collapse = "-")

      mig_df <-
        do(group_by(barge$mig_par, idx), {
          FOmegas_mig <- calcFOmegas_mig(.$sels, .$migs, .$sources, barge)

          if(barge$cholesky){
            det <- map(FOmegas_mig, ~chol_det(.x))
            inv <- map(FOmegas_mig, ~chol_inv(.x))
          } else {
            det <- map(FOmegas_mig, ~det(.x))
            inv <- map(FOmegas_mig, ~ginv(.x))
          }
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

      mig_cle <-
        pmap_dfr(grid_df, function(selSite, idx, site_idx, sels, migs, sources){
          cle <- calcCompLikelihood_par(site_idx, mig_df$det[[idx]], mig_df$inv[[idx]], barge)
          tibble(selected_sites = selSite, sels, migs, sources, cle, locus = barge$locus_name)
        })

      mig_cle <- mutate(mig_cle, gs = as.numeric(NA), times = as.numeric(NA), sel_pops = cpops, model = cmodes)
      return(mig_cle)

    } else if(mode == "standing"){

      cmodes <- "standing"
      cpops <- paste0(barge$selPops, collapse = "-")

      sv_df <-
        do(group_by(barge$sv_par, idx), {
          FOmegas_sv <- calcFOmegas_stdVar(.$sels, .$gs, .$times, barge)

          if(barge$cholesky){
            det <- map(FOmegas_sv, ~chol_det(.x))
            inv <- map(FOmegas_sv, ~chol_inv(.x))
          } else {
            det <- map(FOmegas_sv, ~det(.x))
            inv <- map(FOmegas_sv, ~ginv(.x))
          }

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

      sv_cle <-
        pmap_dfr(grid_df, function(selSite, idx, site_idx, sels, gs, times){
          cle <- calcCompLikelihood_par(site_idx, sv_df$det[[idx]], sv_df$inv[[idx]], barge)
          tibble(selected_sites = selSite, sels, gs, times, cle, locus = barge$locus_name)
        })

      sv_cle <- mutate(sv_cle, migs = as.numeric(NA), sources = as.numeric(NA), sel_pops = cpops, model = cmodes)
      return(sv_cle)

    } else if(mode == "standing_source"){

      cmodes <- "standing_source"
      cpops <- paste0(barge$selPops, collapse = "-")


      sv_df <-
        do(group_by(barge$svsrc_par, idx), {
          FOmegas_sv <- calcFOmegas_stdVar.source(.$sels, .$gs, .$times, .$sources, barge)

          if(barge$cholesky){
            det <- map(FOmegas_sv, ~chol_det(.x))
            inv <- map(FOmegas_sv, ~chol_inv(.x))
          } else {
            det <- map(FOmegas_sv, ~det(.x))
            inv <- map(FOmegas_sv, ~ginv(.x))
          }

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

      svsrc_cle <-
        pmap_dfr(grid_df, function(idx, site_idx,
                                   selSite, sels, gs, times, sources){
          cle <- calcCompLikelihood_par(site_idx, sv_df$det[[idx]], sv_df$inv[[idx]], barge)
          tibble(selected_sites = selSite, sels, gs, times, sources, cle, locus = barge$locus_name)
        })

      svsrc_cle <- mutate(svsrc_cle, migs = as.numeric(NA), sel_pops = cpops, model = cmodes)
      return(svsrc_cle)

    } else if(mode == "multi"){

      cmodes <- paste0(rep(barge$modes, map_dbl(barge$sets, length)), collapse = "-")
      cpops <- paste0(map_chr(barge$sets, paste, collapse = "_"), collapse = "-")

      multi_df <-
        do(group_by(barge$multi_par, idx), {
          FOmegas_multi <- calcFOmegas_mixed(.$sels, .$gs, .$times, .$migs, .$sources, barge)

          if(barge$cholesky){
            det <- map(FOmegas_multi, ~chol_det(.x))
            inv <- map(FOmegas_multi, ~chol_inv(.x))
          } else {
            det <- map(FOmegas_multi, ~det(.x))
            inv <- map(FOmegas_multi, ~ginv(.x))
          }

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

      multi_cle <-
        pmap_dfr(grid_df , function(idx, site_idx,
                                    selSite, sels, gs, times, migs, sources){
          cle <- calcCompLikelihood_par(site_idx, multi_df$det[[idx]], multi_df$inv[[idx]], barge)
          tibble(selected_sites = selSite, sels, gs, times, migs, sources, cle, locus = barge$locus_name, sel_pops = cpops, model = cmodes)
        })

      multi_cle <- mutate(multi_cle,
                          sels = as.numeric(sels),
                          gs = as.numeric(gs),
                          times = as.numeric(times),
                          migs = as.numeric(migs),
                          sources = as.numeric(sources),
                          sel_pops = as.character(sel_pops),
                          model = as.character(model),
                          sel_pops = as.character(sel_pops),
                          model = cmodes)

      return(multi_cle)
    } else{
      stop('Legal mode arugments must be one of "neutral", "independent", "migration", "standing", "standing_source", or "multi".')
    }

  }

