# Functions for calculating composite log-likelihoods under all models
#	of convergent selection: neutral model, independent mutations,
#	standing variant model with source, standing variant model without source,
#	migration model, combinations of convergent selection models (mixed)
#
# Args:
#	numPops: number of populations sampled (both selected and non-selected)
#
#	positions: vector of genomic positions for region
#
#freqs: matrix of allele frequencies for window of interest with
#		dimension numberOfPopulations x numberOfSites
#
# numBins: the number of bins in which to bin alleles a given distance from the proposed
#		selected sites
#   NOTE: must be same as number specified in generating covariance matrices
#
#	*Specify parameter spaces for likelihood calculations*
# *NOTE: must be same as number specified in generating covariance matrices*
#	selSite: vector of positions of proposed selected sites
#	sels: vector of proposed selection coefficients
#	times: vector of proposed time in generations the variant is standing
#		in populations before selection occurs and prior to migration from
#		source population
#	gs: vector of proposed frequencies of the standing variant
#	migs: migration rate (proportion of individuals from source each generation)
#		*Note: cannot be 0
#	sources: vector of proposed source population of the beneficial allele
#		for both migration and standing variant with source models
#		*Note: the source must be a selected population in selPops
#
#
# Inverses and determinants for all relevant models of convergent selection

calcLikelihood_bin <- function(site, selSiteLoc, det_FOmegas, inv_FOmegas, params, neutral = FALSE) {
    # Calculates log-likelihood of data at a given position for models
    #	with a single parameter (independent mutations model)
    #
    # Args:
    #	site: element of vector "positions" of position for log-likelihood to be calculated
    #	selSiteLoc: element of vector "selSite" of proposed location of selected site
    #	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
    #		given model of convergent adaptation
    #	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
    #		given model of convergent adaptation
    #	par1: element of vector of parameter used to specify given model of convergent adaptation
    #		("sels" for independent mutations model)
    #
    # Returns:
    #	log-likelihood of data at a given position

    bin = params$distBins[site, selSiteLoc]
    my.x = as.matrix(params$freqs_MC[ , site])
    my.e = params$epsilons[site] * (1 - params$epsilons[site])

    mu <- params$mu

    if(neutral){

      #likelihood = 1 / (sqrt((2 * pi)^params$k * (det_FOmegas * my.e^ params$rank))) * exp(-1 / 2 * t(my.x - mu) %*% (inv_FOmegas / my.e) %*% (my.x - mu))

      likelihood = -1/2*(params$k*log(2 * pi) + log(det_FOmegas) + params$rank*log(my.e)) + (-1/2 * t(my.x - mu) %*% (inv_FOmegas / my.e) %*% (my.x - mu))

    } else{

      likelihood = -1/2*(params$k*log(2 * pi) + log(det_FOmegas[[bin]]) + params$rank*log(my.e)) + (-1/2 * t(my.x - mu) %*% (inv_FOmegas[[bin]] / my.e) %*% (my.x - mu))

      #likelihood = 1 / (sqrt((2 * pi)^params$k * (det_FOmegas[[bin]] * my.e^ params$rank))) * exp(-1 / 2 * t(my.x - mu) %*% (inv_FOmegas[[bin]] / my.e) %*% (my.x - mu))


    }

    #return(log(likelihood))
    return(likelihood)
}

calcCompLikelihood_par = function(selSiteLoc, det_FOmegas, inv_FOmegas, params, neutral = FALSE) {

  # Calculates composite log-likelihood of all data for models with a
  #	single parameter (independent mutations model)
  #
  # Args:
  #	selSiteLoc: element of vector "selSite" of proposed location of selected site
  #	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
  #		given model of convergent adaptation
  #	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
  #		given model of convergent adaptation
  #	par1: element of vector of parameter used to specify given model of convergent adaptation
  #		("sels" for independent mutations model)
  #
  # Returns:
  #	composite log-likelihood of data under model
  if(neutral){
    all = sapply(1 : length(params$positions), {
      function(i) calcLikelihood_bin(i, selSiteLoc, det_FOmegas, inv_FOmegas, params, neutral = TRUE)
    })
  } else {
  all = sapply(1 : length(params$positions), {
    function(i) calcLikelihood_bin(i, selSiteLoc, det_FOmegas, inv_FOmegas, params)
  })
  }
  return(sum(all))
}
