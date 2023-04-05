#' Map predictions for group overall response
#'
#' Given a set of raster layers for fire history (frequency, severity, time
#' since fire), this function identifies the distribution of group overall
#' response values corresponding to the fire regime for each raster cell and
#' returns requested quantiles of that distribution plus the probability of the
#' overall response being zero. Output is a raster with a band for each quantile
#' plus a band for the probability of a zero value (see details). The default
#' quantiles are \code{c(0.05, 0.5, 0.95)}, i.e. the median and central 90\%
#' bounds. Raster processing is performed using \code{terra} package functions
#' and will be run on parallel cores if these are available. The three input
#' raster layers must be spatially aligned, i.e. have identical extent and map
#' projection. This can be checked with the \code{terra} package function
#' \code{\link[terra]{compareGeom}}. The raster cell size is allowed to differ
#' between layers although it's probably best to avoid this whenever possible.
#'
#' The overall response of a group to a particular fire regime, defined in terms
#' of fire frequency (preceding fifty year period), severity of the most recent
#' fire, and time since last fire (years), is treated as a zero-inflated (ZI-)
#' beta distribution described by three parameters: \code{pzero} (probability of
#' zero value); \code{shape1, shape2} (parameters of the beta distribution that
#' represents values greater than zero). The ZI-beta distribution is derived
#' from combining the three triangular distributions that represent the
#' expert-elicited estimates of the group's response to each of the fire regime
#' components. The zero-inflated beta is technically a hurdle model, so the
#' ZI-beta distribution can be interpreted in two parts: the probability of
#' relative abundance or habitat suitability being zero; and the distribution of
#' likely values if habitat suitability is greater than zero. For cases where
#' the probability of zero is substantial, one or more the quantiles summarizing
#' the distribution will be zero, e.g. if \code{pzero > 0.5} then the median
#' value and any lower quantiles will all be zero.
#'
#' @param group Group
#'
#' @param rfrequency Single-band raster layer (\code{terra::rast} object) with
#'   integer cell values for fire frequency in the preceding fifty year period.
#'   Any values greater than the maximum reference value in the expert data
#'   look-up table (\code{\link{GroupOverallResponse}}) will be clamped to the
#'   maximum reference value (currently 10 fires).
#'
#' @param rseverity Single-band raster layer (\code{terra::rast} object) with
#'   integer cell values for the severity of the last fire expressed as FESM
#'   class values (0: unburnt; 2: low; 3: moderate; 4: high; 5: extreme). Note
#'   that FESM class 1 (reserved) is not generally used and will be treated as
#'   unburnt if it appears in the raster. FESM values will be transformed to the
#'   scale used for severity in the expert-elicited data (0-8), such that FESM
#'   values 3, 4 and 5 map to expert scale values 4, 6, and 8. Any raster values
#'   outside the valid range for FESM will result in the function terminating
#'   with an error message warning message being issued and \code{NA} values
#'   being written to the output rasters.
#'
#' @param rtsf Single-band raster layer (\code{terra::rast} object) with
#'   integer cell values for time (years) since the last fire. Any values
#'   greater than the maximum reference value in the expert data look-up table
#'   (\code{\link{GroupOverallResponse}}) will be clamped to the maximum
#'   reference value (currently 40 years).
#'
#' @param probs Numeric vector of one or more probability values for the overall
#'   response quantiles to return. The default is a named vector:
#'   \code{c(lwr90=0.05, median=0.5, upr90=0.95)} for
#'   the median and central 90\% bounds. Element names will be used for the
#'   corresponding output raster layers. If \code{NULL}, \code{NA} or an empty
#'   vector is passed to this argument, no quantiles will be calculated and the
#'   output raster will have a single layer for the probability of zero. If the
#'   elements in the input vector are not named, default names will be used of
#'   the form \code{'q0.05', q0.5, ...}.
#'
#' @param filename Character string. The default is to based the name on the
#'   value of the \code{group} argument, e.g. \code{'group01_overall.tif'}. If
#'   an empty string or \code{NULL} is provided, the output raster will be
#'   returned as a \code{terra::rast} object linked to a temporary file, which
#'   can be then be saved using the \code{terra::writeRaster} function directly.
#'
#' @param cores The number of parallel cores to use for raster processing.
#'   Default is 1 for sequential processing on a single core. If multiple cores
#'   are available this can greatly speed up processing. However, on systems
#'   with many cores available there will usually be a upper limit of cores to
#'   use, beyond which the competition for disk access to write output data
#'   becomes a blockage. From trial and error, this limit seems to be around 50
#'   cores.
#'
#' @param quiet Logical. If \code{TRUE}, no progress messages will be written to
#'   the console. The default is \code{FALSE}.
#'
#' @return A \code{terra::rast} raster with one layer for the probability of
#'   zero value, plus a layer for each requested quantile. The first layer
#'   (probability of zero value) will be named \code{'pzero'}. Further layers
#'   for requested quantiles will use element names for the provided
#'   \code{'probs'} vector, or default names as described for that argument.
#'
#' @seealso \code{\link{qzibeta}} for calculating quantiles of a ZI-beta
#'   distribution, and \code{\link{get_zibeta_ggdata}} to generate a data frame
#'   to use with \code{ggplot} to visualize distributions.
#'
#' @export
#'
predict_raster <- function(group,
                           rfrequency, rseverity, rtsf,
                           probs = c(0.05, 0.5, 0.95),
                           filename = sprintf("group%02d_overall.tif", group),
                           ncores = 1,
                           quiet = FALSE) {

  if (!("terra" %in% installed.packages()[,"Package"])) {
    stop("The terra package needs to be installed to use this function")
  }

  if (length(group) > 1) {
    warning("More than one group specified but predictions will only be returned for the first.",
            immediate. = TRUE)
    group <- group[1]
  }

  if (!.group_is_defined(group)) {
    msg <- glue::glue("Group {group} is not defined in the expert data table GroupOverallResponse")
  }

  # Check rasters are fully aligned
  if (!terra::compareGeom(rfrequency, rseverity, rtsf,
                          crs = TRUE, ext = TRUE, rowcol = FALSE, res = FALSE,
                          stopOnError = FALSE)) {

    stop("The input fire history rasters have differing extents and/or projections")
  }

  # Check probability values for requested quantiles
  if (length(probs) == 0 || is.null(probs) || (length(probs) == 1 && is.na(probs))) {
    q_names <- character(0)
  } else {
    if (anyNA(probs) || (!all(probs > 0 & probs < 1))) {
      stop("All values for argument 'probs' must be non-NA and in the range 0 < p < 1")
    }

    if (length(unique(probs)) < length(probs)) {
      stop("One or more duplicate values supplied to argument 'probs'")
    }

    # Check probs element names
    q_names <- unique(names(probs))
    if (length(q_names) == 0) {
      # No names provided - make some up
      q_names <- sprintf("q%g", probs)
    } else if (length(q_names) < length(probs)) {
      # Some problem with the supplied names
      msg <- glue::glue("Some, but not all, element names for argument 'probs' are either
                       missing or duplicated")
      stop(msg)
    }
  }

  ncores <- ncores[1]
  if (ncores > 1) {
    ncores.available <- parallel::detectCores()
    if (ncores > ncores.available) {
      msg <- glue::glue("More parallel requested than are available ({ncores.available}).
                               You need a bigger computer!")
      stop(msg)
    }
  }


  # Prepare a look-up table of response distribution parameters covering all
  # fire regimes in the range defined by the expert data
  #
  if (!quiet) {
    msg <- glue::glue("  preparing parameters for group {group}")
    message(msg)
  }

  # limits of expert look-up data for each fire regime component
  LIMS <- apply(GroupOverallResponse[, c("frequency", "severity", "tsf")], MARGIN=2, max)

  group_params <- get_zibeta_parameters(grp = group,
                                        frequency = 0:LIMS['frequency'],
                                        severity = c(0,2,4,6,8),  # Constrained to values than map to FESM
                                        tsf = 0:LIMS['tsf']) %>%

    # Recode 0-8 severity values in expert data to standard FESM values
    dplyr::mutate(severity = ifelse(severity <= 2, severity, severity/2 + 1))

  # If quantiles were requested, add them to the look-up table
  if (length(q_names) > 0) {

    # Old school code is simpler for this bit than dplyr `do` or `reframe`
    qdat <- apply(group_params[, c("pzero", "shape1", "shape2")], MARGIN = 1, function(pars) {
      q <- qzibeta(probs, pzero=pars[1], shape1=pars[2], shape2=pars[3])
      names(q) <- q_names
      q
    })

    group_params <- group_params %>% dplyr::bind_cols( as.data.frame(t(qdat)) )
  }


  # Apply prediction function to the input rasters via terra::app
  #
  if (is.null(filename[1]) || filename[1] == "") {
    filename <- ""
    overwrite <- FALSE
  } else {
    overwrite <- TRUE
  }

  wopts <- list(datatype = "FLT4S")
  if (!quiet) {
    wopts$progress <- 1
  }

  if (!quiet) message("  mapping predictions")

  rfire <- c(rfrequency, rseverity, rtsf)
  names(rfire) <- c("frequency", "severity", "tsf")

  terra::app(rfire, fun = .fn_do_predict,
             group_params = group_params,
             max_frequency = LIMS['frequency'],
             max_tsf = LIMS['tsf'],
             q_names = q_names,
             cores = ncores,
             filename = filename,
             overwrite = overwrite,
             wopt = wopts)
}


# Private helper function for predict_raster to be used with `terra::app`.
# First argument x will be a SpatRaster with three layers for frequency, severity and tsf.
#
.fn_do_predict <- function(x, group_params, max_frequency, max_tsf, q_names) {
  res <- rep(NA_real_, 1 + length(q_names))

  if (!anyNA(x)) {
    # clamp fire frequency values to the limit of the reference data
    x[1] <- min(x[1], max_frequency)

    # clamp time since fire values to the limit of the reference data
    x[3] <- min(x[3], max_tsf)

    # Using base R code here rather than dplyr to run on parallel cores
    # without worrying about exports etc.
    i <- which( with(group_params, frequency == x[1] & severity == x[2] & tsf == x[3]) )

    if (length(i) == 1) {
      res <- c(pzero = group_params$pzero[i])

      if (length(q_names) > 0) {
        res <- c(res, unlist(group_params[i, q_names]))
      }
    }
  }

  res
}
