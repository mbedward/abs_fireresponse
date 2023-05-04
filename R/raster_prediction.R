#' Convert FESM integer fire severity values to scale used for expert-elicited
#' data
#'
#' The FESM scale of fire severity comprises integers in the range 0-5, where 0
#' is unburnt; 2 is low; 3 is moderate; 4 is high; and 5 is extreme. This
#' function transforms these values to the 0-8 scale used in the expert-elicited
#' data for functional group responses, such that FESM values \code{(0,2,3,4,5)}
#' are mapped to expert values \code{(0,2,4,6,8)}.
#'
#' Values of 1 (referred to as 'reserved class' in FESM documentation) do not
#' usually occur in FESM data and will be treated as unburnt (0) by this
#' function.
#'
#' @param sev Vector of integer values for FESM fire severity (0-5).
#'
#' @return A vector of integer severity values on the expert-elicited scale
#'   (0-8).
#'
#' @examples
#' # Convert the commonly used FESM values
#' fesm_to_expert_severity(c(0, 2, 3, 4, 5))
#'
#' # Can also convert FESM value of 1 (not normally used)
#' fesm_to_expert_severity(1)
#'
#' @export
#'
fesm_to_expert_severity <- function(sev) {
  if ( !isTRUE(all.equal(sev, trunc(sev))) ) {
    stop("All input severity values should be integers")
  }

  if (any(sev < 0 | sev > 5)) stop("All input severity values should be in the range 0-5")

  # Recode any input FESM values of 1 as 0 (unburnt)
  sev[sev == 1] <- 0

  ifelse(sev <= 2, sev, 2*(sev-1))
}


#' Convert integer fire severity values on the scale used for expert-elicited
#' data to FESM severity values
#'
#' In the expert-elicited data for functional group responses, fire severity
#' values are coded on a 0-8 scale. This contrasts with the commonly used FESM
#' scale where severity is coded as: 0 (unburnt); 2 (low); 3 (moderate); 4
#' (high); 5 (extreme). Expert scale values \code{(0, 2, 4, 6, 8)} map directly
#' to FESM values \code{(0, 2, 3, 4, 5)}. The remaining expert-scale values
#' represent intermediate severities: 1 is unburnt/low (between FESM 0 and 2); 3
#' is low/moderate (between FESM 2 and 3); 5 is moderate/high (between FESM 3
#' and 4); and 7 is high/extreme (between FESM 4 and 5). The \code{intermediate}
#' argument can be used to control to handle these values, with options to
#' return the higher or lower FESM value, or (by default) to stop with an error
#' if any such values are encountered in the input vector.
#'
#' Values of 1 (referred to as 'reserved class' in FESM documentation) do not
#' usually occur in FESM data and will never be returned by this function.
#'
#' @param sev Vector of integer values for fire severity on the 0-8 scale used
#'   for expert-elicited data.
#'
#' @param intermediate A character value that controls how expert-scale values
#'   \code{(1,3,5,7)} that map to intermediate values on the FESM scale should
#'   be handled. Options are \code{'error'} (default) to stop with an error if
#'   the input vector contains any such values; \code{'up'} to map these values
#'   to the higher corresponding FESM values \code{(2,3,4,5)}; or \code{'down'}
#'   to map these values to the lower corresponding FESM values
#'   \code{(0,2,3,4)}. May be abbreviated.
#'
#' @return A vector of integer severity values on the FESM scale (0-5).
#'
#' @examples
#' # Expert severity codes that map directly to FESM severity codes.
#' expert_to_fesm_severity(c(0, 2, 4, 6, 8))  # returns c(0, 2, 3, 4, 5)
#'
#' # Allow intermediate expert values and map each to the lower
#' # corresponding FESM value
#' expert_to_fesm_severity(c(1, 3, 5, 7), intermediate = "down")
#'
#' # Map intermediate values to the upper corresponding FESM values
#' expert_to_fesm_severity(c(1, 3, 5, 7), intermediate = "up")
#'
#' # If the argument intermediate is not 'up' or 'down', any intermediate
#' # values in the input vector will cause the function to stop with
#' # an error message
#' \dontrun{
#' expert_to_fesm_severity(3)
#' }
#'
#' @export
#'
expert_to_fesm_severity <- function(sev, intermediate = c("error", "up", "down")) {
  if ( !isTRUE(all.equal(sev, trunc(sev))) ) {
    stop("All input severity values should be integers")
  }

  if (any(sev < 0 | sev > 8)) stop("All input severity values should be in the expert-elicited range: 0-8")

  intermediate <- match.arg(intermediate)
  INTERMEDIATE <- c(1,3,5,7)

  if (intermediate == "error") {
    if (any(INTERMEDIATE %in% sev)) {
      stop("Input data contains one or more intermediate values (1,3,5,7).\n",
           "Argument intermediate must be set to 'up' or 'down' to allow such values.")
    }
  } else {
    ii <- sev %in% INTERMEDIATE
    if (intermediate == "up") {
      sev[ii] <- sev[ii] + 1
    } else { # 'down'
      sev[ii] <- sev[ii] - 1
    }
  }

  # Convert to FESM scale values (0, 2:5)
  ifelse(sev == 0, sev, sev/2 + 1)
}


#' Map predictions for group overall response
#'
#' Given a set of raster layers for fire history (frequency, severity, time
#' since fire), this function maps predicted values of overall response for a
#' given functional group. These values can be thought of as habitat quality due
#' to fire regime (0-1 scale).  For each raster cell, the prediction is based on
#' fire frequency over the past 50 years, severity of the most recent fire, and
#' time since last fire. Predictions are mapped as a set of output rasters (e.g.
#' median value plus 5\% and 95\% quantiles for lower and upper 90% bounds on
#' the prediction) as specified using the \code{probs} argument. An additional
#' raster for the probability of zero habitat quality is always produced.
#' Optionally, predictions can be mapped for a future fire scenario: either (a)
#' no fire for a given number of years from the present; (b) a future fire of
#' specified uniform severity over the entire area; or (c) a future fire in
#' which severity varies over the study area. See Details (below) for how such
#' future fire predictions are made. Raster processing is performed using
#' \code{terra} package functions and will be run on parallel cores if these are
#' available. The three input raster layers (plus the optional future fire
#' severity raster if applicable) must be spatially aligned, i.e. have identical
#' extent and map projection. This can be checked with the \code{terra} package
#' function \code{\link[terra]{compareGeom}}. The raster cell size is allowed to
#' differ between layers although it's probably best to avoid this whenever
#' possible.
#'
#' Predictions can be made for future fire scenarios in one of two ways.
#' Firstly, a prediction for group habitat status after a period of N years from
#' the present with no additional fire can be made by setting argument
#' \code{next_time} to an integer value greater than zero, e.g.
#' \code{next_time=5} for a fire-free period of 5 years. For this case, the
#' function simply adds the value \code{N} to the value of each cell in the
#' \code{rtsf} input raster, other than cells with missing (\code{NA}) values.
#'
#' Secondly, a prediction can be made for habitat quality given a future fire in
#' N years time, specified via argument \code{next_time}, and with uniform
#' severity, specified via argument \code{next_severity}. As with the first
#' option, only those cells that have non-missing values for their actual time
#' since fire in the \code{rtsf} layer will be considered. The value of
#' \code{next_time} can be 0, to simulate the additional fire this year, or an
#' integer N greater than 0 for a fire in N years time. With this option,
#' habitat quality is first calculated for the actual fire history of each
#' raster cell (i.e. based on the values in the three input rasters), and then
#' the modifying effect of the additional fire is calculated. This two-step
#' process is intended to avoid unrealistic predictions that could otherwise
#' result. For example, if a cell currently has low or zero habitat quality
#' because it was recently burnt at high severity, a further low severity burn
#' could appear to greatly increase the habitat quality if the calculation was
#' only based on the severity of the most recent fire.
#'
#' @param group Functional group integer ID.
#'
#' @param rfrequency Single-band raster layer (\code{terra::rast} object) with
#'   integer cell values for fire frequency in the preceding fifty year period.
#'   Any values greater than the maximum reference value in the expert data
#'   look-up table (\code{\link{GroupOverallResponse}}) will be clamped to the
#'   maximum reference value (currently 10 fires).
#'
#' @param rseverity Single-band raster layer (\code{terra::rast} object) with
#'   integer cell values for the severity of the last fire. By default, cell
#'   values are assumed to be FESM severity values (0: unburnt; 2: low; 3:
#'   moderate; 4: high; 5: extreme). Note that FESM class 1 (reserved) is not
#'   generally used and will be treated as unburnt (0) if it occurs in the
#'   raster. FESM values will be transformed to the scale used for severity in
#'   the expert-elicited data, such that FESM values \code{(0,2,3,4,5)} map to
#'   expert scale values \code{(0,2,4,6,8)}. If the raster cell values are already on
#'   the scale used for expert-elicited data, this can be indicated via the
#'   \code{severity_scale} argument (see below).
#'
#' @param rtsf Single-band raster layer (\code{terra::rast} object) with
#'   integer cell values for time (years) since the last fire. Any values
#'   greater than the maximum reference value in the expert data look-up table
#'   (\code{\link{GroupOverallResponse}}) will be clamped to the maximum
#'   reference value (currently 40 years).
#'
#' @param next_time Time (integer years) for a future fire scenario. If provided
#'   and \code{next_severity} is \code{NULL}, the value will be added to the
#'   time since fire value for all non-missing cells in raster \code{rtsf} to
#'   simulate an expected fire-free period. If \code{next_severity} is also
#'   provided, then \code{next_time} is taken as the time from present at which
#'   the future fire burning all raster cells with that severity will occur. In
#'   the latter case, group habitat status will be calculated using the two-step
#'   process described in the Details section.
#'
#' @param next_severity Severity of a hypothetical future fire. If this argument
#'   is provided together with a value for \code{next_time}, group habitat
#'   status will be calculated using the two-step process described in the
#'   Details section. Note: this value must be an integer greater than zero
#'   referring to a severity class on the scale specified via the
#'   \code{severity_scale} argument (see below).
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
#' @param severity_scale (character) Either \code{'FESM'} (default) to indicate
#'   that severity raster values are on the FESM (0-5) scale, or \code{'Expert'}
#'   to indicate that values are on the 0-8 scale used for expert-elicited data.
#'   Case-insensitive and may be abbreviated. Note: this argument applies to
#'   both the severity raster \code{rseverity} as well as the value of
#'   \code{next_severity}, if provided.
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
                           next_time = NULL,
                           next_severity = NULL,
                           probs = c(0.05, 0.5, 0.95),
                           filename = sprintf("group%02d_overall.tif", group),
                           severity_scale = c("FESM", "Expert"),
                           ncores = 1,
                           quiet = FALSE) {

  if (!("terra" %in% installed.packages()[,"Package"])) {
    stop("The terra package needs to be installed to use this function")
  }

  # limits of expert look-up data for each fire regime component
  LIMS <- apply(GroupOverallResponse[, c("frequency", "severity", "tsf")], MARGIN=2, max)

  checkmate::assert_int(group, lower=1, upper=max(get_group_ids()))

  .assert_rast_or_int(rfrequency)
  .assert_rast_or_int(rseverity)
  .assert_rast_or_int(rtsf)

  # Check rasters are fully aligned
  if (!terra::compareGeom(rfrequency, rseverity, rtsf,
                          crs = TRUE, ext = TRUE, rowcol = FALSE, res = FALSE,
                          stopOnError = FALSE)) {

    stop("The input fire history rasters have differing extents and/or projections")
  }

  if (!is.null(next_time)) {
    checkmate::assert_int(next_time, lower=0)
  }

  if (!is.null(next_severity)) {
    # In case next_time was NULL for the check above
    if (is.null(next_time)) stop("Using next_severity but argument next_time has not been set")

    checkmate::assert_int(next_severity, lower=1)
  }

  if (is.null(probs)) {
    qnames = character(0)

  } else {
    checkmate::assert_numeric(probs, lower=0, upper = 1, any.missing = FALSE)

    if (length(unique(probs)) < length(probs)) {
      stop("One or more duplicate values supplied to argument 'probs'")
    }

    # Check and/or set probs element names to use for output layers
    if (length(probs) > 0) {
      q_names <- unique(names(probs))        # will be NULL if there are no element names
      num_names <- length(na.omit(q_names))  # will be zero ...

      if (num_names > 0) {
        if (num_names < length(probs)) {
          # Some problem with the supplied names
          msg <- glue::glue("Some, but not all, element names for argument 'probs' are either
                         missing or duplicated. Setting default names for output layers.")
          warning(msg, immediate. = TRUE)

          q_names <- NULL
        }
      }

      if (is.null(q_names)) q_names <- sprintf("q%g", probs)
    }
  }

  # Selected severity scale - allowing for variable case and partial strings
  ptn <- paste0("^", tolower(severity_scale[1]))
  i <- grep(ptn, c("fesm", "expert"))
  if (length(i) != 1) stop("Unknown option for argument severity_scale: ", severity_scale[1])
  is_FESM_scale <- i == 1

  # Valid severity values
  if (is_FESM_scale) {
    # Expert-scale values that map directly to FESM values
    SEVERITIES <- c(0, 2, 4, 6, 8)
  } else {
    SEVERITIES <- 0:LIMS['severity']
  }

  # Now that we know the scale, scale and validate the 'next_severity' value if one
  # was provided
  if (!is.null(next_severity)) {
    if (is_FESM_scale) {
      ok <- next_severity %in% 1:5
    } else {
      ok <- next_severity %in% 1:8
    }

    if ( !ok ) {
      msg <- paste("Invalid next_severity value for")

      if (is_FESM_scale) msg <- paste(msg, "FESM scale:")
      else msg <- paste(msg, "Expert data scale:")

      msg <- paste(msg, next_severity)

      stop(msg)
    }

    # put next_severity on the expert scale for the look-up
    # table steps (below)
    if (is_FESM_scale) next_severity <- fesm_to_expert_severity(next_severity)
  }


  checkmate::assert_int(ncores)
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


  # Look-up table of distributions for functional group responses for
  # possible values of the current/realized fire regime
  group_params <- get_zibeta_parameters(grp = group,
                                        frequency = 0:LIMS['frequency'],
                                        severity = SEVERITIES,
                                        tsf = 0:LIMS['tsf'])

  if (is.null(next_severity)) {
    # Calculate requested quantiles for group status based on the beta
    # distribution parameters in the group_params table.
    #
    if (length(q_names) > 0) {

      qdat <- apply(group_params[, c("pzero", "shape1", "shape2")], MARGIN = 1, function(pars) {
        q <- qzibeta(probs, pzero=pars[1], shape1=pars[2], shape2=pars[3])
        names(q) <- q_names
        q
      })

      group_params <- group_params %>% dplyr::bind_cols( as.data.frame(t(qdat)) )
    }

  } else {
    # A value was provided for next_severity, so calculate requested quantiles
    # and the probability of zero status based on the product of the beta
    # distribution for the actual fire history and the distribution that
    # would result from the additional fire occurring `next_time` years from now.

    # First get ZI-beta parameters for all possible future regime states
    next_params <- group_params %>%
      dplyr::select(group, frequency) %>%

      dplyr::distinct() %>%

      dplyr::mutate(frequency = frequency + 1,
                    severity = next_severity,
                    tsf = next_time) %>%

      # clamp fire frequency to the limit of the expert-data
      dplyr::filter(frequency <= LIMS['frequency'])

    next_params <- get_zibeta_parameters(grp = group,
                                         frequency = next_params$frequency,
                                         severity = next_params$severity,
                                         tsf = next_params$tsf,
                                         cross = FALSE) %>%
      dplyr::mutate(index = dplyr::row_number())

    # Vector of random samples from each 'next' distribution
    NSAMP <- 1e5
    set.seed(123)

    next_rand <- lapply(seq_len(nrow(next_params)), function(i) {
      with(next_params, rzibeta(NSAMP, pzero[i], shape1[i], shape2[i]))
    })

    # Calculate product of prior and next distributions and record
    # pzero and requested quantiles
    res <- lapply(seq_len(nrow(group_params)), function(i) {
      prior_pzero <- group_params$pzero[i]
      prior_shape1 <- group_params$shape1[i]
      prior_shape2 <- group_params$shape2[i]

      r <- rzibeta(NSAMP, prior_pzero, prior_shape1, prior_shape2)

      # Note - clamping frequency to the limit of the expert data here
      i_next <- which(next_params$frequency == min(LIMS['frequency'], group_params$frequency[i] + 1))

      rprod <- r * next_rand[[i_next]]

      res <- data.frame(next_pzero = mean(rprod == 0))

      if (length(q_names) > 0) {
        q <- quantile(rprod, probs = probs, names = FALSE)
        q <- data.frame( t(q) )
        colnames(q) <- paste("next", q_names, sep = "_")
      }

      cbind(res, q)
    })

    res <- do.call(rbind, res)
    group_params <- cbind(group_params, res)
  }

  if (is_FESM_scale) {
    # Recode severity values in the look-up table to standard FESM values so
    # that they will correspond to cell values in the severity raster
    group_params <- group_params %>%
      dplyr::mutate(severity = expert_to_fesm_severity(severity))

    if (!is.null(next_severity)) next_severity <- expert_to_fesm_severity(next_severity)
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
             next_time = next_time,
             next_severity = next_severity,

             cores = ncores,
             filename = filename,
             overwrite = overwrite,
             wopt = wopts)
}


# Private helper function for predict_raster to be used with `terra::app`.
# First argument x will be a vector of integer values for frequency, severity and tsf.
#
.fn_do_predict <- function(x, group_params, max_frequency, max_tsf, q_names, next_time, next_severity) {
  res <- rep(NA_real_, 1 + length(q_names))

  if (is.null(next_time)) next_time <- 0

  if (is.null(next_severity)) {
    pzero_name <- "pzero"
  } else {
    pzero_name <- "next_pzero"
    q_names <- paste("next", q_names, sep = "_")
  }

  if (!anyNA(x)) {
    # clamp fire frequency values to the limit of the reference data
    x[1] <- min(x[1], max_frequency)

    # clamp time since fire values to the limit of the reference data
    # after adding the value of next_time
    x[3] <- min(x[3] + next_time, max_tsf)

    # Using base R code here rather than dplyr to run on parallel cores
    # without worrying about exports etc.
    i <- which( with(group_params, frequency == x[1] & severity == x[2] & tsf == x[3]) )

    if (length(i) == 1) {
      res <- c(pzero = group_params[i, pzero_name])

      if (length(q_names) > 0) {
        res <- c(res, unlist(group_params[i, q_names]))
      }
    }
  }

  res
}


.assert_rast_or_int <- function(x) {
  xname <- checkmate::vname(x)

  if (!(checkmate::test_class(x, "SpatRaster") || checkmate::test_int(x, lower = 0))) {
    stop(xname, "must be either a SpatRaster object or an integer value >= 0")
  }
}
