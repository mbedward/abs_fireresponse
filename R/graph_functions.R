#' Draw response curves and bounding regions for a given group
#'
#' This is a convenience function to graph the response of a specified group to
#' each of the three component variables of fire regime: TSF (time since last
#' fire), severity of the most recent fire and frequency of fires within a
#' preceding fifty year period. For each component, the response curve is drawn
#' with a bounding region, representing the point-wise degree of expert
#' certainty at each discrete value.
#'
#' For each record in the data frame, the lower, most likely and upper estimates
#' for relative abundance are presently interpreted as the minimum, modal and
#' maximum parameters of a triangular distribution.
#'
#' @param the_group Integer group number: a single value between 1 and the
#'   number of groups defined in \code{\link{GroupExpertData}} (currently 18).
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes facet_wrap geom_line geom_point geom_ribbon labs scale_x_continuous
#'
#' @examples
#' # Graph the response curves for group 3
#' draw_response_curves(3)
#'
#' # Change the default labels and add a title
#' library(ggplot2)
#' draw_response_curves(3) +
#'   labs(x = "Value", title = "Functional group 3")
#'
#' @export
#'
draw_response_curves <- function(the_group) {
  dat_gg <- dplyr::filter(GroupExpertData, group == the_group)

  ggplot(dat_gg, aes(x = value)) +
    geom_ribbon(aes(ymin = ra_lwr, ymax = ra_upr), fill = "darkred", alpha = 0.2) +
    geom_line(aes(y = ra_mode), colour = "darkred") +
    geom_point(aes(y = ra_mode), colour = "darkred") +

    scale_x_continuous(breaks = ~round(unique(pretty(.)))) +

    labs(x = "Fire component value", y = "Relative abundance") +

    facet_wrap(~type, scales = "free_x")
}



#' Graph group overall response to fire regime(s)
#'
#' Presently, this just graphs the beta distribution that approximates the
#' overall group response for a given fire regime. If multiple groups and/or
#' fire regimes are specified, the graph for each group is displayed as a
#' separate facet with the curve for each fire regime identified by colour.
#'
#' @param the_group One or more integer identifiers for group. Each value must
#'   be between 1 and the number of groups defined in
#'   \code{\link{GroupExpertData}} (currently 18).
#'
#' @param frequency Integer vector with one or more values of fire frequency in
#'   the preceding fifty year period.
#'
#' @param severity Integer vector with one or more values of severity for the
#'   most recent fire.
#'
#' @param tsf Integer vector with one or more values of time since last fire
#'   (years).
#'
#' @param cross Logical value. If \code{TRUE} (default), a set of fire regimes
#'   will be created for all combinations of the fire component vectors:
#'   frequency, severity and tsf values. In this case, vectors can have
#'   different lengths; e.g. 2 frequency values x 3 severity values x 2 tsf
#'   values will result in 12 fire regimes. If \code{FALSE}, a fire regime
#'   \emph{i} is defined by the \emph{i}th element of each of the three
#'   component vectors. In this case, the vectors must all have the same length.
#'
#' @param linewidth Value to pass to ggplot for line width (default = 1).
#'
#' @param draw_pzero Logical. If \code{TRUE} (default), a separate graph showing
#'   the probability of zero relative abundance under each fire regime is drawn
#'   and placed above the graph of overall response distribution using the
#'   \code{patchwork} package.
#'
#' @param colour_map Character string indicating the colour map option to use
#'   with the viridis palette. Default is \code{"cividis"}. See
#'   \code{\link[ggplot2]{scale_colour_viridis_d}} for details.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap label_both scale_x_continuous labs
#'
#' @export
#'
draw_overall_response <- function(the_group,
                                  frequency, severity, tsf,
                                  cross = TRUE,
                                  linewidth = 1,
                                  draw_pzero = TRUE,
                                  colour_map = "viridis") {

  dat <- get_overall_response_data(the_group, frequency, severity, tsf, cross)

  gg_dens <- .do_draw_overall_response(dat,
                                       linewidth = linewidth,
                                       colour_map = colour_map)
  if (draw_pzero) {
    gg_zero <- .do_draw_prob_zero(dat, colour_map = colour_map)

    patchwork::wrap_plots(gg_zero, gg_dens, ncol = 1, heights = c(2,8))

  } else {
    gg_dens
  }
}


# Private helper function to draw overall response densities
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous scale_colour_viridis_d
#' @importFrom ggplot2 labs facet_wrap label_both
#
.do_draw_overall_response <- function(dat,
                                      linewidth = 1,
                                      colour_map = "cividis") {

  ngroups <- dplyr::n_distinct(dat$group)
  nregimes <- dplyr::n_distinct(dat$regime)

  gg <- ggplot(data = dat, aes(x = relabund, y = density)) +
    geom_line(aes(colour = regime), linewidth=linewidth) +
    scale_colour_viridis_d(option = colour_map, direction = -1, end = 0.8) +

  if (ngroups > 1) {
    gg <- gg + facet_wrap(~group, labeller = label_both, scales = "free_y")
  }

  gg <- gg +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    labs(x = "Relative abundance", y = "Density")

  gg
}


# Private helper function to draw prob zero graph
#' @importFrom ggplot2 ggplot aes geom_point geom_segment scale_x_continuous scale_colour_viridis_d
#' @importFrom ggplot2 labs facet_wrap label_both theme
#
.do_draw_prob_zero <- function(dat, colour_map = "cividis") {
  dat_zero <- dat %>%
    dplyr::group_by(group, regime) %>%
    dplyr::summarize(pzero = dplyr::first(pzero))

  ngroups <- dplyr::n_distinct(dat$group)

  gg <- ggplot(data = dat_zero, aes(y = regime)) +
    geom_segment(x = 0, aes(xend = pzero, yend = regime, colour = regime), linewidth = 2) +
    geom_point(aes(x = pzero, colour = regime), size = 4) +
    scale_colour_viridis_d(option = colour_map, direction = -1, end = 0.8) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    labs(x = "Probability of zero relative abundance", y = "") +
    theme(legend.position = "none")

  if (ngroups > 1) {
    gg <- gg + facet_wrap(~group, labeller = label_both, scales = "free_y")
  }

  gg
}




#' Prepare data to graph group overall response to fire regime(s)
#'
#' This is a convenience function that creates a data frame suitable for
#' graphing the distribution of expected relative abundance values for one or
#' more fire regimes. It is mainly intended as a helper function to be called by
#' the \code{\link{draw_overall_response}} function, but can also be called
#' explicitly.
#'
#' Where a requested combination of frequency, severity and time since fire does
#' not exist in the \code{GroupOverallResponse} look-up table, the overall
#' response will be derived as the weighted mixture of beta distributions for
#' the nearest defined combinations. At the moment, this interpolation between
#' values is relevant to time since fire, since all integer values for frequency
#' and severity between zero and the maximum values considered during
#' expert-elicitation are represented in the look-up table.
#'
#' @param grp One or more integer identifiers for group. Each value must
#'   be between 1 and the number of groups defined in the \code{GroupOverallResponse}
#'   table (currently 18).
#'
#' @param frequency Integer vector with one or more values of fire frequency in
#'   the preceding fifty year period.
#'
#' @param severity Integer vector with one or more values of severity for the
#'   most recent fire.
#'
#' @param tsf Integer vector with one or more values of time since last fire
#'   (years).
#'
#' @param cross Logical value. If \code{TRUE} (default), a set of fire regimes
#'   will be created for all combinations of the fire component vectors:
#'   frequency, severity and tsf values. In this case, vectors can have
#'   different lengths; e.g. 2 frequency values x 3 severity values x 2 tsf
#'   values will result in 12 fire regimes. If \code{FALSE}, a fire regime
#'   \emph{i} is defined by the \emph{i}th element of each of the three
#'   component vectors. In this case, the vectors must all have the same length.
#'
#' @param interpolation If \code{TRUE} (default), interpolation will be performed
#'   for any combinations of frequency, severity and time since fire values that
#'   do not appear in the \code{GroupOverallResponse} look-up table. If
#'   \code{FALSE}, a warning message will be issued for such combinations and no
#'   data returned for them.
#'
#'
#' @param dx Increment for relative abundance values.
#'
#' @return A data frame suitable for use with ggplot, with columns:
#'   group; frequency; severity; tsf; relabund; density; density_adjusted;
#'   regime.
#'
#' @seealso \code{\link{draw_overall_response}}
#'
#' @export
#'
get_overall_response_data <- function(grp,
                                      frequency, severity, tsf,
                                      cross = TRUE,
                                      interpolation = TRUE,
                                      dx = 0.005) {

  # Get distribution parameters for the group and regime(s).
  # This function also does validity checks on the input arguments.
  #
  dat_beta <- get_zibeta_parameters(grp, frequency, severity, tsf, cross, interpolation)

  # Generate data
  dat_ggs <- lapply(seq_len(nrow(dat_beta)), function(j) {
    suppressWarnings(
      dat <- cbind(dat_beta[j, c("group", "frequency", "severity", "tsf", "pzero")],
                   relabund = seq(0, 1, by=dx))
    )

    if (is.na(dat_beta$shape1[j])) {
      # No beta parameters because of high pzero value
      dat$density <- NA
    } else {
      # Beta parameters are defined
      dat <- dat %>% dplyr::mutate(
        density = dzibeta(relabund, dat_beta$pzero[j], dat_beta$shape1[j], dat_beta$shape2[j]) )
    }

    dat
  })

  # Combine individual data frames
  dat_gg <- do.call(rbind, dat_ggs) %>%
    dplyr::arrange(group, frequency, severity, tsf, relabund)

  # Add a column of fire regime labels
  dat_gg <- dat_gg %>%
    dplyr::mutate(regime = glue::glue("f={frequency}, sev={severity}, tsf={tsf}"),

                  # make label a factor to preserve sensible ordering of values when graphed
                  regime = factor(regime, levels = unique(regime)))

  dat_gg
}


#' Get parameters for ZIBeta distributions that approximate group overall
#' response to specified fire regimes
#'
#' Given a set of integer group IDs, and one or more fire regimes defined as
#' combinations of frequency, severity and time since fire values, this function
#' returns the parameters for zero-inflated beta distributions that approximate
#' the group's overall response to each regime. If a regime is defined in the
#' package look-up table \code{\link{GroupOverallResponse}}, the function simply
#' returns the corresponding distribution parameters. For other fire regimes,
#' the parameters will be derived as the weighted mixture of beta distributions
#' for the most similar reference fire regimes.
#'
#' At the moment, interpolation between fire regimes only involves time since
#' fire, since all integer values for frequency and severity between zero and
#' the maximum values considered during expert-elicitation are represented in
#' the look-up table. Extrapolation beyond the range of any of the three fire
#' regime component variables is not supported.
#'
#' @param grp One or more integer identifiers for group. Each value must be
#'   between 1 and the number of groups defined in the
#'   \code{GroupOverallResponse} table (currently 18).
#'
#' @param frequency Integer vector with one or more values of fire frequency in
#'   the preceding fifty year period.
#'
#' @param severity Integer vector with one or more values of severity for the
#'   most recent fire.
#'
#' @param cross Logical value. If \code{TRUE} (default), a set of fire regimes
#'   will be created for all combinations of the fire component vectors:
#'   frequency, severity and tsf values. In this case, vectors can have
#'   different lengths; e.g. 2 frequency values x 3 severity values x 2 tsf
#'   values will result in 12 fire regimes. If \code{FALSE}, a fire regime
#'   \emph{i} is defined by the \emph{i}th element of each of the three
#'   component vectors. In this case, the vectors must all have the same length.
#'
#' @param tsf Integer vector with one or more values of time since last fire
#'   (years).
#'
#' @param interpolation If \code{TRUE} (default), interpolation will be performed
#'   for any combinations of frequency, severity and time since fire values that
#'   do not appear in the \code{GroupOverallResponse} look-up table. If
#'   \code{FALSE}, a warning message will be issued for such combinations and no
#'   data returned for them.
#'
#' @return A data frame suitable for use with ggplot, with columns:
#'   group; pzero; shape1; shape2.
#'
#' @seealso \code{\link{draw_overall_response}} \code{\link{get_overall_response_data}}
#'
#' @export
get_zibeta_parameters <- function(grp,
                                  frequency, severity, tsf,
                                  cross = TRUE,
                                  interpolation = TRUE) {


  # Initial validity checks
  ok <- .group_is_defined(grp)
  if (any(!ok)) {
    msg <- glue::glue("Group(s) not defined in reference data: \\
                       {paste(grp[!ok], collapse = ', ')}")
    stop(msg)
  }

  ok <- .frequency_is_defined(frequency)
  if (any(!ok)) {
    msg <- glue::glue("Fire frequency value(s) not defined in reference data: \\
                     {paste(frequency[!ok], collapse = ', ')}")
    stop(msg)
  }

  ok <- .severity_is_defined(severity)
  if (any(!ok)) {
    msg <- glue::glue("Fire severity value(s) not defined in reference data: \\
                     {paste(severity[!ok], collapse = ', ')}")
    stop(msg)
  }

  # If interpolation is not being used, check that all TSF values are defined
  if (!interpolation) {
    ok <- .tsf_is_defined(tsf)
    if (any(!ok)) {
      msg <- glue::glue("Time since fire value(s) not defined in reference data
                           and interpolation was not requested: \\
                           {paste(tsf[!ok], collapse = ', ')}")
      stop(msg)
    }

  } else {
    # Interpolation is being used, so check that all TSF values are within
    # the valid range
    MaxTSF <- max(GroupOverallResponse$tsf)
    ok <- tsf >= 0 & tsf <= MaxTSF
    if (any(!ok)) {
      msg <- glue::glue("Time since fire value(s) outside range in reference
                           data (0 - {MaxTSF}): \\
                           {paste(tsf[!ok], collapse = ', ')}")
      stop(msg)
    }
  }

  # Regimes for which data are requested
  if (cross[1]) {
    # Orthogonal combinations
    dat_regimes <- expand.grid(group = grp,
                               frequency = unique(frequency),
                               severity = unique(severity),
                               tsf = unique(tsf))
  } else {
    # Check component value vectors have equal length
    n <- length(frequency)
    if (length(severity) != n || length(tsf) != n) {
      stop("'cross' argument is FALSE but frequency, severity \n",
           "       and tsf vectors do not have the same length")
    }

    # Create a regime for each i-th tuple of component values
    dat_regimes <- data.frame(frequency = frequency,
                              severity = severity,
                              tsf = tsf)

    # Remove any repeated regimes
    dat_regimes <- dplyr::distinct(dat_regimes)

    # Cross with group(s)
    dat_regimes <- dplyr::cross_join(data.frame(group = grp), dat_regimes)
  }

  # Get distribution parameters for regimes that are defined in the
  # look-up table.
  dat_beta_def <- dat_regimes %>%
    # Note: using an inner join drops any undefined regimes.
    dplyr::inner_join(GroupOverallResponse,
                      by = c("group", "frequency", "severity", "tsf"))

  # If any group x regime combinations have high pzero and missing beta parameters,
  # set the shape parameters to 1.0 so graph functions can draw something.
  k <- is.na(dat_beta_def$shape1)
  if (any(k)) {
    dat_beta_def$pzero[k] <- 1.0
    dat_beta_def$shape1[k] <- 1.0
    dat_beta_def$shape2[k] <- 1.0
  }

  if (interpolation) {
    # Interpolate parameters for undefined regimes.
    #
    dat_beta_undef <- dplyr::anti_join(dat_regimes, GroupOverallResponse,
                                       by = c("group", "frequency", "severity", "tsf")) %>%

      dplyr::mutate(pzero = NA_real_, shape1 = NA_real_, shape2 = NA_real_)

    for (j in seq_len(nrow(dat_beta_undef))) {
      target_tsf <- dat_beta_undef$tsf[j]
      dat_enclosing <- .tsf_enclosing_values(target_tsf)

      suppressWarnings(
        zibeta_pars <- cbind(dat_beta_undef[j, c("group", "frequency", "severity")],
                             dat_enclosing) %>%

          dplyr::left_join(GroupOverallResponse,
                           by = c("group", "frequency", "severity", "tsf"))
      )

      if(nrow(zibeta_pars) != 2) {
        stop("Bummer! Problem when trying to interpolate for TSF value ", target_tsf)
      }

      nshape <- sum(!is.na(zibeta_pars$shape1))
      if (nshape == 2) {
        # Both neighbouring fire regimes have beta shape parameters
        dat_beta_undef$pzero[j] <- with(zibeta_pars, sum(wt * pzero))
        dat_beta_undef$shape1[j] <- with(zibeta_pars, sum(wt * shape1))
        dat_beta_undef$shape2[j] <- with(zibeta_pars, sum(wt * shape2))

      } else if (nshape == 1) {
        # One regime missing shape parameters - must be high pzero value.
        # Fall back to sampling the component triangular distributions for
        # both enclosing regimes and deriving a weighted mixture.
        N <- 1e4
        tri <- with(zibeta_pars,
                    get_tri_samples(group[1], frequency[1], severity[1], tsf[1], nsamples=N))

        tri2 <- with(zibeta_pars,
                     get_tri_samples(group[2], frequency[2], severity[2], tsf[2], nsamples=N))

        # Weighted mixture
        b <- runif(N) < zibeta_pars$wt[1]
        tri <- rbind(tri[b,], tri2[!b,])

        # Overall response values calculated using default function (with default args)
        overall <- .FUN_OVERALL(tri)

        pars <- .do_find_zibeta_approximation(overall, pzero_threshold = 0.99)

        dat_beta_undef$pzero[j] <- pars['pzero']
        dat_beta_undef$shape1[j] <- pars['shape1']
        dat_beta_undef$shape2[j] <- pars['shape2']

      } else {
        # Neither enclosing regime has shape parameters.
        # Set pzero to 1.0 and both shape parameters to 1.0 so
        # something can be drawn by graph functions
        # (TODO: better way?)
        dat_beta_undef$pzero[j] <- 1.0
        dat_beta_undef$shape1[j] <- 1.0
        dat_beta_undef$shape2[j] <- 1.0
      }
    }
  } else {  # not interpolating
    dat_beta_undef <- NULL
  }

  # Return parameters
  rbind(dat_beta_def, dat_beta_undef) %>%
    dplyr::arrange(group, frequency, severity, tsf)
}


#' Generate data to draw zero-inflated beta distributions with ggplot
#'
#' This function takes parameters for one or more zero-inflated beta
#' distributions and returns a data frame of density or probability values over
#' the unit interval. You can use it to explore the shape of alternative
#' distributions. If any of the arguments \code{pzero, shape1, shape2} have more
#' than one element, and argument \code{cross == TRUE}, data is generated for
#' all combinations. If \code{cross == FALSE}, the arguments \code{pzero,
#' shape1, shape2} must have the same length and data is generated for
#' \code{pzero[i], shape1[i], shape2[i]} with \code{i in 1:length(pzero)}.
#'
#' @param pzero One or more values for probability of zero.
#'
#' @param shape1 One or more values for the first beta parameter.
#' @param shape2 One or more values for the second beta parameter.
#'
#' @param type A character string specifying the type of y variable; either
#'   'density' (default) or 'probability'. May be abbreviated.
#'
#' @param cross (logical; default \code{TRUE}) Whether to generate data for all
#'   orthogonal combinations of the distribution parameter arguments, or only
#'   for row-wise combinations. See Description for details. Ignored if all
#'   parameter arguments have only one element.
#'
#' @return A data frame with columns: pzero, shape1, shape2, x, y and label. The
#'   label column is a factor that identifies the distribution for each data
#'   row. This can be used to facet or colour ggplot elements.
#'
#' @examples
#' library(ggplot2)
#'
#' # Get data for three distributions that differ in the
#' # probability of zero.
#' #
#' dat <- get_zibeta_ggdata(pzero = c(0, 0.3, 0.6), shape1 = 20, shape2 = 30)
#'
#' ggplot(data = dat) +
#'   geom_line(aes(x=x, y=y, colour=label), linewidth=1.5) +
#'   scale_colour_viridis_d(end=0.8, direction = -1)
#'
#' # Get data for two distributions, fully specifying the parameters
#' # for each. Note cross = FALSE tells the function NOT to generate
#' # data for the eight orthogonal combinations of parameters.
#' #
#' dat <- get_zibeta_ggdata(pzero = c(0, 0.3), shape1 = c(20, 30), shape2 = c(30, 50),
#'                          cross = FALSE)
#'
#' @export
#'
get_zibeta_ggdata <- function(pzero, shape1, shape2,
                              type = c("density", "probability"),
                              cross = TRUE) {
  type <- match.arg(type)

  if (cross) {
    # orthogonal combinations of parameter values
    param_sets <- expand.grid(pzero=pzero, shape1=shape1, shape2=shape2)

  } else {
    # row-wise combinations of parameter values
    if (length(pzero) != length(shape1) ||
        length(pzero) != length(shape2)) {
      stop("Argument cross == FALSE, but arguments pzero, shape1 and shape2 have different lengths")
    }

    param_sets <- data.frame(pzero=pzero, shape1=shape1, shape2=shape2)
  }

  nsets <- nrow(param_sets)

  param_sets$id <- 1:nsets

  fn <- switch(type,
               density = dzibeta,
               probability = pzibeta)

  dat <- lapply(param_sets$id, function(id) {
    ps <- param_sets[id,]

    x <- seq(0, 1, length.out = 1001)
    y <- fn(x, ps$pzero, ps$shape1, ps$shape2)
    data.frame(id = id, pzero = ps$pzero, shape1 = ps$shape1, shape2 = ps$shape2, x, y)
  })

  dat <- do.call(rbind, dat)
  dat$label <- sprintf("zibeta(%g, %g, %g)", dat$pzero, dat$shape1, dat$shape2)
  dat$label <- factor(dat$label, levels = unique(dat$label))

  dat
}

