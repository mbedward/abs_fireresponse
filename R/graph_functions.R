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
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap label_both scale_x_continuous labs
#'
#' @export
#'
draw_overall_response <- function(the_group, frequency, severity, tsf,
                                  adjusted_density = TRUE,
                                  linewidth = 1,
                                  draw_pzero = TRUE) {

  dat <- get_overall_response_data(the_group, frequency, severity, tsf)

  gg_dens <- .do_draw_overall_response(dat,
                                       adjusted_density = adjusted_density,
                                       linewidth = linewidth)
  if (draw_pzero) {
    gg_zero <- .do_draw_prob_zero(dat)

    patchwork::wrap_plots(gg_zero, gg_dens, ncol = 1, heights = c(2,8))

  } else {
    gg_dens
  }
}


# Private helper function to draw overall response densities
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous labs facet_wrap label_both
.do_draw_overall_response <- function(dat,
                                      adjusted_density = TRUE,
                                      linewidth = 1) {

  ngroups <- dplyr::n_distinct(dat$group)
  nregimes <- dplyr::n_distinct(dat$regime)

  if (adjusted_density) {
    gg <- ggplot(data = dat, aes(x = relabund, y = density_adjusted))
  } else {
    gg <- ggplot(data = dat, aes(x = relabund, y = density))
  }

  gg <- gg +
    if (nregimes > 1) {
      geom_line(aes(colour = regime), linewidth=linewidth)
    } else {
      geom_line(linewidth=linewidth)
    }

  if (ngroups > 1) {
    gg <- gg + facet_wrap(~group, labeller = label_both, scales = "free_y")
  }

  gg <- gg +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    labs(x = "Relative abundance", y = ifelse(adjusted_density[1], "Adjusted density", "Density"))

  gg
}


# Private helper function to draw prob zero graph
#' @importFrom ggplot2 ggplot aes geom_point geom_segment scale_x_continuous labs facet_wrap label_both theme
.do_draw_prob_zero <- function(dat) {
  dat_zero <- dat %>%
    dplyr::group_by(group, regime) %>%
    dplyr::summarize(pzero = dplyr::first(pzero))

  ngroups <- dplyr::n_distinct(dat$group)

  gg <- ggplot(data = dat_zero, aes(y = regime)) +
    geom_segment(x = 0, aes(xend = pzero, yend = regime, colour = regime), linewidth = 2) +
    geom_point(aes(x = pzero, colour = regime), size = 4) +
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
                                      interpolation = TRUE,
                                      dx = 0.005) {


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
  dat_regimes <- expand.grid(group = grp,
                             frequency = frequency,
                             severity = severity,
                             tsf = tsf)

  # Get zero-inflated beta parameters for regimes that are present in the
  # look-up table.
  dat_beta <- dat_regimes %>%
    # Note: using an inner join drops any undefined regimes (TSF not present in look-up)
    dplyr::inner_join(GroupOverallResponse, by = c("group", "frequency", "severity", "tsf"))

  # Generate data for defined regimes
  dats_def <- lapply(seq_len(nrow(dat_beta)), function(j) {
    suppressWarnings(
      dat <- cbind(dat_beta[j, c("group", "frequency", "severity", "tsf", "pzero")],
                   relabund = seq(0, 1, by=dx))
    )

    if (is.na(dat_beta$shape1[j])) {
      # No beta parameters because of high pzero value
      dat$density <- NA
      dat$adjusted_density <- NA
    } else {
      # Beta parameters are defined
      dat <- dat %>% dplyr::mutate(
        density = dbeta(relabund, dat_beta$shape1[j], dat_beta$shape2[j]),
        density_adjusted = density * (1.0 - pzero))
    }

    dat
  })

  # Generate data for undefined regimes (ie. TSF values needing interpolation).
  # NOTE: we are assuming that the checks above ensured interpolation was requested.
  #
  dat_undef_regimes <- dplyr::anti_join(dat_regimes, GroupOverallResponse,
                                        by = c("group", "frequency", "severity", "tsf"))

  dats_interp <- lapply(seq_len(nrow(dat_undef_regimes)), function(j) {
    target_tsf <- dat_undef_regimes$tsf[j]
    dat_enclosing <- .tsf_enclosing_values(target_tsf)

    suppressWarnings(
    beta_pars <- cbind(dat_undef_regimes[j, c("group", "frequency", "severity")],
                       dat_enclosing) %>%

      dplyr::left_join(GroupOverallResponse, by = c("group", "frequency", "severity", "tsf"))
    )

    if(nrow(beta_pars) != 2) stop("Bummer: error when trying to interpolate")

    # Interpolated value for prob zero
    pzero_interp <- with(beta_pars, sum(wt * pzero))

    suppressWarnings(
      res <- cbind(beta_pars[1, c("group", "frequency", "severity")],
                   tsf = target_tsf,
                   pzero = pzero_interp,
                   relabund = seq(0, 1, by=dx))
    )

    if (anyNA(beta_pars$shape1)) {
      # One or both sets of beta parameters missing because of
      # high pzero value, so can't interpolate
      res$density <- NA_real_
      res$density_adjusted <- NA_real_

    } else {
      shape1_interp <- with(beta_pars, sum(wt * shape1))
      shape2_interp <- with(beta_pars, sum(wt * shape2))

      res <- res %>%
        dplyr::mutate(density = dbeta(relabund, shape1_interp, shape2_interp),
                      density_adjusted = density * (1.0 - pzero))
    }

    res
  })

  # Combine individual data frames
  dat_gg <- do.call(rbind, c(dats_def, dats_interp)) %>%
    dplyr::arrange(group, frequency, severity, tsf, relabund)

  # Add a column of fire regime labels
  dat_gg <- dat_gg %>%
    dplyr::mutate(regime = glue::glue("f={frequency}, sev={severity}, tsf={tsf}"),

                  # make label a factor to preserve sensible ordering of values when graphed
                  regime = factor(regime, levels = unique(regime)))

  dat_gg
}


