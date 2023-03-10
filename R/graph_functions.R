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
                                  linewidth = 1) {
  ngroups <- length(the_group)
  nregimes <- prod(length(frequency), length(severity), length(tsf))

  dat <- get_overall_response_data(the_group, frequency, severity, tsf)

  gg <- ggplot(data = dat, aes(x = ra, y = dens))

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
    labs(x = "Relative abundance", y = "Density")

  gg
}


#' Prepare data to graph group overall response to fire regime(s)
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
#' @return A data frame.
#'
#' @export
#'
get_overall_response_data <- function(the_group, frequency, severity, tsf,
                                      dx = 0.005) {

  dat_combns <- expand.grid(group = the_group,
                            frequency = frequency,
                            severity = severity,
                            tsf = tsf)

  dat_beta <- dat_combns %>%
    dplyr::left_join(GroupOverallResponse, by = c("group", "frequency", "severity", "tsf"))

  if (nrow(dat_beta) == 0) {
    warning("None of the requested combinations of group and fire regime are defined")

    # Return empty data frame
    dat_beta

  } else {
    # Check all combinations were found
    x <- dplyr::distinct(dat_beta, group, frequency, severity, tsf)
    if (nrow(x) != nrow(dat_combns)) {
      warning("One or more requested combinations of group and fire regime are not defined")
    }

    dats <- lapply(seq_len(nrow(dat_beta)), function(i) {
      suppressWarnings(
        dat <- cbind(dat_beta[i, c("group", "frequency", "severity", "tsf")],
                     ra = seq(0, 1, by=dx))
      )

      dat$dens <- dbeta(dat$ra, dat_beta$shape1[i], dat_beta$shape2[i])

      dat
    })

    # Combine individual data frames
    dat_gg <- do.call(rbind, dats)

    # Add a column of fire regime labels
    dat_gg <- dat_gg %>%
      dplyr::mutate(regime = glue::glue("f={frequency}, sev={severity}, tsf={tsf}"),

                    # make label a factor to preserve sensible ordering of values when graphed
                    regime = factor(regime, levels = unique(regime)))

    dat_gg
  }
}


