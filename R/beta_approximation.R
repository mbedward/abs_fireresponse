#' Fit a beta distribution to the overall fire response of a group
#'
#' Note: This is a private function and is mainly intended for generating the
#' package data frame \code{\link{GroupOverallResponse}}
#'
#' @param the_group Integer group number: a single value between 1 and the
#'   number of groups defined in \code{\link{GroupExpertData}} (currently 18).
#'
#' @param frequency Fire frequency in the preceding fifty year period.
#'
#' @param severity Severity of most recent fire.
#'
#' @param tsf Time since last fire (years).
#'
#' @param weights Weights to apply to the component variables when calculating
#'   the overall response as a weighted product. Values must be positive and
#'   will be scaled to sum to 1.0. The default is to weight all components
#'   equally.
#'
#' @param nsamples Number of random samples to generate from each component
#'   triangular distribution.
#'
#' @return A two-element vector of beta parameters
#'
#' @keywords internal
#'
find_beta_approximation <- function(the_group, frequency, severity, tsf, weights = c(1,1,1), nsamples = 1e4) {

  stopifnot(length(the_group) == 1,
            the_group >= 1,
            the_group <= max(GroupExpertData$group))

  dat_samples <- get_tri_samples(the_group,
                                 frequency = frequency, severity = severity, tsf = tsf,
                                 weights = weights, nsamples = nsamples)

  # Fit a beta distribution
  init_pars = c(mean(dat_samples$overall), 1.0)
  suppressWarnings({
    o <- optim(par = init_pars, fn = .fn_beta_ll, gr=NULL, x = dat_samples$overall,
               method = "L-BFGS-B",
               lower = c(0.01, 0.1), upper = c(0.99, 100))
    fitted_pars <- o$par
  })

  shape1 <- fitted_pars[2] * fitted_pars[1]
  shape2 <- fitted_pars[2] * (1 - fitted_pars[1])

  c(shape1=shape1, shape2=shape2)
}


#' Draw random samples from triangular distributions for a given group and fire regime
#'
#' Note: This is a private function and is mainly intended for generating the
#' package data frame \code{\link{GroupOverallResponse}}
#'
#' @param the_group Integer group number: a single value between 1 and the
#'   number of groups defined in \code{\link{GroupExpertData}} (currently 18).
#'
#' @param frequency Fire frequency in the preceding fifty year period.
#'
#' @param severity Severity of most recent fire.
#'
#' @param tsf Time since last fire (years).
#'
#' @param weights Weights to apply to the component variables when calculating
#'   the overall response as a weighted product. All values must be positive.
#'   The default is to weight all components equally. Order of elements is
#'   frequency, severity, tsf.
#'
#' @param nsamples Number of random samples to generate from each component
#'   triangular distribution.
#'
#' @return A data frame with \code{nsamples} rows
#'
#' @keywords internal
#'
get_tri_samples <- function(the_group, frequency, severity, tsf, weights = c(1,1,1), nsamples = 1e4) {
  stopifnot(length(the_group) == 1,
            the_group >= 1,
            the_group <= max(GroupExpertData$group))

  stopifnot(length(weights) == 3 && all(weights > 0))

  dat <- GroupExpertData %>%
    dplyr::filter(group == the_group & (
      (type == "frequency" & value == frequency) |
      (type == "severity" & value == severity) |
      (type == "tsf" & value == tsf)) )

  if (nrow(dat) != 3) {
    msg <- glue::glue("No data for frequency={frequency}, severity={severity}, tsf={tsf}")
    stop(msg)
  }

  # Function to apply small adjustments to tri-parameters as required to ensure a
  # valid distribution to sample
  fn_make_valid <- function(lwr, upr, mode, eps = 1e-3) {
    res <- c(lwr, upr, mode)
    if (isTRUE(all.equal(lwr, upr))) {
      res <- res + c(-eps, 0, -eps/2)
    }
    res
  }

  dat_samples <- dat %>%
    dplyr::group_by(type) %>%
    dplyr::do({
      tri_pars <- fn_make_valid(.$ra_lwr, .$ra_upr, .$ra_mode)
      data.frame(rep = 1:nsamples, r = triangulr::rtri(nsamples, min=tri_pars[1], max=tri_pars[2], mode=tri_pars[3]))
    }) %>%
    dplyr::ungroup() %>%

    # Reshape to wide format with a column for each fire variable sample vector
    tidyr::pivot_wider(names_from = type, values_from = r) %>%
    dplyr::select(-rep) %>%

    # Weighted product
    dplyr::rowwise() %>%
    dplyr::mutate(overall = prod(c(frequency, severity, tsf) * weights)) %>%
    dplyr::ungroup()

  dat_samples
}


#### Non-exported helper functions

# Beta distribution parameterized via mean and dispersion
.dbeta2 <- function(x, mu, phi, ...) dbeta(x, phi*mu, phi*(1-mu), ...)


# Calculate minus the log-likelihood given a vector of two beta parameters and a
# vector of data values
.fn_beta_ll <- function(pars, x) {
  mu <- pars[1]
  phi <- pars[2]
  -sum(.dbeta2(x, mu, phi, log=TRUE))
}

