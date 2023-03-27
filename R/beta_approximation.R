#' Combine sampled values from fire component distributions into overall response values
#'
#' This function can be used with \code{\link{find_zibeta_approximation}} to
#' calculate overall response values based on vectors of values sampled from the
#' fire component distributions. It combines the component values as a weighted
#' product raised to a power \code{lambda} (default = 1/3).
#'
#' @param xcomponents A three-column matrix or data frame of sample values for
#'   fire regime components. Column order must be: frequency, severity, tsf.
#'
#' @param weights A vector of three, positive weight values to apply to the
#'   component variables when calculating the overall response as a weighted
#'   product. Order of elements is frequency, severity, tsf. Input weights will
#'   be scaled so that they sum to 3.0. Zero values are allowed (might be useful
#'   for exploratory purposes). The default is to weight all components equally.
#'
#' @param lambda A single numeric value (default = 1/3) to map the weighted
#'   product to the overall response value as:
#'   \deqn{y_{overall} = (\sum{x_i wt_i})^\lambda }
#'
#' @return A vector of overall response values.
#'
.FUN_overall_product = function(xcomponents, weights = c(1,1,1), lambda = 1/3) {
  stopifnot(ncol(xcomponents) == 3)
  stopifnot(length(weights) == 3)
  stopifnot(all(weights >= 0))
  stopifnot(sum(weights) > 0)

  if (length(lambda) > 1) {
    warning("lambda should be length 1; only using first value")
    lambda <- lambda[1]
  }

  weights <- weights * 3 / sum(weights)

  apply(xcomponents, MARGIN = 1, function(xi) prod(xi * weights)^lambda)
}


#' Combine sampled values from fire component distributions into overall response values
#'
#' This function can be used with \code{\link{find_zibeta_approximation}} to
#' calculate overall response values based on vectors of values sampled from the
#' fire component distributions. Overall response values are calculated using a
#' heuristic adopted by DPE researchers in which overall response is set to zero
#' if any component sample value is less than 0.05; otherwise the overall response
#' is calculated as the weighted sum of component values.
#'
#' @param xcomponents A three-column matrix or data frame of sample values for
#'   fire regime components. Column order must be: frequency, severity, tsf.
#'
#' @param weights A vector of three, positive weight values to apply to the
#'   component variables when calculating the overall response as a weighted
#'   sum. Order of elements is frequency, severity, tsf. Zero values are allowed
#'   (might be useful for exploratory purposes). The default is to weight all
#'   components equally. The provided values will be scaled to sum to 1.0.
#'
#' @param zero_threshold Threshold value for individual component samples. If
#'   any component has a lower value the overall response is set to zero.
#'
#' @return A vector of overall response values.
#'
.FUN_overall_sum_threshold = function(xcomponents, weights = c(1,1,1), zero_threshold = 0.05) {
  stopifnot(ncol(xcomponents) == 3)
  stopifnot(all(weights >= 0))
  stopifnot(sum(weights) > 0)

  weights <- weights / sum(weights)

  # DPE threshold for individual components
  zero_adj <- apply(xcomponents, MARGIN = 1, function(xi) if(min(xi) < zero_threshold) 0 else 1)

  # Overall response
  (as.matrix(xcomponents) %*% weights) * zero_adj
}



#' Derive overall fire response and approximating zero-inflated beta distribution
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
#' @param nsamples Number of random samples to generate from each component
#'   triangular distribution. Default value is 10,000. Minimum allowable value
#'   is 1000.
#'
#' @param FUN A function that will be applied to the matrix of component sample
#'   values to calculate the corresponding overall response values. Default is
#'   the function \link{.FUN_overall_sum_threshold}.
#'
#' @return A numeric vector with three named elements:
#' \describe{
#'   \item{pzero}{probability of zero relative abundance value}
#'   \item{shape1}{first beta parameter for non-zero values (NA if pzero is 1.0)}
#'   \item{shape2}{second beta parameter for non-zero values (NA if pzero is 1.0)}
#' }
#'
find_zibeta_approximation <- function(the_group,
                                      frequency, severity, tsf,
                                      nsamples = 1e4,
                                      FUN = .FUN_overall_sum_threshold) {

  stopifnot(length(the_group) == 1)
  stopifnot(the_group %in% GroupExpertData$group)

  dat_samples <- get_tri_samples(the_group,
                                 frequency = frequency, severity = severity, tsf = tsf,
                                 nsamples = nsamples)

  # Calculate overall response value for each sample row
  overall <- FUN(dat_samples)

  # Proportion of zero sample values
  isz <- overall == 0
  pzero <- mean(isz)

  # Fit a beta distribution to the non-zero values, if there are enough.
  shape1 <- NA_real_
  shape2 <- NA_real_
  if (pzero < 0.99) {
    init_pars = c(mean(overall[!isz]), 1.0)
    suppressWarnings({
      o <- optim(par = init_pars, fn = .fn_beta_ll, gr=NULL, x = overall[!isz],
                 method = "L-BFGS-B",
                 lower = c(0.01, 0.1), upper = c(0.99, 100))
      fitted_pars <- o$par
    })

    shape1 <- fitted_pars[2] * fitted_pars[1]
    shape2 <- fitted_pars[2] * (1 - fitted_pars[1])
  }

  c(pzero = pzero, shape1=shape1, shape2=shape2)
}


#' Sample component triangular distributions for a group and fire regime
#'
#' Draws random samples from each of the triangular distributions representing
#' the bounded expert estimates of the effect of each fire regime component
#' (frequency, severity and time since fire) on the expected relative abundance
#' of a given group. This function is mainly intended for generating the package
#' data frame \code{\link{GroupOverallResponse}} but can also be used for other
#' purposes.
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
#' @param nsamples Number of random samples to generate from each component
#'   triangular distribution.
#'
#' @param seed Seed value to initialize the pseudo-random number generator
#'   (default = 42). This allows reproducible results between sessions. Set to
#'   \code{NULL} or \code{NA} to allow random samples from a given set of
#'   distributions to vary between sessions.
#'
#' @return A data frame with \code{nsamples} rows and a column for each fire
#'   component: frequency, severity and tsf. The group and seed value (if
#'   defined) are added as attributes of the returned data frame.
#'
get_tri_samples <- function(the_group,
                            frequency, severity, tsf,
                            nsamples = 1e4,
                            seed = 42) {

  stopifnot(length(the_group) == 1,
            the_group >= 1,
            the_group <= max(GroupExpertData$group))

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
      if (lwr > 0.0 & upr < 1.0) {
        # squashed but away from 0 or 1, so nudge lwr and upr apart
        res <- res + c(-min(eps, lwr), min(eps, upr), 0)

      } else if (lwr < eps) {
        # squashed at 0, so nudge upr and mode upwards
        res <- res + c(0, eps, 0)

      } else {
        # squashed at 1, so nudge lwr and mode downwards
        res <- res + c(-eps, 0, 0)
      }
    }
    res
  }

  if (!is.null(seed) && !is.na(seed)) set.seed(seed)
  dat_samples <- dat %>%
    dplyr::group_by(type) %>%
    dplyr::do({
      tri_pars <- fn_make_valid(.$ra_lwr, .$ra_upr, .$ra_mode)
      data.frame(rep = 1:nsamples, r = triangulr::rtri(nsamples, min=tri_pars[1], max=tri_pars[2], mode=tri_pars[3]))
    }) %>%
    dplyr::ungroup() %>%

    # Reshape to wide format with a column for each fire variable's sample vector
    tidyr::pivot_wider(names_from = type, values_from = r) %>%
    dplyr::select(-rep)

  # Record group and seed value for this set of samples
  attr(dat_samples, "group") <- the_group
  attr(dat_samples, "seed") <- seed

  dat_samples
}


#### Non-exported helper functions

# Beta distribution density function parameterized via mean and dispersion
.dbeta2 <- function(x, mu, phi, ...) dbeta(x, phi*mu, phi*(1-mu), ...)


# Calculate negative summed log-likelihood given a vector of two beta parameters and a
# vector of data values
.fn_beta_ll <- function(pars, x) {
  mu <- pars[1]
  phi <- pars[2]
  -sum(.dbeta2(x, mu, phi, log=TRUE))
}

