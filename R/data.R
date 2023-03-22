#' Expert-elicited group responses to fire regime variables
#'
#' A data frame with, for each functional group, bounded estimates for the
#' expected relative abundance (on the unit interval) given discrete values for
#' time since last fire (years); severity of last fire; and frequency of fire in
#' the preceding fifty year period.
#'
#' @format A data frame with 666 rows and 6 columns:
#' \describe{
#' \item{group}{Integer identifier for the functional group (1-18)}
#' \item{type}{Type of fire variable, coded as a factor with levels: frequency; severity; tsf}
#' \item{value}{Value of the given fire variable (e.g. years for TSF)}
#' \item{ra_lwr}{Estimate of lowest expected relative abundance}
#' \item{ra_mode}{Estimate of most likely value for relative abundance}
#' \item{ra_max}{Estimate of highest expected relative abundance}
#' }
"GroupExpertData"


#' Overall group responses to specified fire regimes
#'
#' A data frame with parameters for a zero-inflated beta distribution that best
#' approximates the overall response of a functional group to a fire regime, as
#' defined by: frequency of fire in the preceding fifty year period; severity of
#' last fire; and time since last fire (years).
#'
#' @format A data frame with 1683 rows and 7 columns:
#' \describe{
#' \item{group}{Integer identifier for the functional group (1-18)}
#' \item{frequency}{Number of fires in the preceding fifty years}
#' \item{severity}{Severity of the most recent fire}
#' \item{tsf}{Time since last fire (years)}
#' \item{pzero}{Probability of zero relative abundance}
#' \item{shape1}{First parameter of the approximating beta distribution}
#' \item{shape2}{Second parameter of the approximating beta distribution}
#' }
"GroupOverallResponse"
