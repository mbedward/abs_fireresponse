#' Density function for the zero-inflated beta distribution
#'
#' @param x Vector of quantiles.
#' @param pzero Probability of zero value.
#' @param shape1,shape2 Non-negative parameters of the beta distribution for non-zero values.
#' @param log If \code{TRUE}, log density values are returned.
#'
#' @return Vector of density values.
#'
#' @export
#'
dzibeta <- function(x, pzero, shape1, shape2, log=FALSE) {
  d <- ifelse(x == 0, pzero, (1 - pzero) * dbeta(x, shape1, shape2))
  if (log) d <- log(d)
  d
}


#' Distribution function for the zero-inflated beta distribution
#'
#' @param x Vector of quantiles.
#' @param pzero Probability of zero value.
#' @param shape1,shape2 Non-negative parameters of the beta distribution for
#'   non-zero values.
#' @param lower.tail If \code{TRUE} (default), probabilities are \deqn{P[X \le x]},
#'   otherwise \deqn{P[X \gt x]}.
#' @param log.p If \code{TRUE}, log probability values are returned.
#'
#' @return Vector of probability values.
#'
#' @export
#'
pzibeta <- function(x, pzero, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {
  p <- ifelse(x == 0, pzero, pzero + (1 - pzero) * pbeta(x, shape1, shape2))

  if (!lower.tail) p <- 1-p
  if (log.p) p <- log(p)

  p
}

#' Distribution function for the zero-inflated beta distribution
#'
#' @param x Vector of probabilities.
#' @param pzero Probability of zero value.
#' @param shape1,shape2 Non-negative parameters of the beta distribution for
#'   non-zero values.
#' @param lower.tail If \code{TRUE} (default), probabilities are \deqn{P[X \le x]},
#'   otherwise \deqn{P[X \gt x]}.
#' @param log.p If \code{TRUE}, log probability values are returned.
#'
#' @return Vector of probability values.
#'
#' @export
#'
qzibeta <- function(p, pzero, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- exp(p)

  ifelse(p <= pzero,
         qbeta((p-pzero)/(1-pzero), shape1, shape2, lower.tail=TRUE, log.p=FALSE)
         )
}


#' Draw random values from a zero-inflated beta distribution
#'
#' @param n Number of observations. If \code{length(n) > 1}, the length is taken
#'   to be the number required.
#' @param shape1,shape2 Non-negative parameters of the beta distribution for
#'   non-zero values.
#'
#' @return Vector of random values.
#'
rzibeta <- function(n, pzero, shape1, shape2) {
  r <- rep(0, n)
  nonzero <- runif(n) >= pzero
  r[nonzero] <- rbeta(sum(nonzero), shape1, shape2)
  r
}
