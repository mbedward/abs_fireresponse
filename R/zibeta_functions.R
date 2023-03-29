#' Zero-inflated Beta Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the Zero-inflated Beta distribution with parameters \code{pzero},
#' \code{shape1} and \code{shape2}.
#'
#' @param x Vector of quantiles.
#'
#' @param p Vector of probabilities.
#'
#' @param n Number of observations. If \code{length(n) > 1}, the length is taken
#'   to be the number required.
#'
#' @param pzero Probability of zero value.
#'
#' @param shape1,shape2 Non-negative parameters of the beta distribution for non-zero values.
#'
#' @param log,log.p If \code{TRUE}, log density / probability values are returned.
#'
#' @param lower.tail If \code{TRUE} (default), probabilities are \eqn{P[X \le x]},
#'   otherwise \eqn{P[X \gt x]}.
#'
#' @return \code{dbeta} returns the density, \code{pbeta} returns the probability,
#'   \code{qbeta} returns the quantile, and \code{rbeta} generates random values.
#'
#' @export
#'
dzibeta <- function(x, pzero, shape1, shape2, log=FALSE) {
  # Create four column matrix of p values and parameters, allowing default recycling
  # to happen for any variable(s) as per stats::pbeta()
  pars <- suppressWarnings(cbind(x, pzero, shape1, shape2))
  X <- 1; PZERO <- 2; SHAPE1 <- 3; SHAPE2 <- 4;

  d <- apply(pars, MARGIN = 1, function(xpars) {
    if(xpars[X] == 0) xpars[PZERO]
    else (1-xpars[PZERO]) * dbeta(xpars[X], xpars[SHAPE1], xpars[SHAPE2])
  })

  if (log[1]) d <- log(d)

  d
}


#' @rdname dzibeta
#'
#' @export
#'
pzibeta <- function(x, pzero, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {
  # Create four column matrix of p values and parameters, allowing default recycling
  # to happen for any variable(s) as per stats::pbeta()
  pars <- suppressWarnings(cbind(x, pzero, shape1, shape2))
  X <- 1; PZERO <- 2; SHAPE1 <- 3; SHAPE2 <- 4;

  p <- apply(pars, MARGIN = 1, function(xpars) {
    if(xpars[X] == 0) xpars[PZERO]
    else xpars[PZERO] + (1-xpars[PZERO]) * pbeta(xpars[X], xpars[SHAPE1], xpars[SHAPE2])
  })

  if (!lower.tail[1]) p <- 1-p
  if (log.p[1]) p <- log(p)

  p
}

#' @rdname dzibeta
#'
#' @export
#'
qzibeta <- function(p, pzero, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {
  # Create four column matrix of p values and parameters, allowing default recycling
  # to happen for any variable(s) as per stats::qbeta()
  pars <- suppressWarnings(cbind(p, pzero, shape1, shape2))
  P <- 1; PZERO <- 2; SHAPE1 <- 3; SHAPE2 <- 4;

  if (!lower.tail[1]) pars[,P] <- 1 - pars[,P]
  if (log.p[1]) pars[,P] <- exp(pars[,P])

  apply(pars, MARGIN = 1, function(xpars) {
    if(xpars[P] <= xpars[PZERO]) 0
    else qbeta(p = (xpars[P]-xpars[PZERO])/(1-xpars[PZERO]),
               shape1 = xpars[SHAPE1], shape2 = xpars[SHAPE2],
               lower.tail=TRUE, log.p=FALSE)
  })
}


#' @rdname dzibeta
#'
#' @export
#'
rzibeta <- function(n, pzero, shape1, shape2) {
  # Force all parameters to be length(1) to avoid confusion
  stopifnot(length(n) == 1)
  stopifnot(length(pzero) == 1)
  stopifnot(length(shape1) == 1)
  stopifnot(length(shape2) == 1)

  r <- rep(0, n)
  nonzero <- runif(n) >= pzero
  r[nonzero] <- rbeta(sum(nonzero), shape1, shape2)
  r
}
