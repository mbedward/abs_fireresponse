#' Generate data to draw zero-inflated beta distributions with ggplot
#'
#' This function takes parameters for one or more zero-inflated beta
#' distributions and returns a data frame of density or probability values over
#' the unit interval. You can use it to explore the shape of alternative
#' distributions. If any of the arguments \code{pzero, shape1, shape2} have more
#' than one element, data is generated for all combinations.
#'
#' @param pzero One or more values for probability of zero.
#'
#' @param shape1 One or more values for the first beta parameter.
#' @param shape2 One or more values for the second beta parameter.
#'
#' @param type A character string specifying the type of y variable; either
#'   'density' (default) or 'probability'. May be abbreviated.
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
#' dat <- get_ggdata_zibeta(pzero = c(0, 0.3, 0.6), shape1 = 20, shape2 = 30)
#'
#' ggplot(data = dat) +
#'   geom_line(aes(x=x, y=y, colour=label), linewidth=1.5) +
#'   scale_colour_viridis_d(end=0.8, direction = -1)
#'
#' @export
#'
get_ggdata_zibeta <- function(pzero, shape1, shape2,
                              type = c("density", "probability")) {
  type <- match.arg(type)

  param_sets <- expand.grid(pzero=pzero, shape1=shape1, shape2=shape2)
  nsets <- nrow(param_sets)

  o <- rev(order(param_sets$pzero))
  param_sets <- param_sets[o, ]

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

