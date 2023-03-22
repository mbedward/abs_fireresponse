# Private (non-exported) helper functions


# Test if one or more group IDs are present in the GroupOverallResponse table.
#
.group_is_defined <- function(grps) {
  sapply(grps, `%in%`, unique(GroupOverallResponse$group))
}


# Test if one or more integer fire frequency values are present in the
# GroupOverallResponse table.
#
.frequency_is_defined <- function(freqs) {
  sapply(freqs, `%in%`, unique(GroupOverallResponse$frequency))
}


# Test if one or more integer fire severity values are present in the
# GroupOverallResponse table.
#
.severity_is_defined <- function(sevs) {
  sapply(sevs, `%in%`, unique(GroupOverallResponse$severity))
}


# Test if one or more integer time since fire values are present in the
# GroupOverallResponse table.
#
.tsf_is_defined <- function(tsfs) {
  sapply(tsfs, `%in%`, unique(GroupOverallResponse$tsf))
}


# For one or more given values of TSF, find the nearest enclosing values of TSF
# defined in the GroupOverallResponse table, and calculate a weight value for
# each of the enclosing values to use when interpolating response values.
#
# Returns a data.frame with two rows for each desired TSF value, and columns:
#   target (desired TSF)
#   tsf (left/right enclosing reference values)
#   wt (left/right weights to use when combining beta densities)
#
.tsf_enclosing_values <- function(tsfs) {
  TSF <- sort(unique(GroupOverallResponse$tsf))

  res <- lapply(tsfs, function(tsf) {
    d <- TSF - tsf
    ilwr <- which(d == max(d[d <= 0]))
    iupr <- which(d == min(d[d >= 0]))

    vals <- TSF[c(ilwr, iupr)]

    width <- diff(vals)

    if (width > 0) {
      weights <- 1.0 - abs(vals - tsf) / width
    } else {
      # degenerate case where tsf value was defined
      weights <- c(1.0, 0.0)
    }
    data.frame(target = tsf, tsf = vals, wt = weights)
  })

  do.call(rbind, res)
}

