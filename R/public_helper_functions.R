#' Get a vector of all integer ID values defined for functional groups in the expert data
#'
#' This is a convenience function to get a vector of all functional group ID
#' values in ascending order.
#'
#' @export
get_group_ids <- function() { sort(unique(GroupOverallResponse$group)) }
