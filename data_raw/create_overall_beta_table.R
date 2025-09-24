# Generate the GroupOverallResponse table with a row for each group x fire
# regime and columns for group, frequency, severity, tsf, and the parameters of
# the approximating beta distribution (shape1, shape2).
#
# Note: this script assumes that the GroupExpertData table has already been created and resides
# in the 'data' folder as GroupExpertData.rda.

library(dplyr, warn.conflicts = FALSE)
library(here)

library(doFuture)

initls <- ls(all.names = TRUE)


### IMPORTANT - set the number of parallel cores to use here

# On a high-end processing machine with 50 cores it will take about 15 seconds
# to process each group.
# NCORES <- 50

# On laptop or similar use whatever you can find. E.g. On a machine with 8 cores
# and a 2.8GHz processor it will take about 90 seconds to process each group.
NCORES <- 8


source(here::here("R/beta_approximation.R"))

load(here::here("data/GroupExpertData.rda"))
stopifnot(exists("GroupExpertData"))


get_unique_vals <- function(the_type) {
  GroupExpertData %>%
    filter(type == the_type) %>%
    distinct(value) %>%
    pull(value)
}

# Fire regimes defined in expert data
dat_combns <- expand.grid(frequency = get_unique_vals("frequency"),
                          severity = get_unique_vals("severity"),
                          tsf = get_unique_vals("tsf"))


GroupOverallResponse <- lapply(sort(unique(GroupExpertData$group)), function(the_group) {

  msg <- glue::glue("Group {the_group}")
  message(msg)

  plan(multisession, workers = NCORES)

  dat_params <- foreach(i = 1:nrow(dat_combns), .options.future = list(seed = TRUE)) %dofuture% {
    vals <- unlist(dat_combns[i,])

    zib_pars <- find_zoabeta_approximation(the_group,
                                          frequency = vals['frequency'],
                                          severity = vals['severity'],
                                          tsf = vals['tsf'])

    data.frame(group = the_group,
               frequency = vals['frequency'],
               severity = vals['severity'],
               tsf = vals['tsf'],
               pzero = zib_pars['pzero'],
               shape1 = zib_pars['shape1'],
               shape2 = zib_pars['shape2'])
  }

  plan(sequential)

  bind_rows(dat_params)
})

GroupOverallResponse <- bind_rows(GroupOverallResponse) %>%
  arrange(group, frequency, severity, tsf)

rownames(GroupOverallResponse) <- NULL

usethis::use_data(GroupOverallResponse, overwrite = TRUE)

# Clean up all objects that were not here when we started
rm(list = setdiff(ls(all.names = TRUE), initls))

