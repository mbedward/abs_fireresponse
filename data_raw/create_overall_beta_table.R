# Generate the GroupOverallResponse table with a row for each group x fire
# regime and columns for group, frequency, severity, tsf, and the parameters of
# the approximating beta distribution (shape1, shape2).
#
# Note: this script assumes that the GroupExpertData table has already been created and resides
# in the 'data' folder as GroupExpertData.rda.

library(dplyr, warn.conflicts = FALSE)

library(doFuture)
registerDoFuture()

initls <- ls(all.names = TRUE)

# On UOW processing machine
NCORES <- 40

# On laptop
# NCORES <- 4


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

  dat_params <- foreach(i = 1:nrow(dat_combns)) %dopar% {
    vals <- unlist(dat_combns[i,])
    beta_pars <- find_beta_approximation(the_group,
                                         frequency = vals['frequency'],
                                         severity = vals['severity'],
                                         tsf = vals['tsf'])

    data.frame(group = the_group,
               frequency = vals['frequency'],
               severity = vals['severity'],
               tsf = vals['tsf'],
               shape1 = beta_pars['shape1'],
               shape2 = beta_pars['shape2'])
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

