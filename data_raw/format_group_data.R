# Format the raw expert data table as provided by DPE ('all_species.csv')
#
library(dplyr, warn.conflicts = FALSE)
library(stringr)

GroupExpertData <- local ({
  path <- here::here("data_raw/all_species.csv")

  read.csv(path) %>%
    select(group = Group,
           type = Type,
           value = Xaxis,
           ra_lwr = RelativeAbundLow,
           ra_mode = RelativeAbundAve,
           ra_upr = RelativeAbundUp) %>%

    mutate(group = as.integer(str_extract(group, "\\d+")),
           type = factor( tolower(type) ) ) %>%

    arrange(group, type, value)
})

usethis::use_data(GroupExpertData, overwrite = TRUE)

# Clean up
rm(GroupExpertData)
