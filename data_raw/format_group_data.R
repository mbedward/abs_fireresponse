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

# Create a lookup table of fire component labels
FireComponentLookup <- data.frame(
  type = c("frequency", "severity", "tsf"),
  short_label = c("Frequency", "Severity", "TSF"),
  long_label = c("Frequency (50 years)", "Severity", "Time since fire")
)

# Check labels are in sync with data
stopifnot(setequal(FireComponentLookup$type, GroupExpertData$type))

usethis::use_data(FireComponentLookup, overwrite = TRUE)

# Clean up
rm(GroupExpertData)
rm(FireComponentLookup)
