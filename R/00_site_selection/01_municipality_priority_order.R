# Set random seed for reproducibility
set.seed(20230701)

# Original site codes
site_codes <- c(
  "SITE-MU-CD9420",
  "SITE-MU-DA1844",
  "SITE-MU-587BE5",
  "SITE-MU-0FDB54",
  "SITE-MU-9B6174",
  "SITE-MU-9805D9",
  "SITE-MU-47BBAF",
  "SITE-MU-3641E9",
  "SITE-MU-F57831",
  "SITE-MU-0FE18F",
  "SITE-MU-24D6CA",
  "SITE-MU-CAB78C",
  "SITE-MU-D27642",
  "SITE-MU-E09CC0",
  "SITE-MU-1BEAA5",
  "SITE-MU-F1F984"
)

# Perform simple random ordering
random_order <- sample(site_codes, length(site_codes), replace = FALSE)

# Create result data frame with ranks starting from 1
result <- data.frame(Rank = 1:length(random_order), Site = random_order)

# Print the result
print(result)
