# Set random seed for reproducibility
set.seed(20230701)

# Define site codes
sites <- c(
  "SITE-WE-D9F845",
  "SITE-WE-A99B46",
  "SITE-WE-14FA00",
  "SITE-WE-C34692",
  "SITE-WE-B5D72A",
  "SITE-WE-4AE03A",
  "SITE-WE-A8248B",
  "SITE-WE-FE3A6F",
  "SITE-WE-2EEAA7",
  "SITE-WE-2CD088",
  "SITE-WE-C4FB6F",
  "SITE-WE-81A924",
  "SITE-WE-90291F",
  "SITE-WE-9800A7",
  "SITE-WE-2F0F70",
  "SITE-WE-B6D4A7",
  "SITE-WE-E3CDF7",
  "SITE-WE-37E0CD",
  "SITE-WE-0684E3",
  "SITE-WE-07792A",
  "SITE-WE-2891E3",
  "SITE-WE-CD0000",
  "SITE-WE-238F94",
  "SITE-WE-1D3BA3"
)

# Perform simple random ordering
random_order <- sample(sites, length(sites), replace = FALSE)

# Create a data frame with ranking
result <- data.frame(
  Rank = 1:length(sites),
  Site_Code = random_order
)

# Print result
print(result)
