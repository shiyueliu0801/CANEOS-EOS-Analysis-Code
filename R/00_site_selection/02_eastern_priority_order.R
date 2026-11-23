# Set random seed for reproducibility
set.seed(20230701)

# Original site codes
site_codes <- c(
  "SITE-EA-AB9AD4",
  "SITE-EA-A6ED4E",
  "SITE-EA-10FEE6",
  "SITE-EA-A27AAB",
  "SITE-EA-E1570C",
  "SITE-EA-71D14F",
  "SITE-EA-E9A03C",
  "SITE-EA-54B3CF",
  "SITE-EA-2CF172",
  "SITE-EA-436024",
  "SITE-EA-E2BEC6",
  "SITE-EA-55A96F",
  "SITE-EA-D1E285",
  "SITE-EA-73B697",
  "SITE-EA-798185",
  "SITE-EA-DFCAAC",
  "SITE-EA-9CB875",
  "SITE-EA-0CE08E",
  "SITE-EA-BD64F2",
  "SITE-EA-C5369B",
  "SITE-EA-1C81AF",
  "SITE-EA-33DA12",
  "SITE-EA-CD170D",
  "SITE-EA-A5C0BB",
  "SITE-EA-EB2845",
  "SITE-EA-1470AB",
  "SITE-EA-13F8F8",
  "SITE-EA-DBD24B"
)

# Generate random order
random_order <- sample(site_codes)

# Create the result data frame
result <- data.frame(
  Rank = 1:length(random_order),
  SiteCode = random_order
)

# Print the result
print(result)
