# Set random seed for reproducibility
set.seed(20230701)

# Original site codes
site_codes <- c(
  "SITE-CE-04AA06",
  "SITE-CE-FF4C80",
  "SITE-CE-14B960",
  "SITE-CE-4BAD51",
  "SITE-CE-7B82A4",
  "SITE-CE-CA842C",
  "SITE-CE-C4F5C5",
  "SITE-CE-AD8E68",
  "SITE-CE-54A424",
  "SITE-CE-D2B886",
  "SITE-CE-681A7F",
  "SITE-CE-CD91CC",
  "SITE-CE-63A60F",
  "SITE-CE-8C4B06",
  "SITE-CE-669A59",
  "SITE-CE-58B273",
  "SITE-CE-BB7986",
  "SITE-CE-8829AE",
  "SITE-CE-2CDF0E",
  "SITE-CE-3A672B",
  "SITE-CE-04C9C3",
  "SITE-CE-B502FD",
  "SITE-CE-4B4599",
  "SITE-CE-A30AEF",
  "SITE-CE-75F4C3",
  "SITE-CE-0439AD"
)

# Perform simple random ordering
random_order <- sample(site_codes, length(site_codes), replace = FALSE)

# Create result data frame with ranks starting from 1
result <- data.frame(Rank = 1:length(random_order), Site = random_order)

# Print the result
print(result)
