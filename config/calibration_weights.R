# ==============================================================================
# Poststratification Calibration Weights (Macro-Region x Year, 2015-2024)
# ==============================================================================
# Public aggregate parameters published in S2 Table of the CANEOS manuscript.
# Contains no patient-level or hospital-level records.
#
# Columns:
#   Year     : Calendar year (2015-2024)
#   Region   : Macro-region stratum (Municipality, Eastern, Central, Western)
#   Isolates : Observed isolate count in stratum-year cell (n_s)
#   Births   : National live births in 10,000s (n_p, China Statistical Yearbook)
#   Weight   : Within-year poststratification weight (wg)
# ==============================================================================

get_calibration_weights <- function() {
  data.frame(
    Year     = c(2015L, 2015L, 2015L, 2015L, 2016L, 2016L, 2016L, 2016L, 2017L, 2017L, 2017L, 2017L, 2018L, 2018L, 2018L, 2018L, 2019L, 2019L, 2019L, 2019L, 2020L, 2020L, 2020L, 2020L, 2021L, 2021L, 2021L, 2021L, 2022L, 2022L, 2022L, 2022L, 2023L, 2023L, 2023L, 2023L, 2024L, 2024L, 2024L, 2024L),
    Region   = c("Municipality", "Central", "Eastern", "Western", "Municipality", "Central", "Eastern", "Western", "Municipality", "Central", "Eastern", "Western", "Municipality", "Central", "Eastern", "Western", "Municipality", "Central", "Eastern", "Western", "Municipality", "Central", "Eastern", "Western", "Municipality", "Central", "Eastern", "Western", "Municipality", "Central", "Eastern", "Western", "Municipality", "Central", "Eastern", "Western", "Municipality", "Central", "Eastern", "Western"),
    Isolates = c(35L, 28L, 19L, 33L, 28L, 25L, 26L, 35L, 31L, 36L, 22L, 29L, 24L, 33L, 20L, 38L, 27L, 31L, 25L, 30L, 32L, 34L, 32L, 26L, 26L, 32L, 34L, 36L, 26L, 34L, 33L, 30L, 40L, 32L, 44L, 26L, 30L, 19L, 26L, 29L),
    Births   = c(53.47021, 678.11706, 486.13831, 370.26323, 44.3044, 649.0884, 474.03367, 369.77909, 53.29531, 673.69799, 556.24337, 384.31995, 50.63874, 700.57102, 586.36131, 400.49384, 45.10669, 638.53322, 520.98941, 376.76617, 44.47085, 600.63454, 500.75216, 361.63591, 36.099, 477.37387, 379.55462, 304.39193, 32.80068, 412.7499, 341.33246, 270.2861, 29.62378, 367.31106, 305.81676, 247.99419, 28.22791, 339.64799, 288.92376, 239.48804),
    Weight   = c(0.110635435, 1.753868555, 1.852919912, 0.8125443, 0.117344228, 1.925469944, 1.352099795, 0.78351481, 0.121654651, 1.324232277, 1.789136951, 0.937770032, 0.139605812, 1.404656318, 1.939845579, 0.697339222, 0.119375483, 1.471837115, 1.489110187, 0.897405224, 0.114311968, 1.453107132, 1.287179461, 1.144101222, 0.148417631, 1.594675557, 1.193325819, 0.90384572, 0.146781099, 1.412435726, 1.203439564, 1.048245705, 0.110612553, 1.71438343, 1.038084101, 1.424597217, 0.109180071, 2.074251203, 1.289424188, 0.958234349),
    stringsAsFactors = FALSE
  )
}
