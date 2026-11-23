# Explicit computational outputs for the final S1 and S2 Tables.
# Inputs and weights are controlled; no study aggregate is embedded here.
s1 <- reshape(hospital_year[c("HospitalID", "Region", "Year", "Isolates")],
              idvar = c("HospitalID", "Region"), timevar = "Year", direction = "wide")
s1 <- s1[match(site_map$HospitalID, s1$HospitalID),
         c("Region", "HospitalID", paste0("Isolates.", rev(year_set)))]
names(s1)[-(1:2)] <- as.character(rev(year_set))
write_result(s1, "S1_Table_Site_Year_Isolates")
cal <- static_cfg$get_calibration_weights()
cal$Sample_proportion_pct <- 100 * cal$Isolates / ave(cal$Isolates, cal$Year, FUN = sum)
cal$Population_proportion_pct <- 100 * cal$Births / ave(cal$Births, cal$Year, FUN = sum)
write_result(cal, "S2_Table_Calibration")
