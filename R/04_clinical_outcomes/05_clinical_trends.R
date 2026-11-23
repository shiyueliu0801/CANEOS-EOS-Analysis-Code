# S17 clinical outcome trends. Calendar year is continuous.
# Trend intervals use the accepted hospital-clustered finite-df inference.
trend_outcomes <- c(names(source_fields),"Withdrawal")
fit_clinical_trends <- function(codes, suspected_as_event=FALSE) {
  do.call(rbind,lapply(codes,function(id) {
    do.call(rbind,lapply(c(regions,ga_names),function(g) {
      mask <- if(g %in% regions) group_mask(inf,g) else !is.na(inf$GAGroup) & inf$GAGroup==g
      d <- inf[mask, ]
      d$x <- if(suspected_as_event) as.integer(d[[paste0(id,"_state")]]>0L) else d[[id]]
      model <- fit_trend(d)
      diagnostic <- trend_support(d,model)
      cbind(Code=id,Outcome=if(suspected_as_event) paste0(id," (confirmed or suspected)") else labels[[id]],
        Subgroup=g,Role=if(suspected_as_event) "Suspected as event" else "Primary",model,diagnostic)
    }))
  }))
}
trend_display <- function(models) {
  values <- t(vapply(seq_len(nrow(models)),function(i) {
    row <- models[i, ]
    if(row$Reporting=="Report") trend_cells(row) else rep("—",5L)
  },character(5L)))
  names <- c("β per year","SE","Annual OR (95% CI)","5-year OR (95% CI)","P value")
  cbind(models[c("Outcome","Subgroup")],setNames(as.data.frame(values),names))
}
clinical_trend_primary <- fit_clinical_trends(trend_outcomes)
clinical_trend_sensitivity <- fit_clinical_trends(c("IVH","NEC"),suspected_as_event=TRUE)
clinical_trend_models <- rbind(clinical_trend_primary,clinical_trend_sensitivity)
validate_trends(clinical_trend_models,"Clinical trends")
S17 <- trend_display(clinical_trend_primary)
stopifnot(nrow(S17)==117L,nrow(clinical_trend_sensitivity)==18L)
clinical_trend_audit <- list(S17_Model_Audit=clinical_trend_models,
  S17_Suspected_Sensitivity=trend_display(clinical_trend_sensitivity))

