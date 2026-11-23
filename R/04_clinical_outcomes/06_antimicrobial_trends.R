# S16 annual distributions and hospital-clustered treatment trends.

carb_codes <- c("CARB","CARB+GLY/LZD")
model_specs <- list()
add_trend <- function(stage, indicator, code, source, eligible, x, family="Main") {
  stopifnot(length(eligible)==nrow(treatment),!anyNA(eligible),all(is.na(x[!eligible])))
  model_specs[[length(model_specs)+1L]] <<- list(Stage=stage,Indicator=indicator,Code=code,
    SourceField=source,Family=family,StageEligible=eligible,X=x)
}
full <- rep(TRUE,nrow(treatment))
for (i in seq_along(regimen_fields)) {
  codes <- if(i==3L) treatment$LateRegimen else treatment[[regimen_fields[i]]]
  eligible <- if(i==3L) treatment$LateEligible else full
  known <- eligible & !is.na(codes) & codes!="UNK"
  unambiguous <- known & codes!="Other"
  for (j in seq_len(10L)) {
    code <- regimen_codes[j]
    reference <- if(code=="Other") known else unambiguous
    add_trend(stage_labels[i],regimen_labels[j],code,regimen_fields[i],eligible,
      ifelse(reference,as.integer(codes==code),NA_integer_))
  }
  add_trend(stage_labels[i],"Either carbapenem category (supplementary)","ANY_CARB",regimen_fields[i],eligible,
    ifelse(unambiguous,as.integer(codes %in% carb_codes),NA_integer_),"Supplementary")
  if(i==1L) add_trend(stage_labels[i],classic_label,"CLASSIC",classic_field,full,treatment$ClassicFlag)
}
add_trend("Treatment modification",s7_titles[3L],"CHANGE_72",modification_fields[1L],full,
  as.integer(treatment[[modification_fields[1L]]]=="change"))
add_trend("Treatment modification",s7_titles[5L],"CHANGE_7",modification_fields[2L],treatment$LateEligible,
  as.integer(treatment[[modification_fields[2L]]]=="change"))
add_trend("Antifungal prophylaxis",s7_titles[7L],"PROPHYLAXIS",prophylaxis_field,full,
  treatment$Prophylaxis)
annual_masks <- expand.grid(Group=regions,Year=year_set,stringsAsFactors=FALSE)
annual_distribution <- do.call(rbind,lapply(seq_len(nrow(annual_masks)),function(j) {
  yr <- annual_masks$Year[j]
  g <- annual_masks$Group[j]
  d <- treatment[treatment$Year==yr & group_mask(treatment,g), ]
  do.call(rbind,lapply(seq_along(regimen_fields),function(i) {
    x <- if(i==3L) d$LateRegimen else d[[regimen_fields[i]]]
    cohort <- nrow(d)
    eligible <- sum(!is.na(x))
    unknown <- sum(x=="UNK",na.rm=TRUE)
    counts <- as.integer(table(factor(x,levels=regimen_codes)))
    stopifnot(sum(counts)==eligible,eligible==if(i==3L) sum(d$LateEligible) else cohort)
    extra <- i==3L
    events <- c(counts,if(extra) cohort-eligible)
    denom <- c(rep(eligible,length(counts)),if(extra) cohort)
    data.frame(Year=yr,Group=g,Stage=stage_labels[i],SourceField=regimen_fields[i],
      Code=c(regimen_codes,if(extra) "NOT_ASSESSABLE"),
      Regimen=c(regimen_labels,if(extra) "Not assessable"),
      Events=events,Denominator=denom,Proportion=ifelse(denom>0,events/denom,NA_real_),
      CohortN=cohort,AssessableN=eligible,RecordedUnknown=unknown,NotAssessable=cohort-eligible,
      EarlyDeath=if(extra) sum(d$LateStatus=="Not assessable after early death") else 0L,
      DerivedEarlyWithdrawal=if(extra) sum(d$LateStatus=="Not assessable after early withdrawal") else 0L)
  }))
}))
for (g in regions) for (i in seq_along(stage_labels)) {
  x <- annual_distribution[annual_distribution$Group==g & annual_distribution$Stage==stage_labels[i], ]
  n <- late_counts[match(g,late_counts$Group), ]
  stopifnot(sum(x$AssessableN[x$Code=="PEN"])==if(i==3L) n$Assessable else n$CohortN)
  if(i==3L) {
    stopifnot(sum(x$Events[x$Code=="UNK"])==n$RecordedUnknown,
              sum(x$Events[x$Code=="NOT_ASSESSABLE"])==n$NotAssessable)
  }
}
annual_indicators <- do.call(rbind,lapply(seq_along(model_specs),function(i) {
  spec <- model_specs[[i]]
  do.call(rbind,lapply(seq_len(nrow(annual_masks)),function(j) {
    yr <- annual_masks$Year[j]
    g <- annual_masks$Group[j]
    mask <- treatment$Year==yr & group_mask(treatment,g)
    x <- spec$X[mask]
    n <- sum(!is.na(x))
    e <- sum(x==1L,na.rm=TRUE)
    data.frame(Index=i,Year=yr,Group=g,Stage=spec$Stage,Indicator=spec$Indicator,Code=spec$Code,
      CohortN=sum(mask),StageEligibleN=sum(mask & spec$StageEligible),ModelN=n,Events=e,
      Proportion=if(n) e/n else NA_real_)
  }))
}))

trend_diagnostics <- function(d, a) {
  d <- d[!is.na(d$x), ]
  tab <- table(factor(d$Year,levels=year_set),factor(d$x,levels=0:1))
  eligible_string <- paste(rowSums(tab),collapse="/")
  positive_string <- paste(tab[,"1"],collapse="/")
  sep <- NA
  if(length(unique(d$x))==2L) {
    r0 <- range(d$t[d$x==0L])
    r1 <- range(d$t[d$x==1L])
    sep <- r0[2L]<=r1[1L] || r1[2L]<=r0[1L]
  }
  ordinary <- NA_real_
  max_shift <- NA_real_
  failures <- NA_integer_
  if(isFALSE(sep)) {
    get_glm <- function(z) {
      if(length(unique(z$x))!=2L) return(NULL)
      f <- capture(glm(x~t, data=z, family=binomial(), control=controls))
      if (is.null(f$Fit) || !isTRUE(f$Fit$converged) || f$Fit$boundary || nzchar(f$Message) ||
          any(!is.finite(coef(f$Fit)))) return(NULL)
      f$Fit
    }
    naive <- get_glm(d)
    if(!is.null(naive)) ordinary <- sqrt(vcov(naive)["t","t"])
    if(a$Status=="Estimated") {
      stopifnot(!is.null(naive),abs(coef(naive)["t"]-a$Beta)<1e-7,a$CovarianceMaxError<1e-6)
      leave_one <- vapply(unique(d$HospitalID),function(h) {
        f <- get_glm(d[d$HospitalID!=h, ])
        if(is.null(f)) NA_real_ else unname(coef(f)["t"])
      },0)
      failures <- sum(!is.finite(leave_one))
      if(any(is.finite(leave_one))) max_shift <- max(abs(leave_one[is.finite(leave_one)]-a$Beta))
    }
  }
  base <- data.frame(OrdinarySE=ordinary,SandwichToOrdinarySE=a$SE/ordinary,
    SeparationByYear=sep,MaxHospitalDeletionShift=max_shift,HospitalDeletionFailures=failures,
    EligibleByYear=eligible_string,EventsByYear=positive_string)
  event_years <- sum(tab[,"1"]>0L)
  nonevent_years <- sum(tab[,"0"]>0L)
  ratio <- base$SandwichToOrdinarySE
  near_precision <- is.finite(ratio) && ratio<0.5 &&
    (isTRUE(base$HospitalDeletionFailures>0L) || isTRUE(base$MaxHospitalDeletionShift>a$SE))
  problems <- c(if(a$Status!="Estimated") a$Diagnostic,
    if(min(a$Events,a$NonEvents)<10L) "Fewer than 10 events or non-events.",
    if(min(a$EventHospitals,a$NonEventHospitals)<4L) "Fewer than four event or non-event hospitals.",
    if(min(event_years,nonevent_years)<3L) "Fewer than three event or non-event years.",
    if(isTRUE(base$SeparationByYear)) "Complete or quasi-complete separation.",
    if(near_precision) "Sandwich SE below half the ordinary SE with unstable hospital deletion.",
    if(isTRUE(a$Status=="Estimated" && (a$FiveUpper>1e5 || a$FiveLower<1e-5))) "Extremely wide five-year confidence interval.")
  reporting <- if(a$Status=="Estimated" && !length(problems)) "Report" else "Withheld"
  cbind(base,EventYears=event_years,NonEventYears=nonevent_years,
    Reporting=reporting,ReportingReason=paste(unique(problems[nzchar(problems)]),collapse=" "))
}

fit_s16 <- function(i,g) {
  s <- model_specs[[i]]
  mask <- group_mask(treatment,g)
  d <- treatment[mask, ]
  d$x <- s$X[mask]
  a <- fit_trend(d)
  diagnostics <- trend_diagnostics(d,a)
  counts <- annual_indicators[annual_indicators$Index==i & annual_indicators$Group==g, ]
  stopifnot(sum(counts$ModelN)==a$N,sum(counts$Events)==a$Events)
  cbind(Index=i,Stage=s$Stage,Indicator=s$Indicator,Code=s$Code,Subgroup=g,SourceField=s$SourceField,
    CohortN=nrow(d),StageEligibleN=sum(s$StageEligible[mask]),FamilyType=s$Family,
    a,diagnostics)
}
s16_models <- do.call(rbind,lapply(seq_along(model_specs),function(i)
  do.call(rbind,lapply(regions,function(g) fit_s16(i,g)))))
validate_trends(s16_models,"S16 OR, variance and time scaling")
s16_models$Family <- paste(s16_models$FamilyType,ifelse(s16_models$Subgroup=="Overall","pooled","macro-regions"))
s16_models$BH_q <- NA_real_
s16_models$FamilyN <- NA_integer_
for(f in unique(s16_models$Family)) {
  family <- s16_models$Family==f
  use <- family & s16_models$Reporting=="Report"
  s16_models$FamilyN[family] <- sum(family)
  s16_models$BH_q[use] <- p.adjust(s16_models$P[use],method="BH",n=sum(family))
}
final_trends <- s16_models[c("Index","Stage","Indicator","Code","Subgroup","StageEligibleN","N","Events",
  "Reporting")]
display <- t(vapply(seq_len(nrow(s16_models)),function(i) {
  a <- s16_models[i, ]
  if(a$Reporting=="Report") trend_cells(a) else rep("—",5L)
},character(5L)))
final_trends <- cbind(final_trends,setNames(as.data.frame(display),c("Beta","SE","Annual","Five","P")))

classic_sensitivity <- do.call(rbind,lapply(regions,function(g) {
  d <- treatment[group_mask(treatment,g), ]
  d$x <- d$ClassicKnownInitial
  a <- fit_trend(d)
  z <- trend_diagnostics(d,a)
  cbind(Subgroup=g,Restriction="Initial regimen known; positive flag versus unmarked",a,z)
}))
stopifnot(sum(classic_sensitivity$N[-1L])==classic_sensitivity$N[1L],
          all(classic_sensitivity$Events==classic_counts$Positive))

# The submission table contains the 33 displayed regimen/indicator blocks.
# Other regimens and the recorded classic flag remain in the analysis audit.
submission_s16 <- final_trends[!final_trends$Code %in% c("Other", "CLASSIC"), ]
S16 <- data.frame(
  Stage = submission_s16$Stage,
  `Regimen / indicator` = sub(" (supplementary)", "", submission_s16$Indicator, fixed=TRUE),
  Subgroup = submission_s16$Subgroup,
  `Events/N` = paste0(n_fmt(submission_s16$Events), "/", n_fmt(submission_s16$N)),
  `β per year` = submission_s16$Beta, SE = submission_s16$SE,
  `Annual OR (95% CI)` = submission_s16$Annual, `5-year OR (95% CI)` = submission_s16$Five,
  `P value` = submission_s16$P, check.names=FALSE)
stopifnot(length(model_specs)==37L,nrow(S16)==165L,nrow(annual_distribution)==1700L)
s16_results <- list(S16_Antimicrobial_Trends=S16, S16_All_Trends=final_trends,
  S16_Annual_Distributions=annual_distribution, S16_Annual_Indicators=annual_indicators,
  S16_Model_Audit=s16_models, S16_Classic_Sensitivity=classic_sensitivity)
