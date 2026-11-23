# S7 cohort-level regimen distributions and treatment modification.
regimen_codes <- c("PEN","BSP","1GC","2GC","3GC","3GC+PEN","CARB","GLY/LZD","CARB+GLY/LZD","Other","UNK")
regimen_labels <- c("Narrow-spectrum penicillin","Broad-spectrum penicillin",
  "First-generation cephalosporin","Second-generation cephalosporin","Third-generation cephalosporin",
  "Third-generation cephalosporin plus a penicillin-class agent","Carbapenem","Glycopeptide or linezolid",
  "Carbapenem plus glycopeptide or linezolid","Other regimens","Unknown")
stage_labels <- c("Initial antimicrobial regimen","Regimen recorded within 72 h","Regimen recorded during 72 h–7 days")
classic_label <- "Ampicillin or benzylpenicillin plus a third-generation cephalosporin"
normalize_regimen <- function(x) {
  x <- clean(x)
  x[is.na(x) | x=="" | toupper(x) %in% c("UNK","UNKNOWN") | x=="不详"] <- "UNK"
  if (any(!x %in% regimen_codes)) stop("Unknown antimicrobial regimen code.")
  x
}
normalize_change <- function(x) {
  x <- tolower(clean(x))
  x[x %in% c("","unk","unknown","不详")] <- NA_character_
  if (any(!is.na(x) & !x %in% c("change","unchange"))) stop("Unknown treatment modification code.")
  x
}
treatment <- inf[c("HospitalID","Region","Year","t","w")]
for (field in regimen_fields) treatment[[field]] <- normalize_regimen(values(field))
for (field in modification_fields) treatment[[field]] <- normalize_change(values(field))
treatment$Prophylaxis <- categorise(values(prophylaxis_field),"Antifungal prophylaxis")
early_death <- !is.na(inf$D72) & inf$D72==1L
early_withdrawal <- !is.na(inf$C72) & !is.na(inf$D72) & inf$C72==1L & inf$D72==0L
not_assessable <- early_death | early_withdrawal
recorded_late <- clean(values(regimen_fields[3L]))
stopifnot(identical(is.na(recorded_late) | recorded_late=="",not_assessable),
  identical(is.na(treatment[[modification_fields[2L]]]),not_assessable))
treatment$LateEligible <- !not_assessable
treatment$LateRegimen <- ifelse(not_assessable,NA_character_,treatment[[regimen_fields[3L]]])
treatment$LateStatus <- ifelse(early_death,"Not assessable after early death",
  ifelse(early_withdrawal,"Not assessable after early withdrawal",
    ifelse(treatment$LateRegimen=="UNK","Recorded Unknown","Classified")))
late_counts <- do.call(rbind,lapply(regions,function(g) {
  m <- group_mask(treatment,g)
  data.frame(Group=g,CohortN=sum(m),Assessable=sum(m & !not_assessable),
    Classified=sum(m & treatment$LateRegimen!="UNK",na.rm=TRUE),
    RecordedUnknown=sum(m & treatment$LateRegimen=="UNK",na.rm=TRUE),
    NotAssessable=sum(m & not_assessable),EarlyDeath=sum(m & early_death),
    DerivedEarlyWithdrawal=sum(m & early_withdrawal))
}))
flag <- clean(values(classic_field))
stopifnot(all(!is.na(flag) & flag %in% c("","经典方案")))
treatment$ClassicFlag <- as.integer(flag=="经典方案")
stopifnot(all(!treatment$ClassicFlag | treatment[[regimen_fields[1L]]]=="3GC+PEN"))
treatment$ClassicKnownInitial <- ifelse(treatment[[regimen_fields[1L]]]=="UNK",NA_integer_,treatment$ClassicFlag)
classic_counts <- do.call(rbind,lapply(regions,function(g) {
  m <- group_mask(treatment,g)
  x <- treatment$ClassicFlag[m]
  y <- treatment$ClassicKnownInitial[m]
  data.frame(Group=g,Positive=sum(x),Unmarked=sum(x==0L),N=length(x),
    KnownInitialPositive=sum(y,na.rm=TRUE),KnownInitialN=sum(!is.na(y)))
}))
stopifnot(identical(classic_counts$Positive,c(528L,132L,180L,93L,123L)),
  late_counts$Assessable[1L]==1112L,late_counts$RecordedUnknown[1L]==62L,
  late_counts$EarlyDeath[1L]==54L,late_counts$DerivedEarlyWithdrawal[1L]==23L)

s7_titles <- c(stage_labels[1L],"Empiric antibiotic therapy","Modification within 72 h",stage_labels[2L],
  "Modification during 72 h–7 days",stage_labels[3L],"Antifungal prophylaxis")
s7_values <- list(treatment[[regimen_fields[1L]]],treatment$ClassicFlag,
  treatment[[modification_fields[1L]]],treatment[[regimen_fields[2L]]],
  treatment[[modification_fields[2L]]],treatment$LateRegimen,treatment$Prophylaxis)
test_treatment <- function(x, label, analysis) {
  cbind(Characteristic=label,Analysis=analysis,rao_scott(x,treatment,category_first=FALSE))
}
s7_tests <- do.call(rbind,lapply(seq_along(s7_values),function(i) test_treatment(s7_values[[i]],s7_titles[i],"Primary")))
s7_sensitivity <- rbind(
  test_treatment(ifelse(early_death,NA_character_,treatment[[regimen_fields[3L]]]),s7_titles[6L],
    "Exclude early death only; retain early withdrawal blanks as unknown"),
  test_treatment(ifelse(treatment$LateRegimen=="UNK",NA_character_,treatment$LateRegimen),s7_titles[6L],
    "Classified late regimen only"),
  test_treatment(treatment$ClassicKnownInitial,s7_titles[2L],"Known initial regimen only"))
s7_counts <- list()
s7_display <- lapply(seq_along(s7_values),function(i) {
  x <- s7_values[[i]]
  variable <- paste0(s7_titles[i]," (n = ",n_fmt(sum(!is.na(x))),")")
  is_dist <- i %in% c(1L,4L,6L)
  keys <- if(is_dist) c(regimen_codes,if(i==6L) "NOT_ASSESSABLE") else if(i %in% c(2L,7L)) "1" else "change"
  categories <- if(is_dist) c(regimen_labels,if(i==6L) "Not assessable") else if(i==2L) classic_label else variable
  do.call(rbind,lapply(seq_along(keys),function(j) {
    value <- keys[j]
    category <- categories[j]
    r <- do.call(rbind,lapply(regions,function(g) {
      m <- group_mask(treatment,g)
      if(value=="NOT_ASSESSABLE") {
        n <- sum(m)
        events <- sum(m & not_assessable)
      } else {
        n <- sum(m & !is.na(x))
        events <- sum(m & x==value,na.rm=TRUE)
      }
      data.frame(Region=g,N=n,Events=events,Risk=events/n,Display=sprintf("%s (%.1f)",n_fmt(events),100*events/n))
    }))
    s7_counts[[length(s7_counts)+1L]] <<- cbind(Stage=s7_titles[i],Code=value,Category=category,r)
    c(Variable=variable,Category=category,setNames(r$Display,regions),
      "P value"=if(j==1L) fmt_p(s7_tests$P[i]) else "")
  }))
})
S7 <- as.data.frame(do.call(rbind,s7_display),stringsAsFactors=FALSE)
stopifnot(nrow(S7)==38L)
antimicrobial_audit <- list(S7_Counts=do.call(rbind,s7_counts),S7_Pvalue_Audit=rbind(s7_tests,s7_sensitivity),
  S7_Late_Record_Status=late_counts,S7_Classic_Regimen=classic_counts)

