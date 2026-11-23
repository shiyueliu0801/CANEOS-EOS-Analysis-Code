# Hospital-clustered clinical comparisons for Table 4 and S9.
rao_scott <- function(x, data=inf, category_first=TRUE) {
  d <- data[!is.na(x), ]
  d$State <- factor(x[!is.na(x)])
  f <- if (category_first) ~State + Region else ~Region + State
  test <- capture(survey::svychisq(f, design(d), statistic="F"))
  if (is.null(test$Fit)) stop(test$Message)
  z <- test$Fit
  data.frame(N=nrow(d), Hospitals=psu_count(d), F=unname(z$statistic),
    DF1=unname(z$parameter[1L]), DF2=unname(z$parameter[2L]), P=z$p.value,
    Diagnostic=test$Message)
}

plans <- list(
  Region=list(Title="By macro-region",Groups=region_names,Adjuster="Gestational age in weeks"),
  GAGroup=list(Title="By gestational-age group",Groups=c(ga_names[4L],ga_names[1:3]),Adjuster="Macro-region"),
  GA35=list(Title="By 35-week gestational-age threshold",Groups=ga2,Adjuster="Macro-region"))
short_group <- setNames(c(region_names,"Extremely preterm","Very preterm","Moderate to late preterm","Term",ga2),
  c(region_names,ga_names,ga2))
columns <- c("Outcome / contrast","Observed, n/N (%)","Adjusted risk, % (95% CI)",
  "Crude RR (95% CI)","Adjusted RR (95% CI)","Overall P (crude / adjusted)")
heading <- function(x, observed="", test="") c(x,observed,"","","",test)
rr_text <- function(effect, stability, i) {
  if (effect$Group[i]==effect$Reference[i]) return("Reference")
  if (!stability$RRDisplayed[i]) return("—")
  fmt_interval(effect$RR[i],effect$RRLower[i],effect$RRUpper[i])
}
compare_outcome <- function(id, dimension) {
  plan <- plans[[dimension]]
  crude <- fit_binary(inf,dimension,plan$Groups,id)
  adjusted <- if(dimension=="Region") fit_binary(inf,dimension,plan$Groups,id,adjust=TRUE) else
    fit_ga_adjusted(id,dimension,plan$Groups)
  a <- marginal(crude)
  b <- marginal(adjusted)
  stopifnot(identical(a$Group,b$Group), all(abs(a$Risk-a$ObservedRisk)<1e-10),
    all(abs(b$RR-b$Risk/b$Risk[1L])<1e-10))
  cr <- sparse_rr(crude,a,support_by_group(crude$data,plan$Groups))
  ad <- sparse_rr(adjusted,b,support_by_group(adjusted$data,plan$Groups))
  labels_group <- unname(short_group[a$Group])
  labels_group[1L] <- paste0(labels_group[1L]," (ref.)")
  groups <- do.call(rbind,lapply(seq_len(nrow(a)),function(i) c(
    labels_group[i],risk_text(a$Events[i],a$GroupN[i]),
    fmt_interval(100*b$Risk[i],100*b$RiskLower[i],100*b$RiskUpper[i],n=1),
    rr_text(a,cr,i),rr_text(b,ad,i),"")))
  p <- paste(fmt_p(crude$audit$P),fmt_p(adjusted$audit$P),sep=" / ")
  contrasts <- cbind(Code=id,Outcome=labels[[id]],Dimension=dimension,a[c("Group","Reference")],
    CrudeEvents=a$Events,CrudeN=a$GroupN,CrudeRisk=a$Risk,CrudeRR=a$RR,
    CrudeLower=a$RRLower,CrudeUpper=a$RRUpper,CrudeDisplayed=cr$RRDisplayed,CrudeReason=cr$RRReason,
    CrudeRD=a$RD_pp,CrudeRDLower=a$RDLower_pp,CrudeRDUpper=a$RDUpper_pp,
    AdjustedEvents=b$Events,AdjustedGroupN=b$GroupN,StandardisationN=b$StandardN,
    AdjustedRisk=b$Risk,AdjustedRiskLower=b$RiskLower,AdjustedRiskUpper=b$RiskUpper,
    AdjustedRR=b$RR,AdjustedLower=b$RRLower,AdjustedUpper=b$RRUpper,
    AdjustedDisplayed=ad$RRDisplayed,AdjustedReason=ad$RRReason,
    AdjustedRD=b$RD_pp,AdjustedRDLower=b$RDLower_pp,AdjustedRDUpper=b$RDUpper_pp)
  tests <- rbind(cbind(Model="Crude",Covariate="None",crude$audit),
    cbind(Model="Adjusted",Covariate=plan$Adjuster,adjusted$audit))
  list(display=rbind(heading(plan$Title,risk_text(sum(a$Events),sum(a$GroupN)),p),groups),
    contrasts=contrasts,tests=tests)
}
clinical_blocks <- setNames(lapply(outcome_codes,function(id) {
  setNames(lapply(names(plans),function(dimension) compare_outcome(id,dimension)),names(plans))
}),outcome_codes)
S9 <- do.call(rbind,lapply(outcome_codes,function(id) {
  do.call(rbind,c(list(heading(labels[[id]])),lapply(clinical_blocks[[id]],"[[","display")))
}))
S9 <- setNames(as.data.frame(S9,stringsAsFactors=FALSE),columns)
s9_contrasts <- do.call(rbind,lapply(clinical_blocks,function(block) do.call(rbind,lapply(block,"[[","contrasts"))))
s9_tests <- do.call(rbind,lapply(clinical_blocks,function(block) do.call(rbind,lapply(block,"[[","tests"))))
stopifnot(nrow(S9)==182L,nrow(s9_tests)==78L,all(s9_tests$Converged),
  sum(s9_contrasts$Group!=s9_contrasts$Reference & !s9_contrasts$CrudeDisplayed)==4L,
  sum(s9_contrasts$Group!=s9_contrasts$Reference & !s9_contrasts$AdjustedDisplayed)==4L)

# Table 4 reports live cohort counts. Its IVH, NEC and meningitis tests
# distinguish confirmed, suspected and negative diagnoses.
table4_counts <- list()
table4_intervals <- list()
table4_tests <- list()
format_table4 <- function(id, section, x, title, test, multistate=FALSE, ci=FALSE) {
  n <- sum(!is.na(x))
  if (multistate) {
    states <- c(NA_integer_,1L,2L)
    display <- c(paste0(title," (n = ",n_fmt(n),")"),"Confirmed","Suspected")
  } else {
    states <- 1L
    display <- paste0(title," (n = ",n_fmt(n),")")
  }
  table4_tests[[id]] <<- cbind(Code=id,test)
  rows <- lapply(seq_along(states),function(i) {
    state <- states[i]
    r <- do.call(rbind,lapply(regions,function(g) {
      mask <- group_mask(inf,g) & !is.na(x)
      denom <- sum(mask)
      events <- if(is.na(state)) NA_integer_ else sum(x[mask]==state)
      pr <- events/denom
      lo <- hi <- NA_real_
      if (ci) {
        d <- inf[mask, ]
        d$y <- as.integer(x[mask]==state)
        dd <- design(d)
        fit <- survey::svymean(~y,dd)
        df <- survey::degf(dd)
        estimate <- estimate_ci(qlogis(pr),1/(pr*(1-pr)),vcov(fit),df,"logit")
        check <- unname(confint(survey::svyciprop(~y,dd,method="xlogit",df=df)))
        stopifnot(max(abs(estimate[2:3]-check))<1e-10)
        lo <- estimate[2L]
        hi <- estimate[3L]
        table4_intervals[[length(table4_intervals)+1L]] <<- data.frame(Code=id,Region=g,Events=events,
          ValidN=denom,Risk=pr,Lower=unname(lo),Upper=unname(hi),Hospitals=psu_count(d),DF=df)
      }
      value <- if (is.na(state)) "" else if (ci) {
        sprintf("%s / %s (%.1f%% [%.1f–%.1f])",n_fmt(events),n_fmt(denom),100*pr,100*lo,100*hi)
      } else sprintf("%s (%.1f)",n_fmt(events),100*pr)
      data.frame(Region=g,N=denom,Events=events,Risk=pr,Display=value)
    }))
    table4_counts[[length(table4_counts)+1L]] <<- cbind(Code=id,Section=section,
      Characteristic=display[i],State=state,r)
    c(Section=section,Characteristic=display[i],setNames(r$Display,regions),
      "P value"=if(i==1L) fmt_p(test$P) else "")
  })
  do.call(rbind,rows)
}
make_management <- function(fields, section) {
  lapply(names(fields), function(id) {
    x <- categorise(values(fields[[id]]),id)
    format_table4(id,section,x,id,rao_scott(x),id=="Meningitis")
  })
}
make_clinical <- function(codes, section) {
  lapply(codes,function(id) {
    three <- id %in% c("IVH","NEC")
    x <- if(three) inf[[paste0(id,"_state")]] else inf[[id]]
    test <- if(three) rao_scott(x) else clinical_blocks[[id]]$Region$tests[1L,c("N","Hospitals","F","NumDF","ResidualDF","P","Diagnostic")]
    names(test)[names(test)=="NumDF"] <- "DF1"
    names(test)[names(test)=="ResidualDF"] <- "DF2"
    format_table4(id,section,x,if(three) id else labels[[id]],test,three,id %in% disposition_codes)
  })
}
Table4 <- as.data.frame(do.call(rbind,c(
  make_management(respiratory_fields,"Respiratory support within 72 h of birth"),
  make_management(management_fields,"Early clinical features and management"),
  make_clinical(setdiff(outcome_codes,disposition_codes),"Major complications"),
  make_clinical(disposition_codes,"In-hospital mortality and treatment withdrawal")
)),stringsAsFactors=FALSE)
stopifnot(nrow(Table4)==30L,length(table4_intervals)==35L)
clinical_audit <- list(Table4_Counts=do.call(rbind,table4_counts),
  Table4_CI=do.call(rbind,table4_intervals),Table4_Tests=do.call(rbind,table4_tests),
  S9_Contrasts=s9_contrasts,S9_Overall_Tests=s9_tests)

