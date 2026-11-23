groups <- list(Overall="Overall",Region=region,GA=ga,Delivery=del)
is_level <- function(field,level) if(field=="Overall") rep(TRUE,nrow(iso)) else iso[[field]] %in% level
by_groups <- function(fun,spec=groups) do.call(rbind,lapply(names(spec),function(field) {
  do.call(rbind,lapply(spec[[field]],function(level) fun(field,level,is_level(field,level))))
}))
table2 <- by_groups(function(field,level,m) do.call(rbind,lapply(c("Gram-positive","Gram-negative",org_names),function(o) {
  y <- if(o %in% c("Gram-positive","Gram-negative")) iso$Gram==o else iso$Organism==o
  data.frame(Block=field,Level=level,Pathogen=o,estimate(domain(make_y(design,as.integer(y)),m),paste("Table2",o,level)))
})))


overall_test <- function(d,x,y,label,weighted=TRUE) {
  f <- reformulate(c(y,x))
  if(!weighted) d$prob[] <- 1
  z <- attempt(label,survey::svychisq(f,d,statistic="F"))
  v <- if(is.null(z)) rep(NA_real_,4) else c(as.numeric(z$statistic),as.numeric(z$parameter),z$p.value)
  data.frame(Statistic=v[1],Numerator_DF=v[2],Denominator_DF=v[3],P=v[4],Design_DF=survey::degf(d))
}
subgroups <- c(groups[-1],list(GA35=ga35))
omni <- do.call(rbind,lapply(names(subgroups),function(field) {
  m <- !is.na(iso[[field]])
  d <- domain(design,m)
  a <- data.frame(Block=field,Type="All 16 pathogen groups",overall_test(d,field,"Organism",paste("Rao–Scott",field)))
  if(field=="GA35") return(a)
  b <- domain(d,d$variables$Gram %in% c("Gram-positive","Gram-negative"))
  b$variables$Gram <- factor(b$variables$Gram)
  rbind(a,data.frame(Block=field,Type="Gram-positive vs Gram-negative",overall_test(b,field,"Gram",paste("Gram",field))))
}))
coef_rows <- list()
s8 <- do.call(rbind,lapply(names(subgroups),function(field) {
  d0 <- domain(design,!is.na(iso[[field]]))
  ans <- do.call(rbind,lapply(org_names[1:12],function(o) {
    label <- paste("S8",field,o)
    d <- make_y(d0,as.integer(d0$variables$Organism==o))
    y <- d$variables$Endpoint
    f <- reformulate(field,"Endpoint")
    z <- attempt(label,survey::svyglm(f,d,family=quasipoisson(link="log")))
    k <- if(is.null(z)) NULL else attempt(label,survey::regTermTest(z,reformulate(field),method="Wald"))
    v <- if(is.null(k)) rep(NA_real_,4) else c(as.numeric(k$Ftest),k$df,k$ddf,k$p)
    levels_with_zero <- any(tapply(y,d$variables[[field]],sum)==0)
    if(!is.null(z) && isTRUE(z$converged)) {
      b <- coef(z)[-1]; se <- sqrt(diag(vcov(z)))[-1]
      adequate <- is.finite(b)&is.finite(se)&abs(b)<5&se>0&se<5
      coef_rows[[length(coef_rows)+1L]] <<- data.frame(Block=field,Pathogen=o,Contrast=names(b),
        Reference=levels(d$variables[[field]])[1],PR=ifelse(adequate,exp(b),NA),
        Low=ifelse(adequate,exp(b-qt(.975,z$df.residual)*se),NA),
        High=ifelse(adequate,exp(b+qt(.975,z$df.residual)*se),NA),DF=z$df.residual)
    }
    data.frame(Block=field,Pathogen=o,N=length(y),Events=sum(y),F=if(levels_with_zero) NA_real_ else v[1],Numerator_DF=v[2],
      Denominator_DF=v[3],P=if(levels_with_zero) NA_real_ else v[4],Design_DF=survey::degf(d),
      Zero_event_level=levels_with_zero,Status=if(levels_with_zero) "Not estimable: zero-event GA level; boundary Wald excluded" else "Wald model")
  }))
  ans$Q <- p.adjust(ans$P,"BH")
  ans$Family_tests <- sum(!is.na(ans$P))
  ans
}))
s8_effect <- do.call(rbind,coef_rows)
table8_counts <- by_groups(function(field,level,m) do.call(rbind,lapply(org_names[1:12],function(o) {
  data.frame(Block=field,Level=level,Pathogen=o,Events=sum(m&iso$Organism==o),N=sum(m),
    Weighted=weighted.mean(as.integer(iso$Organism[m]==o),iso$wg[m]))
})),subgroups)
# Refit the accepted period contrasts with both designs. The BH family is per domain and mode.
s12 <- by_groups(function(field,level,m) {
  rows <- do.call(rbind,lapply(org_names[1:12],function(o) {
    do.call(rbind,lapply(c("Weighted","Unweighted"),function(mode) {
      d <- domain(if(mode=="Weighted") design else unit,m)
      d <- make_y(d,as.integer(d$variables$Organism==o))
      fit <- fit_model(d,"Period",paste("Pathogen period",level,o,mode),3L)
      data.frame(Block=field,Level=level,Pathogen=o,Analysis=mode,fit)
    }))
  }))
  rows$Q <- ave(rows$P,rows$Analysis,FUN=function(p) p.adjust(p,"BH"))
  rows
})

s12_omni <- by_groups(function(field,level,m) {
  d <- domain(design,m)
  do.call(rbind,lapply(c(TRUE,FALSE),function(w) data.frame(Block=field,Level=level,
    Analysis=if(w) "Weighted" else "Unit weights",overall_test(d,"Period","Organism",paste("S12 omnibus",level),w))))
})

# Boundary-case diagnostic results are separate from the formal Wald estimates and BH families.
targets <- s8$Pathogen[s8$Block=="GA" & s8$Zero_event_level]
sparse_ga <- do.call(rbind,lapply(targets,function(p) {
  d <- domain(make_y(design,as.integer(iso$Organism==p)),!is.na(iso$GA))
  rs <- attempt(paste("Diagnostic Rao-Scott",p),survey::svychisq(~Endpoint+GA,d,statistic="F"))
  observed <- table(d$variables$Endpoint,d$variables$GA)
  data.frame(Pathogen=p,GA=colnames(observed),Events=as.integer(observed["1",]),
    Eligible=as.integer(colSums(observed)),Rao_Scott_F=as.numeric(rs$statistic),
    Numerator_DF=as.numeric(rs$parameter[1]),Denominator_DF=as.numeric(rs$parameter[2]),
    Diagnostic_P=as.numeric(rs$p.value))
}))

# Pathogen-identification row of the mixed-domain missingness table.
# Unit-weight second-order Rao-Scott tests compare all years, macro-regions or outcomes.
pattern <- paste(trim(raw[["死亡"]][prim]),trim(raw[["死亡或放弃治疗"]][prim]),sep="/")
outcomes <- c("Survived to discharge","In-hospital death","Withdrawal of care")
disposition <- setNames(c(outcomes,NA_character_),c("否/否","是/是","否/是","不详/不详"))
stopifnot(all(pattern %in% names(disposition)))
full_missing <- make_y(unit,as.integer(iso$Organism=="Unidentified organism(s)"))
full_missing$variables$Calendar_year <- factor(iso$Year,levels=2015:2024)
full_missing$variables$Outcome <- factor(disposition[pattern][match(iso$InfantID,iso$InfantID[prim])],levels=outcomes)
pathogen_missing <- do.call(rbind,lapply(c("Calendar_year","Region","Outcome"),function(group) {
  d <- domain(full_missing,!is.na(full_missing$variables[[group]]))
  d$variables[[group]] <- droplevels(d$variables[[group]])
  events <- d$variables$Endpoint
  denominators <- table(d$variables[[group]])
  n_event <- tapply(events,d$variables[[group]],sum)
  percentage <- 100*n_event/as.integer(denominators)
  t <- overall_test(d,group,"Endpoint",paste("Pathogen identification missingness",group))
  data.frame(Group=group,N=sum(denominators),Unidentified=sum(n_event),
             Min_pct=min(percentage),Max_pct=max(percentage),t)
}))
