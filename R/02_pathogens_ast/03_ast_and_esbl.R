# Publication scopes contain predefined endpoint labels and eligibility flags only.
read_protocol <- function(name) {
  switch(name,
    "ast_fields.csv"          = static_cfg$get_ast_fields()[, c("Drug", "Field")],
    "ast_reporting_scope.csv" = static_cfg$get_ast_reporting_scope(),
    "intrinsic_resistance.csv"= static_cfg$get_ast_intrinsic_rules(),
    "amr_trend_scope.csv"     = static_cfg$get_amr_trend_scope(),
    stop("Unknown static protocol table: ", name, call. = FALSE)
  )
}
pair_key <- function(pathogen,drug) paste(pathogen,drug,sep="::")
classes <- list(Glycopeptides=c("Vancomycin","Teicoplanin"),
  "Third-generation cephalosporins"=c("Ceftazidime","Cefotaxime"),Carbapenems=c("Imipenem","Meropenem"),
  Aminoglycosides=c("Gentamicin","Amikacin"),Fluoroquinolones=c("Ciprofloxacin","Levofloxacin"),
  Macrolides=c("Azithromycin","Erythromycin"))
layout <- read_protocol("ast_fields.csv")
fields <- setNames(layout$Field,layout$Drug)
mapping <- read_protocol("ast_reporting_scope.csv")
mapping$In_Table3 <- as.logical(mapping$In_Table3)
keys <- pair_key(mapping$Pathogen,mapping$Drug)
matrix_r <- read_protocol("intrinsic_resistance.csv")
matrix_keys <- pair_key(matrix_r$Pathogen,matrix_r$Drug)
stopifnot(!anyDuplicated(keys), !anyDuplicated(matrix_keys), !anyNA(mapping$In_Table3),
          setequal(keys,as.vector(outer(org_names[1:12],names(fields),pair_key))),
          all(matrix_keys %in% keys), setequal(fields,names(raw)[46:74]))
sir_map <- c("敏感"=1L,"中介"=2L,"耐药"=3L,S=1L,I=2L,R=3L)
observed <- lapply(fields,function(field) unname(sir_map[canon(raw[[field]])]))
rebuilt <- observed

# Construct classes after excluding matrix-R components, keeping S/I/R ordered R > I > S.
for (drug in names(classes)) {
  parts <- lapply(classes[[drug]],function(agent) {
    replace(observed[[agent]],pair_key(iso$Organism,agent) %in% matrix_keys,NA_integer_)
  })
  rebuilt[[drug]] <- do.call(pmax,c(parts,list(na.rm=TRUE)))
}
mapping$Matrix_R <- keys %in% matrix_keys
mapping$Allowed_components <- NA_integer_
for (i in which(mapping$Drug %in% names(classes))) {
  agents <- classes[[mapping$Drug[i]]]
  mapping$Allowed_components[i] <- sum(!pair_key(mapping$Pathogen[i],agents) %in% matrix_keys)
}
allowed <- mapping$Allowed_components
single <- !is.na(allowed) & allowed==1L
zero <- !is.na(allowed) & allowed==0L & !mapping$Matrix_R
excluded <- mapping$Matrix_R | single | zero
mapping$Dash_basis <- ifelse(mapping$Matrix_R,"Matrix_R",ifelse(single,"Class_Single_Eligible_Component",
                               ifelse(zero,"Class_Zero_Eligible_Component","")))
values <- lapply(seq_along(keys),function(i) {
  x <- rebuilt[[mapping$Drug[i]]]
  x[excluded[i] | iso$Organism!=mapping$Pathogen[i]] <- NA_integer_
  x
})
ast <- do.call(rbind,lapply(seq_along(keys),function(i) {
  x <- values[[i]]
  eligible <- !is.na(x)
  n <- sum(eligible)
  y <- as.integer(x==3L)
  fit <- estimate(domain(make_y(design,y),eligible),keys[i])
  if (!n) {
    fit$DF <- 0
    fit$CI_status <- if (excluded[i]) paste("Suppressed",mapping$Dash_basis[i]) else "No valid susceptibility results"
  }
  wt <- sum(iso$wg[eligible])
  data.frame(mapping[i,],Organism_total=sum(iso$Organism==mapping$Pathogen[i]),
    n_S=sum(x %in% 1L),n_I=sum(x %in% 2L),n_R=sum(x %in% 3L),N_tested=n,
    Weighted_S=if(n) sum(iso$wg[x %in% 1L])/wt else NA_real_,
    Weighted_I=if(n) sum(iso$wg[x %in% 2L])/wt else NA_real_,fit,
    Display_status=if(excluded[i]) "—" else if(!n) "ND" else if(is.na(fit$Low)) "Numeric (NE)" else "Numeric (CI)")
}))
stopifnot(all(ast$n_S+ast$n_I+ast$n_R==ast$N_tested),all(ast$n_R==ast$Events),
          all(abs(ast$Weighted_S+ast$Weighted_I+ast$Weighted-1)<1e-10,na.rm=TRUE))

# Count testing opportunities only for reportable combination types.
ast_missing <- function(field,levels) {
  do.call(rbind,lapply(levels,function(lv) {
    selected <- is_level(field,lv)
    label <- if (field=="Period") c(Early="2015–2019",Late="2020–2024")[[lv]] else as.character(lv)
    do.call(rbind,lapply(which(!excluded),function(i) {
      x <- values[[i]][selected & iso$Organism==mapping$Pathogen[i]]
      n <- length(x)
      valid <- sum(!is.na(x))
      data.frame(Pathogen=mapping$Pathogen[i],Drug=mapping$Drug[i],
        Type=if(mapping$Drug[i] %in% names(classes)) "Class" else "Single agent",
        Domain=field,Level=label,Organism_total=n,S=sum(x %in% 1L),I=sum(x %in% 2L),
        R=sum(x %in% 3L),Valid=valid,Missing=n-valid,Valid_rate=if(n) valid/n else NA_real_,
        Missing_rate=if(n) (n-valid)/n else NA_real_)
    }))
  }))
}
missing_tables <- list(AST_Missing_Overall=ast_missing("Overall","Overall"),
  AST_Missing_by_Year=ast_missing("Year",2015:2024),
  AST_Missing_by_Period=ast_missing("Period",c("Early","Late")),
  AST_Missing_by_Region=ast_missing("Region",region))
missing_detail <- do.call(rbind,missing_tables)
missing_cols <- c("Organism_total","S","I","R","Valid","Missing")
miss_overall <- missing_tables[[1L]]
for (z in missing_tables[-1L]) {
  totals <- rowsum(z[missing_cols],pair_key(z$Pathogen,z$Drug))
  stopifnot(all(totals[pair_key(miss_overall$Pathogen,miss_overall$Drug),]==as.matrix(miss_overall[missing_cols])))
}
stopifnot(all(miss_overall$Valid==ast$N_tested[!excluded]))
missing_summary <- aggregate(missing_detail[missing_cols],missing_detail[c("Domain","Level","Type")],sum)
missing_summary$Missing_rate <- with(missing_summary,Missing/Organism_total)

# Temporal AMR models retain the endpoint-specific inclusion scope of the accepted trend analysis.
trend_scope <- read_protocol("amr_trend_scope.csv")
trend_scope$Source_row <- as.integer(trend_scope$Source_row)
trend_scope$Frozen_scope <- as.logical(trend_scope$Frozen_scope)
stopifnot(!anyNA(trend_scope$Frozen_scope), !anyDuplicated(trend_scope[c("Source_row","Pathogen")]))
modes <- list(Unweighted=unit,Weighted=design)
s13 <- do.call(rbind,lapply(unique(trend_scope$Source_row),function(k) {
  spec <- trend_scope[trend_scope$Source_row==k,]
  panel <- unique(spec$Trend_group)
  drug <- unique(spec$Drug)
  stopifnot(length(panel)==1L, length(drug)==1L, drug %in% names(observed))
  group <- if(panel %in% c("Gram-positive","Gram-negative")) iso$Gram==panel else iso$Organism==panel
  stopifnot(setequal(as.character(unique(iso$Organism[group])),spec$Pathogen))
  included <- group & iso$Organism %in% spec$Pathogen[spec$Frozen_scope] & !is.na(observed[[drug]])
  endpoint <- as.integer(observed[[drug]]==3L)
  do.call(rbind,lapply(names(modes),function(mode) {
    d <- domain(make_y(modes[[mode]],endpoint),included)
    data.frame(Source_row=k,Pathogen=panel,Drug=drug,Analysis=mode,
      fit_model(d,"Year",paste("AMR trend",panel,drug,mode),5L))
  }))
}))
bh <- function(p) {
  keep <- !is.na(p)
  p[keep] <- p.adjust(p[keep],"BH")
  p
}
family <- interaction(s13$Pathogen,s13$Analysis,drop=TRUE)
s13$Q <- ave(s13$P,family,FUN=bh)
s13$FDR_tests <- ave(!is.na(s13$P),family,FUN=sum)

# ESBL is a site-reported phenotype within the two original pathogen groups.
# Preserve the signed-off abstraction convention: reported positive flag versus no positive flag.
esbl_positive <- as.integer(canon(raw[["ESBL"]]) %in% "是")
esbl_taxa <- c("E. coli","Klebsiella spp.")
esbl_design <- make_y(design,esbl_positive)
s5 <- by_groups(function(field,level,m) {
  do.call(rbind,lapply(esbl_taxa,function(p) {
    fit <- estimate(domain(esbl_design,m & iso$Organism==p),paste("ESBL positivity",p,level))
    data.frame(Pathogen=p,Block=field,Level=level,fit,Not_reported_positive=fit$N-fit$Events)
  }))
})
esbl_trends <- do.call(rbind,lapply(esbl_taxa,function(p) {
  eligible <- iso$Organism==p
  do.call(rbind,lapply(names(modes),function(mode) {
    fit <- fit_model(domain(make_y(modes[[mode]],esbl_positive),eligible),"Year",
                     paste("ESBL trend",p,mode),5L)
    data.frame(Pathogen=p,Drug="ESBL phenotype",Analysis=mode,Q=NA_real_,Organism_total=sum(eligible),fit)
  }))
}))
esbl_bh <- esbl_trends[c("Pathogen","Analysis","P")]
esbl_bh$Exploratory_BH_q <- ave(esbl_bh$P,esbl_bh$Analysis,FUN=bh)
kleb <- iso$Organism=="Klebsiella spp."
yearly <- do.call(rbind,lapply(2015:2024,function(y) {
  m <- kleb & iso$Year==y
  data.frame(Year=y,Eligible=sum(m),Positive=sum(esbl_positive[m]),PSU=length(unique(iso$Site[m])),
    Positive_PSU=length(unique(iso$Site[m & esbl_positive==1L])))
}))
leave_one <- do.call(rbind,lapply(names(modes),function(mode) {
  d <- make_y(modes[[mode]],esbl_positive)
  do.call(rbind,lapply(sort(unique(iso$Site)),function(h) {
    hit <- iso$Site==h
    fit <- fit_model(domain(d,kleb & !hit),"Year",paste("ESBL leave site",mode),5L)
    beta <- esbl_trends$Beta[esbl_trends$Pathogen=="Klebsiella spp." & esbl_trends$Analysis==mode]
    data.frame(Analysis=mode,Excluded_PSU=h,Stratum=as.character(iso$Region[which(hit)[1L]]),
      Removed_Klebsiella=sum(hit & kleb),Removed_Positive=sum(esbl_positive[hit & kleb]),
      Delta_Beta=fit$Beta-beta,fit)
  }))
}))
