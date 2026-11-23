# Full-precision aggregate outputs for verification and the publication tables.
write_table <- function(name, data) {
  readr::write_csv(data,file.path(output_dir,paste0(name,".csv")),na="")
}
numeric_tables <- c(list(Cohort=counts,Table2_Pathogen=table2,Table2_Overall_Tests=omni,
  Table3_S6_AST=ast,S5_ESBL=s5,S8_Counts=table8_counts,S8_Effects=s8_effect,S8_Tests=s8,
  S8_Sparse_Diagnostics=sparse_ga,Pathogen_Period=s12,Pathogen_Period_Overall_Tests=s12_omni,
  AMR_Trend=s13,ESBL_Trend=esbl_trends,ESBL_Exploratory_BH=esbl_bh,Klebsiella_ESBL_by_Year=yearly,
  Klebsiella_ESBL_Leave_One=leave_one,AST_Missing_Opportunity_Summary=missing_summary,
  S20_Pathogen_Identification=pathogen_missing),missing_tables)
for (name in names(numeric_tables)) write_table(name,numeric_tables[[name]])
# Publication cells retain the precision and status conventions used in the tables.
dash <- "—"
pct <- function(x) sprintf("%.1f",100*x)
p_text <- function(x) ifelse(is.na(x),dash,ifelse(x<.001,"<0.001",sprintf("%.3f",x)))
ci_text <- function(value,low,high,missing_ci=dash) {
  lims <- ifelse(is.na(low),missing_ci,paste0(pct(low),"–",pct(high)))
  ifelse(is.na(value),dash,paste0(pct(value)," (",lims,")"))
}
pr_text <- function(value,lo,hi) {
  mapply(function(v,l,h) {
    if (is.na(v) || is.na(l) || is.na(h)) return(dash)
    digits <- if (any(round(c(l,h),2)==1 & c(l,h)!=1)) 3L else 2L
    paste0(sprintf("%.2f",v)," (",formatC(l,digits=digits,format="f"),"–",
           formatC(h,digits=digits,format="f"),")")
  },value,lo,hi,USE.NAMES=FALSE)
}
display_t2 <- table2[c("Block","Level","Pathogen","Events")]
display_t2$Weighted_CI <- with(table2,ci_text(Weighted,Low,High))
write_table("Table2_Display",display_t2)
p2 <- rbind(data.frame(Block=omni$Block,Comparison=omni$Type,P=p_text(omni$P)),
            data.frame(Block=s8$Block,Comparison=s8$Pathogen,P=p_text(s8$P)))
write_table("Table2_Display_P",p2[p2$Block!="GA35",])
t3 <- ast[ast$In_Table3,]
t3_text <- with(t3,ifelse(Display_status %in% c(dash,"ND"),Display_status,
           paste0(n_R,"/",N_tested,"\n",ci_text(Weighted,Low,High,"NE"))))
write_table("Table3_Display",data.frame(t3[c("Drug","Pathogen")],Resistance=t3_text))

display_s5 <- s5[order(match(s5$Pathogen,esbl_taxa)),]
write_table("S5_Display",with(display_s5,data.frame(Pathogen,Block,Level,
  ESBL_positive_eligible=paste0(Events,"/",N),Weighted_positivity_CI=ci_text(Weighted,Low,High))))
display_s6 <- do.call(rbind,lapply(c("S","R","I"),function(category) {
  n <- ast[[paste0("n_",category)]]
  p <- ast[[if(category=="R") "Weighted" else paste0("Weighted_",category)]]
  cell <- ifelse(ast$Display_status %in% c(dash,"ND"),ast$Display_status,
                 paste0(n,"/",ast$N_tested," (",pct(p),"%)"))
  data.frame(ast[c("Drug","Pathogen")],AST=category,Results=cell)
}))
write_table("S6_Display",display_s6)

s8a <- data.frame(Pathogen=org_names[1:12])
for (field in c("Region","GA","Delivery")) {
  z <- s8[s8$Block==field,]
  s8a[[paste0(field,"_q")]] <- p_text(z$Q[match(s8a$Pathogen,z$Pathogen)])
}
ga_counts <- table8_counts[table8_counts$Block=="GA35",]
ga_cell <- with(ga_counts,paste0(Events," (",pct(Weighted),")"))
ga_ref <- s8_effect[s8_effect$Block=="GA35",]
ga_tests <- s8[s8$Block=="GA35",]
s8b <- data.frame(Pathogen=ga_tests$Pathogen,
  Earlier_GA=ga_cell[ga_counts$Level==ga35[1L]],Later_GA=ga_cell[ga_counts$Level==ga35[2L]],
  P=p_text(ga_tests$P),q=p_text(ga_tests$Q))
s8b$PR_CI <- with(ga_ref,pr_text(PR,Low,High))[match(s8b$Pathogen,ga_ref$Pathogen)]
write_table("S8_Panel_A",s8a)
write_table("S8_Panel_B",s8b)

trend_cells <- function(data,keys,scale_years) {
  field <- if (scale_years==5L) c("PR5","Low5","High5") else c("PR","Low","High")
  u <- data[data$Analysis=="Unweighted",]
  w <- data[data$Analysis=="Weighted",]
  id <- function(x) do.call(paste,c(x[keys],sep=" :: "))
  w <- w[match(id(u),id(w)),]
  stopifnot(identical(id(u),id(w)))
  data.frame(u[keys],Unweighted_PR_CI=do.call(pr_text,unname(as.list(u[field]))),
    Unweighted_P=p_text(u$P),Unweighted_q=p_text(u$Q),
    Weighted_PR_CI=do.call(pr_text,unname(as.list(w[field]))),Weighted_P=p_text(w$P),Weighted_q=p_text(w$Q))
}
write_table("Pathogen_Period_Display",trend_cells(s12,c("Block","Level","Pathogen"),1L))
trend_cols <- c("Pathogen","Drug","Analysis","PR5","Low5","High5","P","Q")
write_table("AMR_ESBL_Trend_Display",trend_cells(rbind(s13[trend_cols],esbl_trends[trend_cols]),
                                                  c("Pathogen","Drug"),5L))

missing_range <- function(x) paste0(sprintf("%.1f%%",x$Min_pct),"–",sprintf("%.1f%%",x$Max_pct))
y <- pathogen_missing[pathogen_missing$Group=="Calendar_year",]
r <- pathogen_missing[pathogen_missing$Group=="Region",]
o <- pathogen_missing[pathogen_missing$Group=="Outcome",]
write_table("S20_Pathogen_Identification_Display",data.frame(Variable="Unidentified pathogen",
  Missing=sprintf("%d/%d (%.1f%%)",y$Unidentified,y$N,100*y$Unidentified/y$N),
  Annual_range=missing_range(y),Annual_P=p_text(y$P),Regional_range=missing_range(r),Regional_P=p_text(r$P),
  Outcome_range=missing_range(o),Outcome_P=p_text(o$P)))
writeLines(capture.output(sessionInfo()),file.path(output_dir,"sessionInfo.txt"))

# Final-manuscript numbered aliases (Table 2, Table 3, S4, S5, S7, S12, S13 Tables)
static_cfg$copy_required(file.path(output_dir, "S5_Display.csv"),  file.path(output_dir, "S4_Table_ESBL_Display.csv"), overwrite = TRUE)
static_cfg$copy_required(file.path(output_dir, "S6_Display.csv"),  file.path(output_dir, "S5_Table_AST_Full_Display.csv"), overwrite = TRUE)
static_cfg$copy_required(file.path(output_dir, "S8_Panel_A.csv"),  file.path(output_dir, "S7_Table_PanelA_Display.csv"), overwrite = TRUE)
static_cfg$copy_required(file.path(output_dir, "S8_Panel_B.csv"),  file.path(output_dir, "S7_Table_PanelB_Display.csv"), overwrite = TRUE)
static_cfg$copy_required(file.path(output_dir, "Pathogen_Period_Display.csv"), file.path(output_dir, "S12_Table_Pathogen_Period_Display.csv"), overwrite = TRUE)
static_cfg$copy_required(file.path(output_dir, "AMR_ESBL_Trend_Display.csv"),  file.path(output_dir, "S13_Table_AMR_ESBL_Trends_Display.csv"), overwrite = TRUE)
