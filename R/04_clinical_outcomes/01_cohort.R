# Input rows are isolate-level. A blank ID denotes a continuation row.
# Preserve the verified primary-infant row and retain one observation per infant.
read_text_csv <- function(path) {
  result <- readr::read_csv(path, col_types=readr::cols(.default=readr::col_character()),
    locale=readr::locale(encoding="UTF-8"), name_repair="minimal", na=character(),
    trim_ws=FALSE, show_col_types=FALSE, progress=FALSE)
  if (nrow(readr::problems(result))) stop("CSV parsing failed: ", basename(path))
  result
}
clean <- function(x) trimws(as.character(x), whitespace="[\\h\\v]")
canon <- function(x) {
  x <- clean(x)
  x[tolower(x) %in% c("", "na", "n/a", "unk", "unknown", "不详", "未知", "缺失")] <- NA_character_
  x
}
categorise <- function(x, id) {
  x <- canon(x)
  categories <- if (id %in% c("IVH","NEC","Meningitis")) c("疑诊","可疑") else character()
  if (any(!is.na(x) & !x %in% c("否","是",categories))) stop("Unexpected outcome category: ",id)
  ifelse(is.na(x), NA_integer_, ifelse(x=="否", 0L, ifelse(x=="是",1L,2L)))
}

region_names <- c("Municipality","Eastern","Central","Western")
regions <- c("Overall",region_names)
region_map <- setNames(region_names,c("一层","二层","三层","四层"))
year_set <- 2015:2024
ga_names <- c("Extremely preterm infants","Very preterm infants","Moderate to late preterm infants","Term infants")
ga_labels_source <- c("超早产儿","极早产儿","中晚期早产儿","足月儿")
ga2 <- c("≥35 weeks","<35 weeks")
source_fields <- c(D72="生后72小时内是否死亡",C72="生后72小时内是否死亡或放弃",
  D7="生后7天内是否死亡",C7="生后7天内是否死亡或放弃",Death="死亡",Composite="死亡或放弃治疗",
  BPD="BPD（支气管肺发育不良）",Shock="休克",LOS="院内晚发败血症",
  ROP="ROP（早产儿视网膜病变）",IVH="IVH（早产儿脑室内出血）",NEC="NEC（新生儿坏死性小肠结肠炎）")
labels <- c(D72="Death within 72 h of birth",D7="Death within 7 days of birth",
  Death="All-cause in-hospital death",C72="Death or withdrawal of treatment within 72 h of birth",
  C7="Death or withdrawal of treatment within 7 days of birth",Withdrawal="Withdrawal of treatment overall",
  Composite="Death or withdrawal of treatment overall",BPD="BPD",Shock="Shock",
  LOS="In-hospital late-onset sepsis",ROP="ROP",IVH="Confirmed IVH",NEC="Confirmed NEC")
outcome_codes <- names(labels)
disposition_codes <- outcome_codes[1:7]
respiratory_fields <- c(
  "Conventional mechanical ventilation"="是否使用常频",
  "High-frequency ventilation"="是否使用高频",
  CPAP="是否使用CPAP（持续正压治疗）", NIPPV="是否使用NIPPV（无创正压通气）",
  HFNC="是否使用HFNC（经鼻高流量湿化氧疗）")
management_fields <- c(
  "Glucocorticoid therapy"="是否使用GC", "Surfactant administration"="是否使用PS",
  "Seizures within 72 h of birth"="生后72h内出现惊厥", Meningitis="化脓性脑膜炎",
  RDS="RDS(新生儿呼吸窘迫综合症）", "Vasoactive therapy within 72 h of birth"="生后72h内使用血管活性药物")
regimen_fields <- c("First use after Birth","Use within 72 h","Use within 72 h–7 days")
modification_fields <- c("72h内是否更改","72h-7d是否更改")
prophylaxis_field <- "是否预防性使用抗真菌药物"
classic_field <- "是否经典方案"

raw <- read_text_csv(case_path)
clinical_fields <- c("胎龄w",unname(source_fields),unname(respiratory_fields),
  unname(management_fields),regimen_fields,modification_fields,prophylaxis_field,classic_field)
key_fields <- c("医院名称","入院时间","分层")
required <- c("ID",key_fields,clinical_fields,"胎龄(w)","胎龄分类","早产儿分类")
if (!all(required %in% names(raw))) stop("Missing case column(s): ",paste(setdiff(required,names(raw)),collapse=", "))
id <- canon(raw[["ID"]])
primary <- which(!is.na(id))
secondary <- which(is.na(id))
stopifnot(nrow(raw)==1196L, length(primary)==1189L, !anyDuplicated(id[primary]), length(secondary)==7L)
parent <- primary[findInterval(secondary,primary)]
stopifnot(length(parent)==length(secondary), all(parent==secondary-1L))
for (field in key_fields) {
  x <- canon(raw[[field]])
  stopifnot(!anyNA(x[primary]), all(is.na(x[secondary]) | x[secondary]==x[parent]))
}
for (field in setdiff(clinical_fields, classic_field)) {
  x <- canon(raw[[field]])
  # Continuation fields may be blank. Nonblank duplicates must agree.
  has_secondary <- !is.na(x[secondary])
  if (!isTRUE(all(!has_secondary | x[secondary]==x[parent]))) {
    stop("Check the continuation row for ",field)
  }
}

site <- read_text_csv(site_path)
stopifnot(all(c("S4 center","Hospital","S4 stratum") %in% names(site)))
site <- site[grepl("^Site [0-9]+$",clean(site[["S4 center"]])),]
hospital_id <- clean(site[["S4 center"]])
site_name <- clean(site[["Hospital"]])
stopifnot(nrow(site)==30L, !anyDuplicated(hospital_id), !anyDuplicated(site_name))
# Match named or coded hospitals. Codes suffice when depositing deidentified inputs.
hospital <- clean(raw[["医院名称"]][primary])
map_row <- match(hospital,site_name)
map_row[is.na(map_row)] <- match(hospital[is.na(map_row)],hospital_id)
region <- unname(region_map[clean(raw[["分层"]][primary])])
year <- as.integer(clean(raw[["入院时间"]][primary]))
stopifnot(!anyNA(map_row),!anyNA(region),all(region==clean(site[["S4 stratum"]][map_row])),
  all(year %in% year_set), length(unique(hospital_id[map_row]))==30L)
inf <- data.frame(HospitalID=hospital_id[map_row],Region=factor(region,levels=region_names),
  Year=year,t=year-2015L,w=1)
stopifnot(identical(as.integer(table(inf$Region)),c(299L,279L,301L,310L)))

values <- function(field) raw[[field]][primary]
inf$GAWeeks <- as.numeric(canon(values("胎龄w")))
inf$GAGroup <- cut(inf$GAWeeks,c(-Inf,28,32,37,Inf),right=FALSE,labels=ga_names)
inf$GA35 <- factor(ifelse(inf$GAWeeks<35,ga2[2L],ga2[1L]),levels=ga2)
stopifnot(sum(is.na(inf$GAWeeks))==3L,
  identical(as.integer(inf$GAGroup),match(canon(values("早产儿分类")),ga_labels_source)),
  identical(as.integer(inf$GAGroup),as.integer(canon(values("胎龄分类")))),
  identical(ifelse(is.na(inf$GA35),NA_integer_,3L-as.integer(inf$GA35)),as.integer(canon(values("胎龄(w)")))))
for (code in names(source_fields)) {
  state <- categorise(values(source_fields[[code]]),code)
  inf[[paste0(code,"_state")]] <- state
  inf[[code]] <- as.integer(state==1L)
}
for (window in list(c("D72","C72"),c("D7","C7"),c("Death","Composite"))) {
  death <- inf[[window[1L]]]
  composite <- inf[[window[2L]]]
  stopifnot(identical(is.na(death),is.na(composite)),all(death<=composite,na.rm=TRUE))
}
inf$Withdrawal <- as.integer(inf$Death==0L & inf$Composite==1L)
stopifnot(sum(inf$Death,na.rm=TRUE)==98L, sum(inf$Withdrawal,na.rm=TRUE)==128L,
  sum(inf$Composite,na.rm=TRUE)==226L)
group_mask <- function(data, group) if(group=="Overall") rep(TRUE,nrow(data)) else data$Region==group
n_fmt <- function(x) format(x,big.mark=",",trim=TRUE,scientific=FALSE)
risk_text <- function(events, n) {
  if (n==0L) return("—")
  sprintf("%s/%s (%.1f%%)",n_fmt(events),n_fmt(n),100*events/n)
}

