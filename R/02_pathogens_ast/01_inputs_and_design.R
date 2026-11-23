trim <- function(x) trimws(as.character(x), whitespace = "[\\h\\v]")
blank <- function(x) is.na(x) | tolower(trim(x)) %in% c("", "na", "n/a", "unknown", "unk", "不详", "未知", "缺失")
canon <- function(x) {
  x <- trim(x)
  x[blank(x)] <- NA_character_
  x
}
read_csv <- function(path, encoding = "UTF-8") {
  x <- readr::read_csv(path, locale = readr::locale(encoding = encoding), show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character()), na = character(), trim_ws = FALSE, name_repair = "minimal")
  if (nrow(readr::problems(x))) stop("Malformed input: ", path)
  x
}
pick <- function(x, cols) {
  i <- match(cols, names(x))
  if (anyNA(i)) stop("Input fields missing: ", paste(cols[is.na(i)], collapse = ", "))
  as.data.frame(setNames(lapply(i, function(j) trim(x[[j]])), cols), optional = TRUE)
}
region <- c("Municipality", "Eastern", "Central", "Western")
region_map <- setNames(region, c("一层", "二层", "三层", "四层"))

# Study input and survey design. Original case rows are kept in their supplied order.
raw <- read_csv(inputs[["Cases"]])
taxon <- "菌群整理归类（归类到属）"
key <- c("ID", "医院名称", "入院时间", "分层")
linked <- pick(raw, c(key, taxon, "Gram"))
prim <- which(!blank(linked$ID))
extra <- which(blank(linked$ID))
stopifnot(ncol(raw)==95L, length(prim)==1189L, length(extra)==7L,
          !anyDuplicated(linked$ID[prim]), identical(trim(raw[[44]]), trim(raw[[76]])),
          all(vapply(linked[extra, key], function(v) all(blank(v)), TRUE)),
          all(extra > 1L), all((extra-1L) %in% prim))
linked[extra,key] <- linked[extra-1L,key]
stopifnot(identical(as.integer(table(table(linked$ID))),c(1182L,7L)))


site <- read_csv(inputs[["Sites"]])
site <- site[grepl("^Site [0-9]+$", site[["S4 center"]]), ]
ind <- match(linked[["医院名称"]], site[["Hospital"]])
stopifnot(nrow(site) == 30L && !anyNA(ind) && !anyDuplicated(site[["Hospital"]]) &&
  !anyDuplicated(site[["S4 center"]]) && all(site[["S4 stratum"]][ind] == region_map[linked[["分层"]]]) &&
  all(as.integer(linked[["入院时间"]]) %in% 2015:2024))

orgs <- c("大肠埃希菌", "GBS", "CoNS", "克雷伯菌属", "肠球菌属", "产单核细胞李斯特菌", "链球菌属",
  "金黄色葡萄球菌", "肠杆菌属", "不动杆菌属", "假单胞菌属", "沙雷菌属", "其他革兰氏阳性菌",
  "其他革兰氏阴性菌", "真菌", "菌不详")
org_names <- c("E. coli", "GBS", "CoNS", "Klebsiella spp.", "Enterococcus spp.", "L. monocytogenes",
  "Other Streptococcus spp.", "S. aureus", "Enterobacter spp.", "Acinetobacter spp.", "Pseudomonas spp.",
  "Serratia spp.", "Other Gram-positive bacteria", "Other Gram-negative bacteria", "Fungi", "Unidentified organism(s)")
ga <- c("Extremely preterm infants", "Very preterm infants", "Moderate to late preterm infants", "Term infants")
del <- c("Caesarean section", "Vaginal delivery")
ga35 <- c("<35 weeks", "≥35 weeks")
mapped <- function(v, source, dest) unname(setNames(dest, source)[trim(v)])
iso <- data.frame(InfantID = linked$ID, Site = site[["S4 center"]][ind],
  Region = factor(site[["S4 stratum"]][ind], levels = region),
  Year = as.integer(linked[["入院时间"]]), Organism = factor(mapped(linked[[taxon]], orgs, org_names), levels = org_names),
  Gram = mapped(linked$Gram, c("G+", "G-", "真菌", "不详"), c("Gram-positive","Gram-negative","Fungi","Unidentified")))
iso$Period <- factor(ifelse(iso$Year >= 2020, "Late", "Early"), levels = c("Early","Late"))

# A continuation row inherits infant-level fields from its primary row.
sg <- c("早产儿分类","分娩方式","胎龄(w)")
current_cov <- lapply(pick(raw,sg), function(v) {
  v <- canon(v)
  v[extra] <- v[extra-1L]
  v
})


iso$GA <- factor(mapped(current_cov[[1]], c("超早产儿","极早产儿","中晚期早产儿","足月儿"), ga), levels = ga)
iso$Delivery <- factor(mapped(current_cov[[2]], c("剖宫产","自然分娩"), del), levels = del)
iso$GA35 <- factor(mapped(current_cov[[3]], c("1","2"), ga35), levels = ga35)
counts <- data.frame(Region = region, Infants = as.integer(table(iso$Region[prim])),
  Isolates = as.integer(table(iso$Region)), Hospitals = as.integer(table(site[["S4 stratum"]])[region]))
stopifnot(!anyNA(iso$Organism), !anyNA(iso$Gram), identical(counts$Isolates,c(299L,281L,304L,312L)),
  identical(counts$Infants,c(299L,279L,301L,310L)), identical(as.integer(table(iso$GA)),c(69L,200L,306L,618L)),
  identical(as.integer(table(iso$Delivery)),c(552L,641L)))

weights <- static_cfg$get_calibration_weights()
wk <- function(y, s) paste(y, s, sep = "/")
key_w <- wk(weights$Year, weights$Region)
key_i <- wk(iso$Year, as.character(iso$Region))
iw <- match(key_i, key_w)
n_cell <- as.integer(table(factor(key_i, levels = key_w)))
pop <- as.numeric(weights$Births)
wg <- as.numeric(weights$Weight)
frac <- function(x) x / ave(x, weights$Year, FUN = sum)
stopifnot(!anyNA(iw), length(key_w) == 40L, !anyDuplicated(key_w), all(wg > 0),
          all(n_cell == as.integer(weights$Isolates)),
          max(abs(wg - frac(pop) / frac(n_cell))) < 1e-8)
iso$wg <- wg[iw]


options(survey.lonely.psu="fail",survey.adjust.domain.lonely=FALSE)
design <- survey::svydesign(ids=~Site,strata=~Region,weights=~wg,nest=TRUE,data=iso)
stopifnot(survey::degf(design)==26)
domain <- function(d, keep) {
  stopifnot(length(keep)==nrow(d$variables))
  use_rows <- !is.na(keep)&keep
  subset(d,use_rows)
}
make_y <- function(d, y) {
  stopifnot(length(y)==nrow(d$variables))
  d$variables$Endpoint <- y
  d
}
diag_rows <- list()
attempt <- function(label, expr) {
  warn <- character()
  fit <- tryCatch(withCallingHandlers(expr,warning=function(w) {
    warn <<- c(warn,conditionMessage(w))
    invokeRestart("muffleWarning")
  }),error=function(e) {warn <<- c(warn,conditionMessage(e));NULL})
  if (length(warn)) diag_rows[[length(diag_rows)+1L]] <<- data.frame(Analysis=label,Detail=paste(unique(warn),collapse=" | "))
  fit
}
estimate <- function(d,label) {
  z <- d$variables
  n <- nrow(z)
  e <- sum(z$Endpoint)
  df <- survey::degf(d)
  ans <- data.frame(N=n,Events=e,Unweighted=if(n) e/n else NA_real_,
    Weighted=if(n) weighted.mean(z$Endpoint,z$wg) else NA_real_,
    SE=NA_real_,Low=NA_real_,High=NA_real_,DF=df,PSU=length(unique(z$Site)),
    Method="Taylor logit-t; actual domain design df",CI_status="No evaluable results")
  if (!n) return(ans)
  m <- attempt(label,survey::svymean(~Endpoint,d))
  if (!is.null(m)) {
    stopifnot(abs(as.numeric(coef(m))-ans$Weighted)<1e-10)
    ans$SE <- as.numeric(survey::SE(m))
  }
  logit <- qlogis(ans$Weighted)
  ans$CI_status <- "The Taylor logit-t interval is undefined at an exact boundary"
  if (!is.finite(logit)) return(ans)
  ans$CI_status <- "Insufficient design df or nonpositive variance"
  if (df<=0 || !is.finite(ans$SE) || ans$SE<=0) return(ans)
  scores <- setNames(numeric(nrow(site)),site[["S4 center"]])
  totals <- tapply(z$wg*(z$Endpoint-ans$Weighted)/sum(z$wg),z$Site,sum)
  scores[names(totals)] <- totals
  manual_var <- sum(vapply(split(scores,site[["S4 stratum"]]),function(x) length(x)*var(x),0))
  stopifnot(abs(manual_var-ans$SE^2)<1e-10)
  manual <- plogis(logit+c(-1,1)*qt(.975,df)*ans$SE/(ans$Weighted*(1-ans$Weighted)))
  ans$Low <- manual[1]; ans$High <- manual[2]
  ans$CI_status <- if (min(e,n-e)<5) "Estimated; fewer than five events or non-events" else "Estimated"
  ans
}

fit_model <- function(d, x, label, min_event) {
  a <- d$variables
  e <- sum(a$Endpoint)
  df <- survey::degf(d)
  n <- nrow(a)
  out <- data.frame(N=n,Events=e,PSU=length(unique(a$Site)),Design_DF=df,Model_DF=NA_real_,
    Beta=NA_real_,SE=NA_real_,Statistic=NA_real_,P=NA_real_,PR=NA_real_,Low=NA_real_,High=NA_real_,
    PR5=NA_real_,Low5=NA_real_,High5=NA_real_,
    Reference="Student t, residual model df",
    Status="Sparse or no outcome/covariate variation")
  if(n<10L || e<min_event || e==n || length(unique(a[[x]]))<2 ||
      (x=="Year" && n-e<min_event)) return(out)
  f <- reformulate(x,"Endpoint")
  z <- attempt(label,survey::svyglm(f,d,family=quasipoisson(link="log")))
  out$Status <- "Fit failed or unstable"
  if(is.null(z) || !isTRUE(z$converged)) return(out)
  b <- unname(coef(z)[2])
  se <- sqrt(vcov(z)[2,2])
  mdf <- z$df.residual
  if(!is.finite(b)||!is.finite(se)||abs(b)>=5||se<=0||se>=5||is.na(mdf)||mdf<=0) return(out)
  t <- b/se
  p <- 2*pt(-abs(t),mdf)
  lim <- b+c(-1,1)*qt(.975,mdf)*se
  model_summary <- coef(summary(z, df.resid=mdf))[2,]
  stopifnot(isTRUE(all.equal(unname(model_summary),c(b,se,t,p),tolerance=1e-10)),
    (p<.05)==(lim[1]>0 || lim[2]<0))
  out[1,c("Model_DF","Beta","SE","Statistic","P","PR","Low","High","PR5","Low5","High5",
    "Status")] <- list(mdf,b,se,t,p,
    exp(b),exp(lim[1]),exp(lim[2]),exp(5*b),exp(5*lim[1]),exp(5*lim[2]),"Estimated")
  out
}
unit <- survey::svydesign(ids=~Site,strata=~Region,weights=~1,nest=TRUE,data=iso)
