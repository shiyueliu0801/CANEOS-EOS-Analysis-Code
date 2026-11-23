#!/usr/bin/env Rscript
# CANEOS S1 Fig, 4B-6 legend-only refinement (v1.2.2).
# Usage: Rscript --vanilla plot_S1_Fig.R [INPUT_DIR [OUTPUT_TIFF]]
suppressPackageStartupMessages({library(sf); library(ggplot2); library(cowplot)})

aea <- st_crs("+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 +datum=WGS84 +units=m")
wgs <- st_crs(4326)
blue <- "#00A9E6"
orange <- "#FFAA00"
edge <- "#233440"
reef <- "#416C80"
px <- c(3509L,2481L)
inches <- px/300
bubble_max_mm <- 13
inset_limits <- c(104.4,1.50,122.45,25.60)

project <- function(x) st_transform(x,aea)
box <- function(v,crs) st_bbox(setNames(v,c("xmin","ymin","xmax","ymax")),crs=crs)
flat <- function(g) st_set_crs(g,NA)
map_theme <- theme_void() + theme(plot.margin=margin(0,0,0,0),
    panel.background=element_rect(fill="white",colour=NA))
coord_map <- function(b,crs) coord_sf(crs=crs,xlim=b[c(1,3)],ylim=b[c(2,4)],
                                      expand=FALSE,datum=NA)
island_points <- function(d,symbol) d[d$symbol==symbol,]
# Rendering derivatives only: the GeoPackage and its five reference points
# remain byte-identical to v1.1. No buffer, snapping or vertex editing.
prepare_display <- function(x) {
  hk <- x$sar[x$sar$GU_A3=="HKG",]
  x$macau <- x$sar[x$sar$GU_A3=="MAC",]
  stopifnot(nrow(hk)==1L,nrow(x$macau)==1L)
  country <- st_union(x$main_china)
  x$render_country <- st_union(c(country,st_geometry(hk)))
  shared <- st_intersection(st_boundary(st_geometry(hk)),st_boundary(country))
  x$hk_land_border <- suppressWarnings(st_collection_extract(shared,"LINESTRING"))
  stopifnot(sum(as.numeric(st_length(x$hk_land_border)))>0,
            abs(as.numeric(st_area(x$render_country))-sum(as.numeric(st_area(country)))-
                sum(as.numeric(st_area(hk))))<1)
  refs <- x$east$refs
  x$east$major <- refs[refs$label %in% c("Diaoyu Dao","Huangwei Yu","Chiwei Yu"),]
  pair <- refs[match(c("Beixiao Dao","Nanxiao Dao"),refs$label),]
  stopifnot(nrow(pair)==2L,!anyNA(pair$label))
  xy <- st_coordinates(project(pair))
  # A single small-scale GROUP SYMBOL at the projected reference midpoint.
  # This does not relocate either island or create a connecting land polygon.
  x$east$pair_symbol <- st_sf(members="Beixiao Dao | Nanxiao Dao",geometry=
    st_sfc(st_point(colMeans(xy)),crs=aea))
  x
}

draw_main_map <- function(x,b) {
  ggplot() +
    geom_sf(data=x$render_country,fill=blue,colour=edge,linewidth=.25) +
    geom_sf(data=x$macau,fill=blue,colour="#294D5D",linewidth=.12) +
    geom_sf(data=x$hk_land_border,colour="#294D5D",linewidth=.17) +
    geom_sf(data=x$province$geometry,colour="#294D5D",linewidth=.17) +
    geom_sf(data=x$east$major,shape=21,fill=blue,colour=edge,size=.40,stroke=.06) +
    geom_sf(data=x$east$pair_symbol,shape=21,fill=blue,colour=edge,size=.20,stroke=.05) +
    geom_sf(data=x$mar$east,colour=edge,linewidth=.44) +
    geom_sf(data=x$bubbles,aes(size=infant_n),shape=21,fill=orange,colour="#6E512E",stroke=.25,alpha=.98) +
    # Extend the reference scale to 200 without changing observed bubble sizes.
    scale_size_area(max_size=bubble_max_mm*sqrt(200/max(x$bubbles$infant_n)),
                    limits=c(0,200),breaks=c(20,50,100,200),
                    name="Infants") +
    coord_map(b,aea) + map_theme
}

draw_scs_inset <- function(x,b,crs) {
  ggplot() +
    geom_sf(data=x$context_land,fill="#EDF1F2",colour=NA) +
    geom_sf(data=x$context_coast,colour="#A3B3BD",linewidth=.15) +
    geom_sf(data=x$render_country,fill=blue,colour=edge,linewidth=.19) +
    geom_sf(data=x$macau,fill=blue,colour="#294D5D",linewidth=.12) +
    geom_sf(data=x$hk_land_border,colour="#294D5D",linewidth=.14) +
    geom_sf(data=x$province$geometry,colour="#294D5D",linewidth=.14) +
    geom_sf(data=x$south$reef,fill=NA,colour=reef,linewidth=.20) +
    geom_sf(data=x$south$land,fill=blue,colour=reef,linewidth=.16) +
    geom_sf(data=island_points(x$south$points,"land"),shape=16,colour=reef,size=.80,stroke=0) +
    geom_sf(data=island_points(x$south$points,"water"),shape=23,fill="white",colour=reef,size=1.05,stroke=.24) +
    geom_sf(data=island_points(x$south$points,"small water"),shape=23,fill="white",colour=reef,size=.40,stroke=.14) +
    geom_sf(data=x$mar$south,colour=edge,linewidth=.39) +
    coord_map(b,crs) + map_theme
}

contains <- function(g,b) {
  region <- flat(st_as_sfc(b))
  all(lengths(st_covered_by(flat(g),region))==1L)
}
fit <- function(b,slot,unit_ratio,align=c(.5,.5)) {
  ratio <- (b[3]-b[1])/(b[4]-b[2])*unit_ratio/(px[1]/px[2])
  h <- min(slot[4],slot[3]/ratio)
  w <- h*ratio
  setNames(c(slot[1]+(slot[3]-w)*align[1],slot[2]+(slot[4]-h)*align[2],w,h),c("x","y","w","h"))
}
place <- function(p,b) draw_plot(p,x=b["x"],y=b["y"],width=b["w"],height=b["h"])
rect_grob <- function(col,lwd) grid::rectGrob(gp=grid::gpar(fill=NA,col=col,lwd=lwd))

composition <- function(x) {
  # Preserve every China POV polygon in the main panel at or north of Hainan.
  main <- st_bbox(x$source_china)
  main[c(1,2)] <- main[c(1,2)]-c(43000,26000)
  main[c(3,4)] <- main[c(3,4)]+c(47000,38000)
  stopifnot(contains(x$main_china,main),contains(st_geometry(project(x$mar$east)),main),
            contains(st_geometry(project(x$east$refs)),main),contains(st_geometry(x$bubbles),main))
  inset_wgs <- box(inset_limits,wgs)
  icrs <- wgs
  ib <- inset_wgs
  required <- list(st_geometry(st_transform(x$mar$south,icrs)),st_geometry(st_transform(x$south$points,icrs)),
                   st_geometry(st_transform(x$south$reef,icrs)),st_geometry(st_transform(x$south$land,icrs)))
  stopifnot(all(vapply(required,contains,logical(1),b=ib)))
  p <- draw_main_map(x,main)
  markers <- tail(ggplot_build(p)$data,1)[[1]]
  stopifnot(nrow(markers)==20L,all(abs(markers$size^2/x$bubbles$infant_n-bubble_max_mm^2/max(x$bubbles$infant_n))<1e-10))
  u <- ggplotGrob(p+guides(size=guide_legend(title="Unique infants",title.position="top",direction="vertical"))+
      theme(legend.position="right",legend.title=element_text(size=10.5,face="bold",colour=edge),
            legend.text=element_text(size=9.4,colour=edge),legend.key.height=grid::unit(2,"mm"),
            legend.key.width=grid::unit(bubble_max_mm+1.5,"mm"),legend.margin=margin(1,2,1,2,"mm")))
  legend <- u$grobs[[which(u$layout$name=="guide-box-right")]]
  # Keep the previous title's layout footprint so shorter text cannot move keys.
  guide_i <- which(legend$layout$name=="guides")
  title_i <- which(legend$grobs[[guide_i]]$layout$name=="title")
  title_grob <- legend$grobs[[guide_i]]$grobs[[title_i]]
  stopifnot(length(guide_i)==1L,length(title_i)==1L,
            length(title_grob$children)==1L,title_grob$children[[1]]$label=="Unique infants")
  legend$grobs[[guide_i]]$grobs[[title_i]]$children[[1]]$label <- "Infants"
  mb <- fit(main,c(.018,.021,.805,.957),1,align=c(0,.5))
  scale <- cos(mean(ib[c(2,4)])*pi/180)
  sb <- fit(ib,c(.781,.042,.197,.383),scale,align=c(1,0))
  north <- grid::grobTree(grid::textGrob("N",x=.5,y=.14,gp=grid::gpar(fontsize=10.5,col=edge)),
    grid::polygonGrob(c(.5,.17,.5),c(.93,.26,.51),gp=grid::gpar(fill="white",col=edge,lwd=.8)),
    grid::polygonGrob(c(.5,.83,.5),c(.93,.26,.51),gp=grid::gpar(fill=edge,col=edge,lwd=.8)))
  figure <- ggdraw() +
    place(p+theme(legend.position="none"),mb) +
    place(draw_scs_inset(x,ib,icrs),sb) +
    draw_grob(rect_grob(edge,.7),x=sb["x"],y=sb["y"],width=sb["w"],height=sb["h"]) +
    draw_grob(north,x=.939,y=.866,width=.033,height=.087) +
    draw_grob(grid::rectGrob(gp=grid::gpar(fill="white",col=NA)),x=.052,y=.046,width=.149,height=.224) +
    draw_grob(legend,x=.059,y=.048,width=.135,height=.220) +
    draw_grob(rect_grob("#17232C",.77),x=.010,y=.013,width=.98,height=.974)
  panel <- function(name,b,slot,crs) data.frame(panel=name,t(slot),xmin=b[1],ymin=b[2],xmax=b[3],ymax=b[4],crs=crs,row.names=NULL)
  list(fig=figure,main=main,inset=ib,inset_crs=icrs,main_box=mb,inset_box=sb,
       panels=rbind(panel("main",main,mb,"Albers"),panel("inset",ib,sb,"wgs84")))
}



build_provincial_counts <- function(cohort_csv, hosp_prov_csv, anchors_df) {
  raw <- readr::read_csv(cohort_csv, col_types = readr::cols(.default = readr::col_character()),
                         locale = readr::locale(encoding = "UTF-8"), na = character(),
                         trim_ws = FALSE, show_col_types = FALSE, name_repair = "minimal")
  id_clean <- trimws(raw[["ID"]])
  primary <- which(nzchar(id_clean))
  extra   <- which(!nzchar(id_clean))
  hosp    <- trimws(raw[["医院名称"]])
  hosp[extra] <- hosp[extra - 1L]
  hmap <- readr::read_csv(hosp_prov_csv, col_types = readr::cols(.default = readr::col_character()),
                          locale = readr::locale(encoding = "UTF-8"), show_col_types = FALSE)
  idx <- match(hosp, trimws(hmap$Hospital))
  if (anyNA(idx)) stop("Hospital(s) missing from hospital_province_map: ", paste(unique(hosp[is.na(idx)]), collapse = ", "))
  iso_df <- data.frame(
    Hospital      = hosp,
    is_primary    = seq_len(nrow(raw)) %in% primary,
    province_cn   = trimws(hmap$province_cn[idx]),
    province_en   = trimws(hmap$province_en[idx]),
    ne_adm1_code  = trimws(hmap$ne_adm1_code[idx]),
    ne_iso_3166_2 = trimws(hmap$ne_iso_3166_2[idx]),
    stringsAsFactors = FALSE
  )
  keys <- unique(iso_df[c("province_cn", "province_en", "ne_adm1_code", "ne_iso_3166_2")])
  keys <- keys[order(keys$province_cn), , drop = FALSE]
  res <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    m <- iso_df$ne_adm1_code == keys$ne_adm1_code[i]
    data.frame(
      province_cn   = keys$province_cn[i],
      province_en   = keys$province_en[i],
      ne_adm1_code  = keys$ne_adm1_code[i],
      ne_iso_3166_2 = keys$ne_iso_3166_2[i],
      hospital_n    = as.integer(length(unique(iso_df$Hospital[m]))),
      infant_n      = as.integer(sum(iso_df$is_primary[m])),
      isolate_n     = as.integer(sum(m)),
      stringsAsFactors = FALSE
    )
  }))
  ai <- match(res$ne_adm1_code, anchors_df$ne_adm1_code)
  stopifnot(!anyNA(ai))
  res$x_m <- as.numeric(anchors_df$x_m[ai])
  res$y_m <- as.numeric(anchors_df$y_m[ai])
  rownames(res) <- NULL
  res
}

read_plot_data <- function(gpkg, provincial_df) {
  if (!file.exists(gpkg)) stop("Missing S1_geography.gpkg: ", gpkg)
  d <- provincial_df
  fields <- c("province_cn","province_en","ne_adm1_code","ne_iso_3166_2",
              "hospital_n","infant_n","isolate_n","x_m","y_m")
  stopifnot(identical(names(d),fields),nrow(d)==20L,!anyNA(d),
            !anyDuplicated(d$ne_adm1_code),!anyDuplicated(d$ne_iso_3166_2),
            !anyDuplicated(d$province_cn),!anyDuplicated(d$province_en),
            all(d$hospital_n>0L),all(d$infant_n>0L),all(d$isolate_n>=d$infant_n),
            sum(d$hospital_n)==30L,sum(d$infant_n)==1189L,sum(d$isolate_n)==1196L,
            identical(d$infant_n[d$ne_iso_3166_2=="CN-YN"],79L),
            all(is.finite(d$x_m)),all(is.finite(d$y_m)),!anyDuplicated(d[c("x_m","y_m")]))
  bubbles <- st_as_sf(d,coords=c("x_m","y_m"),remove=FALSE,crs=aea)
  layers <- c("china_pov_nearshore","china_pov_display","china_sar","province_internal_lines","inset_land","inset_coast",
              "maritime_lines","east_island_points","south_reef_outlines","south_island_polygons","south_feature_points")
  stopifnot(setequal(st_layers(gpkg)$name,layers))
  g <- setNames(lapply(layers,function(n) {
    s <- st_read(gpkg,layer=n,quiet=TRUE)
    projected <- n %in% c("china_pov_nearshore","china_pov_display","china_sar","province_internal_lines",
                         "south_reef_outlines","south_island_polygons")
    crs <- if (projected) aea else wgs
    stopifnot(nrow(s)>0L,!any(st_is_empty(s)),st_crs(s)==crs)
    st_set_crs(s,crs)
  }),layers)
  m <- g$maritime_lines
  stopifnot(nrow(m)==10L,identical(m$segment_id,sprintf("SU%02d",1:10)),
            identical(m$segment_id[m$panel=="main"],"SU04"),
            identical(m$panel,ifelse(m$segment_id=="SU04","main","inset")),
            nrow(g$east_island_points)==5L,nrow(g$china_pov_nearshore)==74L,nrow(g$china_pov_display)==71L,nrow(g$province_internal_lines)==75L,
            nrow(g$south_reef_outlines)==175L,nrow(g$south_island_polygons)==7L,
            identical(as.integer(table(factor(g$south_feature_points$symbol,
                      levels=c("land","water","small water")))),c(2L,2L,2L)))
  x <- list(source_china=st_geometry(g$china_pov_nearshore),main_china=st_geometry(g$china_pov_display),
       inset_china=st_geometry(g$china_pov_display),
       sar=g$china_sar,province=list(geometry=st_geometry(g$province_internal_lines)),bubbles=bubbles,
       east=list(refs=g$east_island_points),mar=list(east=m[m$panel=="main",],south=m[m$panel=="inset",]),
       context_land=g$inset_land,context_coast=g$inset_coast,south=list(reef=g$south_reef_outlines,
       land=g$south_island_polygons,points=g$south_feature_points))
  prepare_display(x)
}

render_s1_figure <- function(cohort_csv, hosp_prov_csv, anchors_df, gpkg_path, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sf_use_s2(FALSE)
  prov_df <- build_provincial_counts(cohort_csv, hosp_prov_csv, anchors_df)
  readr::write_csv(prov_df, file.path(out_dir, "S1_Provincial_Counts_Dynamic.csv"))
  x <- read_plot_data(gpkg_path, prov_df)
  p <- composition(x)
  out_tiff <- file.path(out_dir, "S1_Fig_300dpi.tiff")
  out_png  <- file.path(out_dir, "S1_Fig_300dpi.png")
  out_pdf  <- file.path(out_dir, "S1_Fig_vector.pdf")
  ggsave(out_tiff, p$fig, device = ragg::agg_tiff, width = inches[1], height = inches[2], dpi = 300,
         compression = "lzw", bg = "white")
  ggsave(out_png,  p$fig, device = ragg::agg_png,  width = inches[1], height = inches[2], dpi = 300, bg = "white")
  ggsave(out_pdf,  p$fig, device = grDevices::cairo_pdf, width = inches[1], height = inches[2], bg = "white")
  message("[PASS] Rendered S1 Fig: ", out_tiff, " (3509 x 2481 at 300 dpi; 20 provinces, 30 hospitals, 1189 infants, 1196 isolates)")
  invisible(list(tiff = out_tiff, pdf = out_pdf, counts = prov_df))
}
