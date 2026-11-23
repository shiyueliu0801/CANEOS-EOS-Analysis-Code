if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
#!/usr/bin/env Rscript
# CANEOS S1 Fig: rebuild the figure's 4B-5 geometry from pinned source components.
# Usage: Rscript --vanilla prepare_geography.R RAW_DIR [OUTPUT_GPKG]
# Upstream versions, direct download URLs and terms: see README.md.
suppressPackageStartupMessages(library(sf))
aea <- st_crs("+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 +datum=WGS84 +units=m")
wgs <- st_crs(4326)
project <- function(x) st_transform(x,aea)
km <- function(x) as.numeric(st_length(x))/1000
area <- function(x) as.numeric(st_area(x))/1e6
box <- function(v,crs) st_bbox(setNames(v,c("xmin","ymin","xmax","ymax")),crs=crs)
surface <- function(x) st_point_on_surface(st_geometry(x))
points_wgs <- function(d) st_as_sf(d,coords=c("lon","lat"),crs=wgs,remove=FALSE)
flat <- function(g) st_set_crs(g,NA)
window_intersects <- function(g,b) lengths(st_intersects(flat(g),flat(st_as_sfc(b))))>0L
coords_wgs <- function(g) st_coordinates(st_transform(g,wgs))[,1:2,drop=FALSE]
extent_rows <- function(g) {
  d <- as.data.frame(t(vapply(g,function(x) as.numeric(st_bbox(x)),numeric(4))))
  setNames(d,c("west","south","east","north"))
}

raw_hashes <- c(
  "china.dbf" = "15ad6a0580212a5630352c5c89d1f0e2c50a4cbed85b5f8399b18b3daf3c8170",
  "china.prj" = "a02a27b1d1982c8516d83398e85a3c8b1aef1713c13ef4d84d7bde17430c07c4",
  "china.shp" = "41875d9d8a7d994ab7326fdd6b366f10dfdb3fdcd208b894604a258a789a6b47",
  "china.shx" = "c5b9d3a66f3791f9419a2d4d1b983b90af1340de3e959e29f9c17bd418dd1087",
  "china_nine_dotted_line.dbf" = "8cb8f18ff1ab00304a847c202d60293e0966662b90ea3d77ab29e3ccf5e3fd22",
  "china_nine_dotted_line.prj" = "a02a27b1d1982c8516d83398e85a3c8b1aef1713c13ef4d84d7bde17430c07c4",
  "china_nine_dotted_line.shp" = "c5ce7a0a1d2d07e1b327fe2ae5f595c5a3cdd99aa6eb68eb6a32259139e6a538",
  "china_nine_dotted_line.shx" = "64f83a53eb42db8ee708ef8ded4b5393c3f9e4675efdf5d7cf67ab175fbabf14",
  "ne_10m_admin_0_countries_chn.dbf" = "c8856a38848583056b9365eb01b1f142233e959c3213d17c40b3286ab3e11abb",
  "ne_10m_admin_0_countries_chn.prj" = "a02a27b1d1982c8516d83398e85a3c8b1aef1713c13ef4d84d7bde17430c07c4",
  "ne_10m_admin_0_countries_chn.shp" = "2b5dc7fb48c76f8aa0a3007d73eb816ce9b0a2b1f4a10c5339e1d73dee52df8e",
  "ne_10m_admin_0_countries_chn.shx" = "6892b529f63f14a01d10709293351095846a39493f13efccce8395e774049946",
  "ne_10m_admin_1_states_provinces.cpg" = "3ad3031f5503a4404af825262ee8232cc04d4ea6683d42c5dd0a2f2a27ac9824",
  "ne_10m_admin_1_states_provinces.dbf" = "445a8a9bea889634faf0af18081830df0b05b8471fc6af8dc42aecdd7a71bba1",
  "ne_10m_admin_1_states_provinces.prj" = "a02a27b1d1982c8516d83398e85a3c8b1aef1713c13ef4d84d7bde17430c07c4",
  "ne_10m_admin_1_states_provinces.shp" = "c6f5c8b4b1320d9417033762419c6df1eb423989cd880fba78ea0b1e3522cbe4",
  "ne_10m_admin_1_states_provinces.shx" = "37a9e2bc79ed31d3bdea3cb62d928f77281a1c88d645cd33430231c75dbcf350",
  "ne_10m_admin_1_states_provinces_lines.cpg" = "3ad3031f5503a4404af825262ee8232cc04d4ea6683d42c5dd0a2f2a27ac9824",
  "ne_10m_admin_1_states_provinces_lines.dbf" = "b524fd1baf59386cd4bca0d6df642ddb02be9cdc4d69f1f1d109885acf21d9ed",
  "ne_10m_admin_1_states_provinces_lines.prj" = "a02a27b1d1982c8516d83398e85a3c8b1aef1713c13ef4d84d7bde17430c07c4",
  "ne_10m_admin_1_states_provinces_lines.shp" = "c89b8545fdb03a764df9de2efeadf0d6735169d94b48a315b5e233780e073ab5",
  "ne_10m_admin_1_states_provinces_lines.shx" = "ba65a7072e06b5454bca8c0a527978333571f3ce45094f006a60f18f8509db10",
  "ne_10m_coastline.cpg" = "3ad3031f5503a4404af825262ee8232cc04d4ea6683d42c5dd0a2f2a27ac9824",
  "ne_10m_coastline.dbf" = "9ccc214342fe400bf8c7d91d7a5b276b0457b0ada03e8d4be16ac5ba13037f3b",
  "ne_10m_coastline.prj" = "a02a27b1d1982c8516d83398e85a3c8b1aef1713c13ef4d84d7bde17430c07c4",
  "ne_10m_coastline.shp" = "459a4a97c09db19aadf5244026612de9d43748be27f83a360242b99f7fabb3c1",
  "ne_10m_coastline.shx" = "f873afee7f56779ce52253f740ec251c2f12244aea911dc40f0a85d75de8d5f2",
  "ne_10m_land.cpg" = "3ad3031f5503a4404af825262ee8232cc04d4ea6683d42c5dd0a2f2a27ac9824",
  "ne_10m_land.dbf" = "e422520372e65c6fb8220888580e5331d04699bcd71c27e53ad42e241a531cbd",
  "ne_10m_land.prj" = "a02a27b1d1982c8516d83398e85a3c8b1aef1713c13ef4d84d7bde17430c07c4",
  "ne_10m_land.shp" = "4cad3a49bc75c1a4c2f3d7efae04f2f8e63151c96764b2658effabf524331fa6",
  "ne_10m_land.shx" = "f3b02a190c3020eeafd0f48de982619826b57ecd853cbc68e4fc049652891a04"
)
shape_files <- c(
  countries="ne_10m_admin_0_countries_chn.shp",admin="ne_10m_admin_1_states_provinces.shp",
  lines="ne_10m_admin_1_states_provinces_lines.shp",land="ne_10m_land.shp",coast="ne_10m_coastline.shp",
  islands="china.shp",maritime="china_nine_dotted_line.shp")

check_sources <- function(raw_dir) {
  p <- file.path(raw_dir, names(raw_hashes))
  missing <- !file.exists(p)
  if (any(missing)) stop("Missing source components: ", paste(names(raw_hashes)[missing], collapse=", "))
  found <- vapply(p, function(f) digest::digest(file=f, algo="sha256"), character(1), USE.NAMES=FALSE)
  bad <- found != unname(raw_hashes)
  if (any(bad)) stop("Wrong source version or changed component: ", paste(names(raw_hashes)[bad], collapse=", "))
  message("SHA-256 verified for all ", length(p), " source components.")
}

# Complete source tables previously supplied as separate CSVs. The official index
# is diagnostic only. IDs 13 and 14 use 44.4 minutes N per the official listing.
read_embedded <- function(s) read.csv(text=s, check.names=FALSE, na.strings="")
south_rules <- read_embedded(r"{fename,expected_rows,new_to_audit,mode,min_area,lon,lat,reference_url,note
Dongsha Qundao,7,FALSE,outline,1,,,https://www.taiwan.nps.gov.tw/home/zh-tw/quarterly/7918/6621,Atoll outlines represent mixed reef and lagoon. Tiny source records omitted at inset scale.
Xisha Qundao,18,FALSE,outline,1,,,https://scsio.cas.cn/xisha/zyts/202503/t20250317_7558785.html,The provider uses a group name for reef and island records.
Zhongsha Qundao,16,FALSE,outline,2,,,https://www.ey.gov.tw/state/4447F4A951A1EC45/094b1d53-de8d-4393-bde6-ab092969cce4,Submarine features receive outlines with no dry-land fill.
Nansha Qundao,162,FALSE,outline,2,,,https://www.mfa.gov.cn/web/wjb_673085/zzjg_673183/bjhysws_674671/bhflfg/dtdmxgfl/202303/P020230313584749968923.pdf,Group-level source class does not identify individual dry land.
Huangyan Dao,3,FALSE,outline,1,,,https://www.gsjn.gov.cn/zfxxgk/zfbmxxgk/zrzyj/fdzdgknr/zrzyxx/chdlxxgl/art/2023/art_7e0a51d9da3940eab95575cf98963df0.html,Source polygons include reef-scale footprint.
Zengmu Ansha,3,FALSE,water_point,,112.283333,3.966667,,Inherited approximate shoal reference. Other same-name positions remain unconfirmed.
Yongxing Dao,1,TRUE,reference_point,,112.3376944444,16.8346722222,https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002422,The source footprint is far from the island. Plot the published geographic reference point.
Ganquan Dao,11,TRUE,reference_point,,111.585,16.506111,https://cil.nus.edu.sg/wp-content/uploads/2017/08/Table-of-Shallow-Features-Names.pdf,Only the independently located island is drawn. Other same-name records may be reefs or different entities.
Zhongjian Dao,1,TRUE,island_outline,,,,https://www.mfa.gov.cn/wjb_673085/zzjg_673183/bjhysws_674671/bhflfg/dtdmxgfl/200904/t20090409_9866755.shtml,Source location and extent agree at draft scale with published baseline point cluster.
Bei Jiao,2,TRUE,water_point,,111.516667,17.116667,https://www.mfa.gov.cn/wjb_673085/zzjg_673183/bjhysws_674671/bhflfg/dtdmxgfl/200904/t20090409_9866755.shtml,Reef complex represented once at inset scale. Nearby smaller source component remains in QC.
Wan an Tan,1,TRUE,hold,,109.775,7.5,https://geodoi.ac.cn/DOIPaper/HTML/2018139402/201802139414.htm,Submerged bank with source and reference position discrepancy. Group-scale draft does not require its individual footprint.
Nan an Jiao,11,TRUE,hold,,,,https://geodoi.ac.cn/DOIPaper/HTML/2022139401/202201139416.html,Repeated name on widely separated features. Individual identity or reef footprint unresolved.
Beiwei Tan,1,TRUE,small_water_point,,116.0032,21.0785,https://www.mee.gov.cn/ywdt/gs/wqgs_1/201905/W020190515614462913511.pdf,Submerged bank; use a smaller open symbol to keep the two neighbouring banks visually separate.
Nanwei Tan,1,TRUE,small_water_point,,115.95,20.966667,https://www.mee.gov.cn/ywdt/gs/wqgs_1/201905/W020190515614462913511.pdf,Submerged bank; use a smaller open symbol to keep the two neighbouring banks visually separate.
Qizhou Liedao,6,TRUE,island_outline,,,,https://oeoc.jou.edu.cn/info/1121/1016.htm,Compact named nearshore group agrees with the published planning locator. Outlines retained only in the inset.}")
east_reference <- read_embedded(r"{fename,label,lon,lat,official_area_km2,reference_url
Tiaoyu Tao,Diaoyu Dao,123.4733333333,25.7433333333,3.91,https://www.mfa.gov.cn/diaoyudao/chn/xwdt/201209/t20120915_8523866.htm
Huangwei Yu,Huangwei Yu,123.6816666667,25.9233333333,0.91,https://www.mfa.gov.cn/diaoyudao/chn/xwdt/201209/t20120915_8523866.htm
Chihwei Yu,Chiwei Yu,124.5583333333,25.9216666667,0.065,https://www.mfa.gov.cn/diaoyudao/chn/xwdt/201209/t20120915_8523866.htm
Peihsiao Tao,Beixiao Dao,123.5416666667,25.73,,https://www.mfa.gov.cn/diaoyudao/chn/xwdt/201209/t20120915_8523866.htm
Nanhsiao Tao,Nanxiao Dao,123.55,25.7233333333,,https://www.mfa.gov.cn/diaoyudao/chn/xwdt/201209/t20120915_8523866.htm}")
east_index <- read_embedded(r"{official_id,lon_deg,lon_min,lat_min_at_25_deg_north
1,123,28.4,44.6
2,123,29.2,45
3,123,27.6,44
4,123,28.6,44.4
5,123,29,44.4
6,123,29.3,44.4
7,123,29.3,44.4
8,123,29.3,44.4
9,123,29.3,44.4
10,123,29.3,44.4
11,123,29.4,44.5
12,123,29.4,44.5
13,123,29.4,44.4
14,123,29.4,44.4
15,123,29.5,44.5
16,123,40.9,55.4
17,123,40.7,55.8
18,123,40.7,55.6
19,123,40.7,55.6
20,123,41,55.7
21,123,41.1,55.8
22,123,41.1,55.6
23,123,41.1,55.6
24,123,41.1,55.6
25,123,41.2,55.6
26,123,41.2,55.6
27,123,40.6,55.5
28,123,40.6,55.5
29,123,40.5,55.5
30,123,40.6,55.4
31,123,40.6,55.3
32,123,40.7,55.3
33,123,40.7,55.2
34,123,41.4,55.3
35,123,41.3,55.6
36,123,41.2,55.2
37,124,33.5,55.3
38,124,33.5,55.5
39,124,33.7,55.5
40,124,33.5,55.4
41,124,33.6,55.4
42,124,33.3,55.3
43,124,33.4,55.3
44,124,33.4,55.3
45,124,33.4,55.3
46,124,33.3,55.3
47,124,33.2,55.2
48,123,32.5,43.8
49,123,32.5,43.9
50,123,32.5,43.9
51,123,32.7,43.7
52,123,33,43.4
53,123,32.7,43.5
54,123,32.6,43.5
55,123,32.7,43.4
56,123,32.6,43.4
57,123,32.6,43.4
58,123,33.4,43.4
59,123,33.4,43.3
60,123,33.3,43.3
61,123,33.4,43.3
62,123,33.4,43.2
63,123,33.3,43.3
64,123,32.6,46.9
65,123,32.6,46.8
66,123,32.6,46.8
67,123,32.4,46.8
68,123,32.5,46.8
69,123,34,45.3
70,123,30.4,44.1
71,123,30.4,44.1}")

inputs <- function(raw_dir) {
  p <- setNames(file.path(raw_dir, shape_files), names(shape_files))
  read_ne <- function(key) st_read(p[[key]],quiet=TRUE)
  a <- read_ne("admin")
  main_admin <- a[which(a$adm0_a3=="CHN" & a$geonunit=="China"),]
  stopifnot(nrow(main_admin)==31L)
  countries <- read_ne("countries")
  country <- project(countries[which(countries$ADMIN=="China" & countries$ADM0_A3_CN=="CHN"),])
  sar <- project(countries[which(countries$GU_A3 %in% c("HKG","MAC")),])
  stopifnot(nrow(country)==1L,nrow(sar)==2L,all(st_is_valid(country)))
  l <- read_ne("lines")
  z <- st_read(p[["islands"]],options="ENCODING=UTF-8",quiet=TRUE)
  mar <- st_read(p[["maritime"]],quiet=TRUE)
  stopifnot(nrow(z)==1367L,nrow(mar)==10L,
            all(c("\u6d77\u5357\u7701","\u5e7f\u4e1c\u7701","\u53f0\u6e7e\u7701") %in% z$OWNER), "\u6c38\u5174\u5c9b" %in% z$FCNAME)

  list(country=country, sar=sar, admin=main_admin, lines=l, land=read_ne("land"),
       coast=read_ne("coast"), islands=z, maritime=mar,
       south_rules=south_rules, east_reference=east_reference, east_index=east_index)
}

build_province_lines <- function(x) {
  i <- which(x$lines$ADM0_A3=="CHN")
  s <- project(x$lines[i,])
  s$source_row <- i
  names_ok <- s$NAME_L %in% x$admin$name & s$NAME_R %in% x$admin$name & s$NAME_L!=s$NAME_R
  names_ok[is.na(names_ok)] <- FALSE
  classes <- c("Admin-1 boundary","Admin-1 statistical boundary","Admin-1 statistical meta bounds")
  class_ok <- s$FEATURECLA %in% classes &
    (is.na(s$FCLASS_CN) | !grepl("Unrecognized|International|Claim|Disputed",s$FCLASS_CN,ignore.case=TRUE))
  g <- st_geometry(s)
  country <- st_geometry(x$country)
  clip <- function(y) {
    a <- st_intersection(y,country)
    if(any(st_geometry_type(a)=="GEOMETRYCOLLECTION")) a <- suppressWarnings(st_collection_extract(a,"LINESTRING"))
    a[which(!st_is_empty(a) & st_dimension(a)==1L)]
  }
  selection <- which(names_ok & class_ok)
  clipped <- lapply(selection,function(k) clip(g[k]))
  drawn <- do.call(c,clipped[lengths(clipped)>0L])
  n <- numeric(nrow(s))
  n[selection] <- vapply(clipped,function(y) sum(km(y)),numeric(1))
  stopifnot(nrow(s)==76L,sum(n>0)==75L,sum(km(st_difference(drawn,country)))<.002)
  xy <- coords_wgs(st_line_sample(st_cast(g[which(n==0)],"LINESTRING"),sample=.5) |> st_cast("POINT"))
  lookup <- st_is_within_distance(g[which(n==0)],g,dist=100)[[1]]
  nearest <- setdiff(lookup,which(n==0))
  fields <- c("source_row","ADM0_A3","FEATURECLA","FCLASS_CN","NAME_L","NAME_R","ADM0_NAME",
              "NAME","ADM1_NAME","TYPE","LABELRANK","MIN_ZOOM","NE_ID")
  fields <- intersect(fields,names(s))
  d <- data.frame(st_drop_geometry(s)[,fields],length_m=1000*km(g),drawn_m=1000*n,
    reason=ifelse(!names_ok,"Exclude: no two distinct named provinces",ifelse(!class_ok,"Exclude: source class/view",
      ifelse(n==0,"Exclude: outside China POV","Draw: approved line class and named sides, clipped to China POV"))))
  d$excluded_midpoint_lon <- d$excluded_midpoint_lat <- NA_real_
  d[d$drawn_m==0,c("excluded_midpoint_lon","excluded_midpoint_lat")] <- xy
  d$excluded_nearby_source_rows_within_100m <- NA_character_
  d$excluded_nearby_source_rows_within_100m[d$drawn_m==0] <- paste(i[nearest],collapse="|")
  list(geometry=drawn,audit=d)
}

describe_shapes <- function(z,idx) {
  g <- project(z[idx,])
  stopifnot(all(st_is_valid(g)))
  ll <- coords_wgs(surface(g))
  cent <- coords_wgs(st_centroid(st_geometry(g)))
  d <- data.frame(source_row=idx,st_drop_geometry(z[idx,]),geometry=as.character(st_geometry_type(g)),
    extent_rows(st_geometry(z[idx,])),centroid_lon=cent[,1],centroid_lat=cent[,2],
    surface_lon=ll[,1],surface_lat=ll[,2],area_km2=area(g))
  list(g=g,d=d)
}

overlap_pairs <- function(g,d) {
  g <- st_geometry(g)
  hits <- st_intersects(g)
  equal <- st_equals(g)
  pairs <- do.call(rbind,lapply(seq_along(hits),function(i) {
    j <- hits[[i]][hits[[i]]>i]
    if(!length(j)) return(NULL)
    data.frame(i=i,j=j,equals=j %in% equal[[i]],area_km2=vapply(j,function(k) {
      a <- st_intersection(g[i],g[k])
      if(length(a)) sum(area(a)) else 0
    },numeric(1)))
  }))
  if(is.null(pairs)) pairs <- data.frame(i=integer(),j=integer(),equals=logical(),area_km2=numeric())
  pairs <- pairs[pairs$equals | pairs$area_km2>1e-8,]
  transform(pairs,bound_a_id_a=d$BOUND_A_ID[i],bound_a_id_b=d$BOUND_A_ID[j],
    fename_a=d$FENAME[i],fename_b=d$FENAME[j],
    fraction_smaller=area_km2/pmin(d$area_km2[i],d$area_km2[j]))
}

build_south <- function(x) {
  z <- x$islands
  # Broad nearshore/sea window, then only specific named Guangdong/Hainan rows.
  coarse <- which(window_intersects(st_geometry(z),box(c(106,2,123,22),wgs)) &
    z$OWNER %in% c("\u6d77\u5357\u7701","\u5e7f\u4e1c\u7701") & !is.na(z$FCNAME) & z$FCNAME!=z$OWNER)
  zt <- describe_shapes(z,coarse)
  d <- zt$d
  take <- which(d$surface_lon>106 & d$surface_lon<123 & d$surface_lat>2 & d$surface_lat<22)
  d <- d[take,]
  g <- zt$g[take,]
  rule <- x$south_rules
  stopifnot(nrow(d)==244L,length(unique(d$FENAME))==15L,!anyDuplicated(rule$fename),
            setequal(d$FENAME,rule$fename))
  j <- match(d$FENAME,rule$fename)
  stopifnot(all(as.integer(table(factor(d$FENAME,levels=rule$fename)))==rule$expected_rows))
  d$new_to_audit <- rule$new_to_audit[j]
  stopifnot(sum(d$new_to_audit)==35L)
  pairs <- overlap_pairs(g,d)
  d$exact_duplicate_ids <- vapply(seq_len(nrow(d)),function(i) {
    paste(c(pairs$bound_a_id_b[pairs$i==i & pairs$equals],pairs$bound_a_id_a[pairs$j==i & pairs$equals]),collapse="|")
  },character(1))
  d$overlap_ids <- vapply(seq_len(nrow(d)),function(i) {
    paste(c(pairs$bound_a_id_b[pairs$i==i],pairs$bound_a_id_a[pairs$j==i]),collapse="|")
  },character(1))
  d$symbol <- "held"
  d$reason <- rule$note[j]
  d$figure_lon <- d$figure_lat <- NA_real_
  d$distance_to_rule_reference_km <- NA_real_
  d$reference_url <- rule$reference_url[j]
  plotted_ll <- data.frame(lon=numeric(),lat=numeric(),symbol=character(),id=integer())
  for(k in seq_len(nrow(rule))) {
    group <- which(j==k)
    mode <- rule$mode[k]
    if(!is.na(rule$lon[k])) {
      anchor <- project(points_wgs(rule[k,c("lon","lat")]))
      d$distance_to_rule_reference_km[group] <- as.numeric(st_distance(g[group,],anchor))/1000
    }
    if(mode=="hold") next
    if(mode=="outline") {
      use <- group[d$area_km2[group]>=rule$min_area[k]]
      # An exact duplicate never needs a second stroke.
      use <- setdiff(use,pairs$j[pairs$equals & pairs$i %in% use & pairs$j %in% use])
      d$symbol[use] <- "unfilled reef/archipelago outline"
      d$reason[use] <- "Visible-scale source outline; wet/dry footprint uncertain, no land fill"
      d$reason[group[d$area_km2[group]<rule$min_area[k]]] <- "Below scale threshold for this group; retained in QC"
      next
    }
    if(mode=="island_outline") {
      d$symbol[group] <- "island polygon"
      next
    }
    chosen <- if(is.na(rule$lon[k])) group[which.max(d$area_km2[group])] else
      group[which.min(d$distance_to_rule_reference_km[group])]
    if(mode=="reference_point") {
      d$symbol[chosen] <- "independent reference point"
      d$reason[chosen] <- paste(rule$note[k],"Reference position only; source polygon excluded.")
      xy <- c(rule$lon[k],rule$lat[k])
      sym <- "land"
    } else {
      d$symbol[chosen] <- "open reef/shoal point"
      d$reason[chosen] <- paste(rule$note[k],"Only one source surface point used.")
      xy <- c(d$surface_lon[chosen],d$surface_lat[chosen])
      sym <- if(mode=="small_water_point") "small water" else "water"
    }
    d[chosen,c("figure_lon","figure_lat")] <- matrix(xy,nrow=1L)
    plotted_ll <- rbind(plotted_ll,data.frame(lon=xy[1],lat=xy[2],symbol=sym,id=d$BOUND_A_ID[chosen]))
  }
  stopifnot(all(c("Dongsha Qundao","Xisha Qundao","Zhongsha Qundao","Nansha Qundao","Huangyan Dao","Zengmu Ansha") %in%
                d$FENAME[d$symbol!="held"]))
  list(audit=d,pairs=pairs,all=g,reef=g[d$symbol=="unfilled reef/archipelago outline",],
       land=g[d$symbol=="island polygon",],points=points_wgs(plotted_ll))
}

build_east <- function(x) {
  z <- x$islands
  b <- box(c(123.20,25.55,124.85,26.12),wgs)
  t <- describe_shapes(z,which(window_intersects(st_geometry(z),b)))
  d <- t$d
  g <- t$g
  stopifnot(nrow(d)==11L,setequal(d$BOUND_A_ID[d$FENAME=="Taiwan Sheng"],c(2141,2142,2143,2145)))
  refs <- points_wgs(x$east_reference)
  distance <- as.matrix(st_distance(g,project(refs)))/1000
  near <- max.col(-distance,ties.method="first")
  same <- match(d$FENAME,refs$fename)
  d$nearest_reference <- refs$label[near]
  d$nearest_reference_km <- distance[cbind(seq_len(nrow(d)),near)]
  index <- x$east_index
  stopifnot(identical(index$official_id,1:71))
  index$lon <- index$lon_deg+index$lon_min/60
  index$lat <- 25+index$lat_min_at_25_deg_north/60
  stopifnot(max(abs(index$lon[c(1,16,37,48,52)]-refs$lon))<1e-8,
            max(abs(index$lat[c(1,16,37,48,52)]-refs$lat))<1e-8)
  candidates <- as.matrix(st_distance(g,project(points_wgs(index))))/1000
  d$nearest_official_id_all_71 <- apply(candidates,1,function(a) paste(index$official_id[abs(a-min(a))<1e-7],collapse="|"))
  d$nearest_official_id_distance_km <- apply(candidates,1,min)
  d$named_reference <- refs$label[same]
  d$official_lon <- refs$lon[same]
  d$official_lat <- refs$lat[same]
  d$official_area_km2 <- refs$official_area_km2[same]
  d$polygon_to_named_reference_km <- distance[cbind(seq_len(nrow(d)),same)]
  d$centre_to_named_reference_km <- as.matrix(st_distance(st_centroid(st_geometry(g)),project(refs)))[
    cbind(seq_len(nrow(d)),same)]/1000
  d$area_vs_official <- d$area_km2/d$official_area_km2
  d$decision <- ifelse(is.na(same),"Held: generic Taiwan Sheng attribute; identity unconfirmed",
    "Source polygon held: use named reference point instead of uncertain island outline")
  d$reference_url <- refs$reference_url[near]
  pairs <- overlap_pairs(g,d)
  both <- do.call(rbind,lapply(which(table(d$FENAME)==2L),function(k) {
    name <- names(table(d$FENAME))[k]
    i <- which(d$FENAME==name)
    ov <- pairs[pairs$i==i[1] & pairs$j==i[2],]
    data.frame(fename=name,id_a=d$BOUND_A_ID[i[1]],id_b=d$BOUND_A_ID[i[2]],
      polygon_gap_km=as.numeric(st_distance(g[i[1],],g[i[2],]))/1000,
      centroid_separation_km=as.numeric(st_distance(st_centroid(st_geometry(g[i[1],])),
                                                   st_centroid(st_geometry(g[i[2],]))))/1000,
      overlap_km2=sum(ov$area_km2),same_geometry=nrow(ov)>0L && any(ov$equals))
  }))
  d$same_name_overlap_km2 <- both$overlap_km2[match(d$FENAME,both$fename)]
  d$same_name_polygon_gap_km <- both$polygon_gap_km[match(d$FENAME,both$fename)]
  d$same_name_centroid_separation_km <- both$centroid_separation_km[match(d$FENAME,both$fename)]
  list(audit=d,all=g,refs=refs,box=b)
}

build_maritime <- function(x) {
  m <- x$maritime
  stopifnot(all(st_geometry_type(m)=="LINESTRING"),all(st_is_valid(m)),!anyDuplicated(m$BOUND_L_ID))
  m$segment_id <- sprintf("SU%02d",seq_len(nrow(m)))
  east <- which(m$BOUND_L_ID==9827)
  stopifnot(identical(east,4L),st_bbox(m[east,])["xmin"]>122.5)
  list(east=m[east,],south=m[-east,],all=m)
}

prepare <- function(input) {
  sf_use_s2(FALSE)
  x <- inputs(input)
  x$province <- build_province_lines(x)
  x$south <- build_south(x)
  x$east <- build_east(x)
  x$mar <- build_maritime(x)
  # Same source context for both projections. Crop an intentionally larger window in WGS84.
  b <- box(c(100,-3,128,29),wgs)
  context <- function(y) suppressWarnings(st_crop(y[window_intersects(st_geometry(y),b),],b))
  x$context_land <- context(x$land)
  x$context_coast <- context(x$coast)
  # Use the source's nearshore polygon components and keep its offshore pieces out of the main panel.
  # The complete China POV geometry is retained above for clipping the internal province lines.
  c <- st_cast(st_geometry(x$country),"POLYGON")
  main_or_coastal <- coords_wgs(st_point_on_surface(c))[,2] >= 18
  x$main_china <- x$inset_china <- c[main_or_coastal]
  x
}


make_layers <- function(x) {
  province_rows <- x$province$audit$source_row[x$province$audit$drawn_m>0]
  stopifnot(length(province_rows)==length(x$province$geometry))
  take <- function(y, cols) y[,cols,drop=FALSE]
  mar <- x$mar$all
  mar$panel <- ifelse(mar$segment_id=="SU04","main","inset")
  # Display-only substitution of the three tiny NE polygon components in the
  # official reference-point cluster. Preserve all original components in the
  # source layer; do not change vertices or infer individual island identities.
  source <- st_as_sf(x$main_china)
  source$source_component <- seq_len(nrow(source))
  bb <- st_transform(st_as_sfc(box(c(123.2,25.55,124.85,26.12),wgs)),aea)
  replace <- which(lengths(st_intersects(st_geometry(source),bb))>0L)
  stopifnot(identical(replace,72:74), all(as.numeric(st_area(source[replace,]))/1e6<4))
  display <- source[-replace,]
  refs <- x$east$refs
  refs$symbol_mm <- c(.45,.45,.45,.06,.06)
  refs$display_rule <- c(rep("Official location point replaces cluster-scale source outline",3),
    rep("Small point at unchanged reference coordinate; no displacement",2))
  layers <- list(
    china_pov_nearshore=source,
    china_pov_display=display,
    china_sar=take(x$sar,"GU_A3"),
    province_internal_lines=st_sf(source_row=province_rows,geometry=x$province$geometry),
    inset_land=take(x$context_land,character()),
    inset_coast=take(x$context_coast,character()),
    maritime_lines=take(mar,c("BOUND_L_ID","segment_id","panel")),
    east_island_points=take(refs,c("fename","label","reference_url","symbol_mm","display_rule")),
    south_reef_outlines=take(x$south$reef,c("BOUND_A_ID","FENAME")),
    south_island_polygons=take(x$south$land,c("BOUND_A_ID","FENAME")),
    south_feature_points=take(x$south$points,c("id","symbol"))
  )
  stopifnot(identical(vapply(layers[c("province_internal_lines","maritime_lines","east_island_points",
    "south_reef_outlines","south_island_polygons","south_feature_points")],nrow,integer(1)),
    c(province_internal_lines=75L, maritime_lines=10L, east_island_points=5L,
      south_reef_outlines=175L, south_island_polygons=7L, south_feature_points=6L)))
  layers
}

build_geopackage <- function(raw_dir, output) {
  if (file.exists(output)) stop("Output already exists: ", output)
  if (!dir.exists(dirname(output))) stop("Output directory does not exist: ", dirname(output))
  check_sources(raw_dir)
  layers <- make_layers(prepare(raw_dir))
  temp <- tempfile("s1_geography_",fileext=".gpkg")
  on.exit({ invisible(gc()); unlink(temp) }, add=TRUE)
  for (name in names(layers)) {
    layer <- layers[[name]]
    stopifnot(nrow(layer)>0L, !is.na(st_crs(layer)), !any(st_is_empty(layer)))
    st_write(layer,temp,layer=name,quiet=TRUE,delete_layer=FALSE)
    message(name, ": ", nrow(layer), " features; ", paste(unique(st_geometry_type(layer)),collapse="/"))
  }
  stopifnot(identical(st_layers(temp)$name,names(layers)))
  if (!file.copy(temp,output,overwrite=FALSE)) stop("Could not save ",output)
  message("Rebuilt GeoPackage: ", output)
}

if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if (!length(args) || length(args)>2L)
    stop("Usage: Rscript --vanilla prepare_geography.R RAW_DIR [OUTPUT_GPKG]")
  print(sessionInfo())
  print(sf_extSoftVersion())
  build_geopackage(args[1], if (length(args)==2L) args[2] else "S1_geography.gpkg")
}
