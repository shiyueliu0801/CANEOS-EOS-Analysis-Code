# CANEOS Figure 3 v3 | Minimal Public Reproducible Script
# Rebuilds the locked Figure 3 v3 Sankey diagram (7.50 x 4.00 in) strictly from
# the aggregated three-phase trajectory counts and palette CSVs in ./data.
# Does not require or access infant-level raw records or external workspaces.
# Usage: Rscript --vanilla Figure3_Reproduce.R [output_dir] [data_dir]

if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(scales)
})

check <- function(condition, reason) {
  if (!isTRUE(condition)) stop(reason, call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")) else getwd()
out_dir  <- if (exists("fig3_out_dir", inherits = TRUE)) fig3_out_dir else if (length(args) >= 1L) args[1] else script_dir
data_dir <- if (exists("fig3_data_dir", inherits = TRUE)) fig3_data_dir else if (length(args) >= 2L) args[2] else file.path(script_dir, "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
check(dir.exists(out_dir), paste0("Cannot create output directory: ", out_dir))
check(dir.exists(data_dir), paste0("Data directory does not exist: ", data_dir))
out <- function(x) file.path(out_dir, x)

# 1. Dependency, graphics device, and Times New Roman font checks
font_family <- "Times New Roman"
check(requireNamespace("ggtext", quietly = TRUE), "The ggtext package is required for italic n labels.")
check(requireNamespace("systemfonts", quietly = TRUE), "The systemfonts package is required for strict font checks.")
check(requireNamespace("ragg", quietly = TRUE), "The ragg package is required for validated TIFF/PNG raster export.")
check(capabilities("cairo"), "Cairo PDF capability is required for embedded-font vector PDF export.")

for (i in seq_len(4)) {
  f <- suppressWarnings(systemfonts::match_font(
    font_family,
    italic = i %in% c(3L, 4L),
    bold = i %in% c(2L, 4L)
  ))
  info <- systemfonts::font_info(path = f$path, index = f$index)
  check(all(info$family == font_family), "Times New Roman was substituted by another font.")
}

# 2. Load and validate aggregated inputs
traj_path <- file.path(data_dir, "Figure3_Three_Phase_Trajectories.csv")
check(file.exists(traj_path), paste0("Missing required trajectory file: ", traj_path))
pal_df <- if (exists("static_cfg", inherits = TRUE)) {
  static_cfg$get_figure3_palette()
} else {
  source(file.path(dirname(dirname(dirname(script_dir))), "config", "field_mappings.R"), local = TRUE, encoding = "UTF-8")
  get_figure3_palette()
}

terminal <- "Death or treatment withdrawal ≤72 h"
drug_levels <- c("CARB+GLY/LZD", "GLY/LZD", "CARB", "3GC+PEN", "3GC",
                 "1GC/2GC", "BSP", "PEN", "Other/UNK")
lvl <- c(drug_levels, terminal)

check(nrow(pal_df) == length(lvl) && identical(pal_df$category, lvl),
      "Palette categories or order do not match the locked Figure 3 v3 specification.")
check(all(grepl("^#[0-9A-Fa-f]{6}$", pal_df$hex_colour)),
      "Invalid hex colour format in Figure3_Palette.csv.")
pal <- setNames(pal_df$hex_colour, pal_df$category)

trajectories <- suppressMessages(read_csv(
  traj_path,
  col_types = cols(T1 = col_character(), T2 = col_character(), T3 = col_character(), n = col_integer()),
  locale = locale(encoding = "UTF-8"),
  show_col_types = FALSE, trim_ws = FALSE, na = c("")
))
check(nrow(problems(trajectories)) == 0L, "CSV parsing error in Figure3_Three_Phase_Trajectories.csv.")
check(all(c("T1", "T2", "T3", "n") %in% names(trajectories)),
      "Missing required columns (T1, T2, T3, n) in Figure3_Three_Phase_Trajectories.csv.")
check(!anyNA(trajectories$T1) && !anyNA(trajectories$T2) && !anyNA(trajectories$T3) && !anyNA(trajectories$n),
      "Missing values detected in Figure3_Three_Phase_Trajectories.csv.")
check(all(trajectories$n > 0L), "All aggregated cell counts n must be positive integers.")
check(all(trajectories$T1 %in% lvl) && all(trajectories$T2 %in% lvl) && all(trajectories$T3 %in% lvl),
      "Unexpected regimen/status category in Figure3_Three_Phase_Trajectories.csv.")
check(sum(trajectories$n) == 1189L,
      sprintf("Total infant count in trajectory table is %d; expected 1189.", sum(trajectories$n)))

# Reconstruct complete three-phase table by expanding aggregated counts
idx <- rep(seq_len(nrow(trajectories)), trajectories$n)
plot_df <- tibble(
  T1 = factor(trajectories$T1[idx], levels = lvl),
  T2 = factor(trajectories$T2[idx], levels = lvl),
  T3 = factor(trajectories$T3[idx], levels = lvl)
)

# Stage, adjacent-flow, terminal-count, and balance conservation checks
flow12 <- trajectories |> group_by(T1, T2) |> summarise(n = sum(n), .groups = "drop") |> rename(from = T1, to = T2)
flow23 <- trajectories |> group_by(T2, T3) |> summarise(n = sum(n), .groups = "drop") |> rename(from = T2, to = T3)
node_counts <- bind_rows(lapply(c("T1", "T2", "T3"), function(v) {
  tibble(phase = v, category = lvl, n = as.integer(table(factor(plot_df[[v]], levels = lvl))))
}))
check(all(tapply(node_counts$n, node_counts$phase, sum) == 1189L) &&
        sum(flow12$n) == 1189L && sum(flow23$n) == 1189L,
      "Stage or adjacent-flow totals do not sum to 1189.")
check(node_counts$n[node_counts$phase == "T1" & node_counts$category == terminal] == 0L &&
        node_counts$n[node_counts$phase == "T2" & node_counts$category == terminal] == 0L &&
        node_counts$n[node_counts$phase == "T3" & node_counts$category == terminal] == 77L,
      "Unexpected terminal early-exit count across phases (expected 0, 0, 77).")

for (cat_i in lvl) {
  t1_n <- node_counts$n[node_counts$phase == "T1" & node_counts$category == cat_i]
  t2_n <- node_counts$n[node_counts$phase == "T2" & node_counts$category == cat_i]
  t3_n <- node_counts$n[node_counts$phase == "T3" & node_counts$category == cat_i]
  check(sum(flow12$n[flow12$from == cat_i]) == t1_n, paste0("T1 outflow mismatch for ", cat_i))
  check(sum(flow12$n[flow12$to == cat_i]) == t2_n && sum(flow23$n[flow23$from == cat_i]) == t2_n,
        paste0("T2 inflow/outflow mismatch for ", cat_i))
  check(sum(flow23$n[flow23$to == cat_i]) == t3_n, paste0("T3 inflow mismatch for ", cat_i))
}
check(all(abs(tapply(node_counts$n / 1189 * 100, node_counts$phase, sum) - 100) < 1e-10),
      "Unrounded phase percentages do not sum to 100%.")

# 3. Geometry helpers (unchanged from locked Figure 3 v3 pipeline)
segment_colour <- "destination"
node_width <- 0.044
axis_gap <- 6
curve_n <- 80
ribbon_alpha <- 0.42
ribbon_alpha_thin <- 0.58

make_axis_nodes <- function(x, axis_name, x_pos, levels_all = lvl, gap = axis_gap) {
  x_chr <- as.character(x)
  tb <- table(factor(x_chr, levels = levels_all), useNA = "no")
  nodes <- tibble(cat = names(tb), n = as.numeric(tb)) |>
    filter(n > 0) |>
    mutate(order_id = match(cat, levels_all)) |>
    arrange(order_id)

  total_height <- sum(nodes$n) + gap * (nrow(nodes) - 1)
  cursor_top <- total_height

  out <- vector("list", nrow(nodes))
  for (i in seq_len(nrow(nodes))) {
    n_i <- nodes$n[i]
    ymax_i <- cursor_top
    ymin_i <- ymax_i - n_i
    out[[i]] <- tibble(
      axis = axis_name,
      x = x_pos,
      cat = nodes$cat[i],
      n = n_i,
      ymin = ymin_i,
      ymax = ymax_i,
      ymid = (ymin_i + ymax_i) / 2,
      order_id = nodes$order_id[i]
    )
    cursor_top <- ymin_i - gap
  }
  bind_rows(out)
}

make_flow_positions <- function(df_pair, from_var, to_var, from_nodes, to_nodes, levels_all = lvl) {
  flows <- df_pair |>
    count(.data[[from_var]], .data[[to_var]], name = "n") |>
    rename(from = 1, to = 2) |>
    filter(!is.na(from), !is.na(to), n > 0) |>
    mutate(
      from = as.character(from),
      to   = as.character(to),
      from_order = match(from, levels_all),
      to_order   = match(to, levels_all)
    )

  source_out <- vector("list", 0)
  for (src in levels_all) {
    sub <- flows |>
      filter(from == src) |>
      arrange(to_order, from_order, desc(n))
    if (nrow(sub) == 0) next

    ymax0 <- from_nodes$ymax[from_nodes$cat == src]
    cursor_top <- ymax0
    tmp <- vector("list", nrow(sub))
    for (i in seq_len(nrow(sub))) {
      n_i <- sub$n[i]
      ymax_i <- cursor_top
      ymin_i <- ymax_i - n_i
      tmp[[i]] <- sub[i, ] |>
        mutate(ymin_from = ymin_i, ymax_from = ymax_i)
      cursor_top <- ymin_i
    }
    source_out <- c(source_out, tmp)
  }
  flows2 <- bind_rows(source_out)

  target_out <- vector("list", 0)
  for (dst in levels_all) {
    sub <- flows2 |>
      filter(to == dst) |>
      arrange(from_order, to_order, desc(n))
    if (nrow(sub) == 0) next

    ymax0 <- to_nodes$ymax[to_nodes$cat == dst]
    cursor_top <- ymax0
    tmp <- vector("list", nrow(sub))
    for (i in seq_len(nrow(sub))) {
      n_i <- sub$n[i]
      ymax_i <- cursor_top
      ymin_i <- ymax_i - n_i
      tmp[[i]] <- sub[i, ] |>
        mutate(ymin_to = ymin_i, ymax_to = ymax_i)
      cursor_top <- ymin_i
    }
    target_out <- c(target_out, tmp)
  }
  bind_rows(target_out)
}

make_one_ribbon <- function(x1, x2, ymin1, ymax1, ymin2, ymax2, id, fill_cat, alpha_val, n_points = curve_n) {
  t <- seq(0, 1, length.out = n_points)
  s <- 3 * t^2 - 2 * t^3
  x <- x1 + (x2 - x1) * t
  y_upper <- ymax1 + (ymax2 - ymax1) * s
  y_lower <- ymin1 + (ymin2 - ymin1) * s

  tibble(
    group = id,
    x = c(x, rev(x)),
    y = c(y_upper, rev(y_lower)),
    fill_cat = fill_cat,
    alpha_val = alpha_val
  )
}

make_ribbons <- function(flow_pos, x1, x2, colour_by = c("destination", "source"), thin_cutoff = 8) {
  colour_by <- match.arg(colour_by)
  out <- vector("list", nrow(flow_pos))

  for (i in seq_len(nrow(flow_pos))) {
    fill_cat_i <- if (colour_by == "destination") flow_pos$to[i] else flow_pos$from[i]
    alpha_i <- ifelse(flow_pos$n[i] <= thin_cutoff, ribbon_alpha_thin, ribbon_alpha)
    out[[i]] <- make_one_ribbon(
      x1 = x1, x2 = x2,
      ymin1 = flow_pos$ymin_from[i], ymax1 = flow_pos$ymax_from[i],
      ymin2 = flow_pos$ymin_to[i],   ymax2 = flow_pos$ymax_to[i],
      id = i,
      fill_cat = fill_cat_i,
      alpha_val = alpha_i,
      n_points = curve_n
    )
  }
  bind_rows(out)
}

# 4. Build nodes, ribbons, and v3 presentation layout
display_terminal <- "No subsequent regimen†"
figure_width_in <- 7.5
figure_height_in <- 4.0
pt_mm <- function(pt) pt / (72.27 / 25.4)

node1 <- make_axis_nodes(plot_df$T1, "Initial empirical phase", 1)
node2 <- make_axis_nodes(plot_df$T2, "Early reassessment phase (0–72 h)", 2)
node3 <- make_axis_nodes(plot_df$T3, "Late-treatment status (72 h–7 d)", 3)
node_df <- bind_rows(node1, node2, node3)

flow12_pos <- make_flow_positions(plot_df, "T1", "T2", node1, node2, lvl)
flow23_pos <- make_flow_positions(plot_df, "T2", "T3", node2, node3, lvl)
ribbon12 <- make_ribbons(flow12_pos, 1 + node_width / 2, 2 - node_width / 2, segment_colour) |>
  mutate(segment = "T1_to_T2")
ribbon23 <- make_ribbons(flow23_pos, 2 + node_width / 2, 3 - node_width / 2, segment_colour) |>
  mutate(segment = "T2_to_T3", group = group + max(ribbon12$group))
ribbon_df <- bind_rows(ribbon12, ribbon23)

check(sum(flow12_pos$n) == 1189L && sum(flow23_pos$n) == 1189L,
      "Ribbon flow totals are not conserved.")
check(all(flow12_pos$ymax_from - flow12_pos$ymin_from == flow12_pos$n) &&
        all(flow23_pos$ymax_to - flow23_pos$ymin_to == flow23_pos$n),
      "Ribbon widths do not equal infant counts.")

xlim_range <- c(0.69, 3.53)
ylim_lo <- -20
ylim_hi <- 1385
panel_h_pt <- 228
y_per_pt <- (ylim_hi - ylim_lo) / panel_h_pt

label_df <- node_df |>
  filter(x %in% c(1, 3)) |>
  mutate(
    is_other = (cat == "Other/UNK"),
    display_name = if_else(cat == terminal, display_terminal, cat),
    percent = n / 1189 * 100,
    statistic = if_else(is_other, "", sprintf("%.1f%%; <i>n</i> = %d", percent, n)),
    half_h_pt = if_else(is_other, 3.8, 7.9),
    label_x = if_else(x == 1, x - node_width / 2 - 0.025, x + node_width / 2 + 0.025),
    label_hjust = if_else(x == 1, 1, 0),
    label_center = ymid
  )

min_gap_pt <- 2.1
for (axis_i in c(1, 3)) {
  ii <- which(label_df$x == axis_i)
  ii <- ii[order(label_df$ymid[ii])]
  req_dists <- (label_df$half_h_pt[ii[-length(ii)]] + label_df$half_h_pt[ii[-1]] + min_gap_pt) * y_per_pt
  offsets <- c(0, cumsum(req_dists))
  pos <- stats::isoreg(seq_along(ii), label_df$ymid[ii] - offsets)$yf + offsets
  top_cap <- max(node_df$ymax[node_df$x == axis_i]) - label_df$half_h_pt[ii[length(ii)]] * y_per_pt + 2 * y_per_pt
  if (pos[length(pos)] > top_cap) {
    pos[length(pos)] <- top_cap
    for (k in rev(seq_len(length(pos) - 1))) {
      pos[k] <- min(pos[k], pos[k + 1] - req_dists[k])
    }
  }
  label_df$label_center[ii] <- pos
}

line_off_y <- 4.05 * y_per_pt
label_df <- label_df |>
  mutate(
    name_y = if_else(is_other, label_center, label_center + line_off_y),
    statistic_y = label_center - line_off_y,
    label_shift = label_center - ymid,
    node_edge_x = if_else(x == 1, x - node_width / 2, x + node_width / 2)
  )
leaders <- label_df |> filter(abs(label_shift) > 7.0 * y_per_pt)

title_df <- tibble(
  x = c(1, 2, 3),
  y = 1352,
  title = c("Initial empirical phase", "Early reassessment phase", "Later treatment status"),
  subtitle = c("Start of therapy", "(0–72 h after birth)", "(72 h–7 d after birth)"),
  subtitle_y = 1290
)

p <- ggplot() +
  geom_polygon(data = ribbon_df, aes(x = x, y = y, group = interaction(segment, group), fill = fill_cat, alpha = alpha_val), colour = NA) +
  geom_rect(data = node_df, aes(xmin = x - node_width / 2, xmax = x + node_width / 2, ymin = ymin, ymax = ymax, fill = cat),
            colour = "white", linewidth = 0.45) +
  geom_segment(data = leaders, aes(x = node_edge_x, y = ymid, xend = label_x, yend = label_center),
               linewidth = 0.16, colour = "#B3B3B3", inherit.aes = FALSE) +
  geom_text(data = label_df, aes(x = label_x, y = name_y, label = display_name, hjust = label_hjust),
            family = font_family, size = pt_mm(7.5), fontface = "bold", colour = "#404040") +
  ggtext::geom_richtext(data = label_df |> filter(!is_other), aes(x = label_x, y = statistic_y, label = statistic, hjust = label_hjust),
                        family = font_family, size = pt_mm(7.0), fontface = "plain", vjust = 0.5,
                        fill = NA, label.color = NA, label.padding = grid::unit(rep(0, 4), "pt"),
                        label.margin = grid::unit(rep(0, 4), "pt"), colour = "#606060") +
  geom_text(data = title_df, aes(x = x, y = y, label = title), family = font_family,
            size = pt_mm(8.0), fontface = "bold", colour = "#000000") +
  geom_text(data = title_df, aes(x = x, y = subtitle_y, label = subtitle), family = font_family,
            size = pt_mm(7.0), fontface = "plain", colour = "#606060") +
  scale_fill_manual(values = pal, limits = drug_levels, drop = FALSE, name = NULL) +
  scale_alpha_identity() +
  coord_cartesian(xlim = xlim_range, ylim = c(ylim_lo, ylim_hi), expand = FALSE, clip = "off") +
  theme_void(base_size = 7.5, base_family = font_family) +
  theme(text = element_text(family = font_family),
        legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
        legend.text = element_text(family = font_family, face = "plain", size = 7.5, colour = "#606060"),
        legend.key.width = grid::unit(7.0, "pt"), legend.key.height = grid::unit(7.0, "pt"),
        legend.spacing.x = grid::unit(5.5, "pt"),
        legend.margin = margin(4, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(5, 6, 11, 6, unit = "pt")) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE), alpha = "none")

# 5. Export final Figure 3 v3 TIFF (300 dpi, 2250 x 1200 px = 7.50 x 4.00 in / 19.05 x 10.16 cm),
# vector PDF (540 x 288 pt), and 300 dpi PNG preview
ggsave(out("Figure3_Sankey_Final_v3_300dpi.tiff"), plot = p, width = figure_width_in, height = figure_height_in,
       units = "in", dpi = 300, compression = "lzw", bg = "white", device = "tiff")
ggsave(out("Figure3_Sankey_Final_v3.pdf"), plot = p, width = figure_width_in, height = figure_height_in,
       units = "in", bg = "white", device = cairo_pdf, family = font_family)
ggsave(out("Figure3_Sankey_Final_v3_300dpi.png"), plot = p, width = figure_width_in, height = figure_height_in,
       units = "in", dpi = 300, bg = "white", device = "png")

message(sprintf("Reproduced Figure 3 v3 in %s: %d infants; terminal n=%d (%.1f%%).",
                out_dir, nrow(plot_df),
                node_counts$n[node_counts$phase == "T3" & node_counts$category == terminal],
                100 * node_counts$n[node_counts$phase == "T3" & node_counts$category == terminal] / nrow(plot_df)))

