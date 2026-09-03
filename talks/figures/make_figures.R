## make_figures.R
## Generates the conceptual figures used in the talk slides.
## Run from /Users/macleand/Desktop/peppwR/talks/figures/ :
##   Rscript make_figures.R
## Outputs PNGs in the current directory, named per the slide they sit on.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(pwr)
})

set.seed(42)

# House style ------------------------------------------------------------

accent <- "#C8102E"   # chilli red, matches the peppwR hex
mid    <- "grey45"
soft   <- "grey80"

theme_talk <- function() {
  theme_minimal(base_size = 16) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      axis.title       = element_text(color = mid),
      axis.text        = element_text(color = mid),
      plot.title       = element_text(color = mid, size = 16, face = "plain"),
      plot.subtitle    = element_text(color = mid, size = 13)
    )
}

theme_anatomy <- function() {
  theme_void(base_size = 16) +
    theme(
      plot.margin = margin(20, 20, 20, 20),
      legend.position = "none"
    )
}

save_png <- function(plot, file, width = 10, height = 5.5, dpi = 150) {
  ggsave(file, plot, width = width, height = height, dpi = dpi, bg = "white")
  message("wrote ", file)
}

# Slide 6 — four levers, anatomically ------------------------------------

{
  dx <- seq(-3, 8, length.out = 600)
  ctrl <- dnorm(dx, 0, 1)
  trt  <- dnorm(dx, 4, 1)
  df <- data.frame(x = c(dx, dx),
                   y = c(ctrl, trt),
                   g = rep(c("control", "treatment"), each = length(dx)))

  crit <- qnorm(0.95)
  tail <- subset(df, g == "control" & x >= crit)

  peak_ctrl <- dnorm(0, 0, 1)
  peak_trt  <- dnorm(4, 4, 1)

  p_four <- ggplot(df, aes(x, y, colour = g)) +
    geom_area(data = tail, aes(x, y), inherit.aes = FALSE,
              fill = mid, alpha = 0.35) +
    geom_line(linewidth = 1.4) +
    # vertical centre lines grounding both peaks
    annotate("segment", x = 0, xend = 0, y = 0, yend = peak_ctrl,
             colour = mid, linetype = "dashed", linewidth = 0.4) +
    annotate("segment", x = 4, xend = 4, y = 0, yend = peak_trt,
             colour = accent, linetype = "dashed", linewidth = 0.4) +
    scale_colour_manual(values = c(control = mid, treatment = accent)) +
    # delta arrow between the two peaks
    annotate("segment", x = 0, xend = 4, y = 0.47, yend = 0.47,
             arrow = arrow(ends = "both", length = unit(0.2, "in")),
             colour = "black", linewidth = 0.9) +
    annotate("text", x = 2, y = 0.52, label = "delta",
             parse = TRUE, size = 13) +
    # sigma half-width on treatment curve
    annotate("segment", x = 4, xend = 5, y = 0.22, yend = 0.22,
             arrow = arrow(ends = "both", length = unit(0.16, "in")),
             colour = "black", linewidth = 0.8) +
    annotate("text", x = 4.5, y = 0.28, label = "sigma",
             parse = TRUE, size = 11) +
    # alpha tail label
    annotate("text", x = 2.5, y = 0.06, label = "alpha",
             parse = TRUE, size = 11) +
    # group labels (large)
    annotate("text", x = -1.8, y = 0.42, label = "control", colour = mid, size = 8) +
    annotate("text", x = 5.8, y = 0.42, label = "treatment", colour = accent, size = 8) +
    coord_cartesian(ylim = c(0, 0.58)) +
    theme_anatomy()

  save_png(p_four, "fig_four_levers.png", width = 10, height = 5.5)
}

# Slide 7 — overlap low / high power -------------------------------------

{
  dx <- seq(-4, 8, length.out = 600)
  low <- data.frame(
    x = c(dx, dx),
    y = c(dnorm(dx, 0, 1.5), dnorm(dx, 1.5, 1.5)),
    g = rep(c("control", "treatment"), each = length(dx)),
    panel = "low power"
  )
  high <- data.frame(
    x = c(dx, dx),
    y = c(dnorm(dx, 0, 1), dnorm(dx, 4, 1)),
    g = rep(c("control", "treatment"), each = length(dx)),
    panel = "high power"
  )
  df <- bind_rows(low, high)
  df$panel <- factor(df$panel, levels = c("low power", "high power"))

  p <- ggplot(df, aes(x, y, colour = g)) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = c(control = mid, treatment = accent)) +
    facet_wrap(~ panel, nrow = 1) +
    coord_cartesian(ylim = c(0, 0.45)) +
    theme_anatomy() +
    theme(
      strip.text = element_text(colour = mid, size = 16, face = "plain",
                                margin = margin(b = 10))
    )

  save_png(p, "fig_overlap_low_high.png", width = 11, height = 4.5)
}

# Slide 9 — what p actually is -------------------------------------------

{
  dx <- seq(-4, 4, length.out = 600)
  df <- data.frame(x = dx, y = dnorm(dx))
  crit <- 1.9
  tail <- subset(df, x >= crit)

  p <- ggplot(df, aes(x, y)) +
    geom_area(data = tail, fill = mid, alpha = 0.45) +
    geom_line(colour = mid, linewidth = 1.1) +
    geom_vline(xintercept = crit, colour = accent, linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = crit + 0.05, y = 0.31, label = "observed treatment mean",
             colour = accent, hjust = 0, size = 4.5) +
    annotate("text", x = 0, y = 0.43, label = "null distribution of the control mean",
             colour = mid, size = 4.5) +
    annotate("text", x = 2.7, y = 0.04, label = "this area = p",
             colour = "black", size = 4.5) +
    coord_cartesian(ylim = c(0, 0.47)) +
    theme_anatomy()

  save_png(p, "fig_pvalue_tail.png", width = 10, height = 5)
}

# Slide 12 — classical power curve ---------------------------------------

{
  ns <- 2:100
  build <- function(d) {
    sapply(ns, function(n) pwr.t.test(n = n, d = d, sig.level = 0.05,
                                      type = "two.sample",
                                      alternative = "two.sided")$power)
  }
  df <- bind_rows(
    data.frame(n = ns, power = build(0.2), D = "0.2 (small)"),
    data.frame(n = ns, power = build(0.5), D = "0.5 (medium)"),
    data.frame(n = ns, power = build(0.8), D = "0.8 (large)")
  )

  cross <- data.frame(
    D = c("0.2 (small)", "0.5 (medium)", "0.8 (large)"),
    n = c(393, 64, 26),
    power = 0.8
  )

  p <- ggplot(df, aes(n, power, colour = D)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = mid) +
    geom_vline(xintercept = 3, linetype = "dotted", colour = accent) +
    geom_line(linewidth = 1.1) +
    annotate("text", x = 3, y = 0.05, label = "n = 3",
             colour = accent, hjust = -0.15, size = 4) +
    scale_colour_manual(values = c("0.2 (small)" = soft,
                                   "0.5 (medium)" = mid,
                                   "0.8 (large)" = accent)) +
    scale_x_continuous(breaks = c(3, 26, 50, 64, 100), limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = "n per group", y = "power",
         colour = "Cohen's D") +
    theme_talk() +
    theme(legend.position = c(0.78, 0.25),
          legend.background = element_rect(fill = "white", colour = NA))

  save_png(p, "fig_classical_power_curve.png", width = 10, height = 5.5)
}

# Slide 14 — tech vs bio reps schematic ----------------------------------

{
  # left panel: 1 subject -> 4 measurements
  left_subjects <- data.frame(x = 0.5, y = 0.75)
  left_tubes <- data.frame(x = seq(0.1, 0.9, length.out = 4), y = 0.25)
  left_edges <- data.frame(
    x = 0.5, xend = left_tubes$x,
    y = 0.7,  yend = 0.3
  )

  p_left <- ggplot() +
    geom_segment(data = left_edges, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = soft, linewidth = 0.8) +
    geom_point(data = left_subjects, aes(x, y),
               size = 30, shape = 21, fill = accent, colour = "white", stroke = 1.5) +
    geom_point(data = left_tubes, aes(x, y),
               size = 14, shape = 22, fill = mid, colour = "white", stroke = 1) +
    annotate("text", x = 0.5, y = 0.96,
             label = "1 biological replicate",
             colour = mid, size = 7) +
    annotate("text", x = 0.5, y = 0.08,
             label = "4 measurements", colour = mid, size = 6) +
    annotate("text", x = 0.5, y = -0.05,
             label = expression(n[bio] == 1 ~ "—" ~ n == 1 ~ ", not 4"),
             colour = "black", size = 7) +
    coord_cartesian(xlim = c(0, 1), ylim = c(-0.15, 1.05)) +
    theme_anatomy()

  # right panel: 4 subjects -> 4 measurements
  right_subjects <- data.frame(x = seq(0.1, 0.9, length.out = 4), y = 0.75)
  right_tubes    <- data.frame(x = seq(0.1, 0.9, length.out = 4), y = 0.25)
  right_edges    <- data.frame(
    x = right_subjects$x, xend = right_tubes$x,
    y = 0.7, yend = 0.3
  )

  p_right <- ggplot() +
    geom_segment(data = right_edges, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = soft, linewidth = 0.8) +
    geom_point(data = right_subjects, aes(x, y),
               size = 24, shape = 21, fill = accent, colour = "white", stroke = 1.5) +
    geom_point(data = right_tubes, aes(x, y),
               size = 14, shape = 22, fill = mid, colour = "white", stroke = 1) +
    annotate("text", x = 0.5, y = 0.96,
             label = "4 biological replicates",
             colour = mid, size = 7) +
    annotate("text", x = 0.5, y = 0.08,
             label = "4 measurements", colour = mid, size = 6) +
    annotate("text", x = 0.5, y = -0.05,
             label = expression(n[bio] == 4 ~ "—" ~ n == 4),
             colour = "black", size = 7) +
    coord_cartesian(xlim = c(0, 1), ylim = c(-0.15, 1.05)) +
    theme_anatomy()

  p_schematic <- p_left + p_right + plot_layout(ncol = 2)
  save_png(p_schematic, "fig_tech_vs_bio_reps.png", width = 11, height = 5)
}

# Slide 15 — ordinal merge schematic -------------------------------------

{
  raw <- data.frame(score = c(0, 2, 3, 1, 4, 2))
  ranked <- raw |>
    mutate(rank = rank(score))
  by_rep <- data.frame(rep = c("rep 1", "rep 2"),
                       mean_rank = c(mean(ranked$rank[1:3]),
                                     mean(ranked$rank[4:6])))

  # We render this in Quarto as hand-set tables; PNG fallback here for safety.
  tbl_theme <- theme_void(base_size = 18) +
    theme(plot.margin = margin(20, 10, 20, 10))

  mk_table <- function(df, title) {
    df_long <- df |>
      mutate(row = row_number())
    ggplot(df_long) +
      geom_text(aes(x = 1, y = -row,
                    label = apply(df, 1, paste, collapse = "  ")),
                size = 6, family = "mono") +
      annotate("text", x = 1, y = 0.5, label = title, fontface = "bold", size = 5) +
      coord_cartesian(xlim = c(0, 2)) +
      tbl_theme
  }

  # The R fallback is ugly. The qmd should render this as native tables.
  # Just emit a placeholder PNG noting that.
  message("Note: slide 15 ordinal-merge schematic is rendered as native Quarto tables in slides.qmd; no PNG produced.")
}

# Slide 18 — HR re-analysis stylised mock --------------------------------

{
  set.seed(7)
  # left panel: ordinal-as-numeric violin, tech reps pooled
  pool <- data.frame(
    condition = rep(c("A", "B"), each = 60),
    score = c(
      sample(0:6, 60, replace = TRUE,
             prob = c(0.05, 0.1, 0.2, 0.3, 0.2, 0.1, 0.05)),
      sample(0:6, 60, replace = TRUE,
             prob = c(0.02, 0.05, 0.1, 0.2, 0.3, 0.2, 0.13))
    )
  )

  p_left <- ggplot(pool, aes(condition, score, fill = condition)) +
    geom_violin(colour = NA, alpha = 0.85, width = 0.7) +
    geom_jitter(width = 0.18, height = 0, size = 0.6,
                colour = "grey25", alpha = 0.5) +
    scale_fill_manual(values = c(A = mid, B = accent)) +
    scale_y_continuous(breaks = 0:6, limits = c(-0.5, 6.5)) +
    labs(title = "as published",
         subtitle = "tech reps pooled into n; ordinal scores treated as numeric",
         x = NULL, y = "HR score") +
    theme_talk() +
    theme(legend.position = "none")

  # right panel: count plot, biological n only (collapse to two reps)
  by_rep <- pool |>
    mutate(rep = rep(c("R1", "R2"), times = 60)) |>
    count(condition, rep, score)

  p_right <- ggplot(by_rep, aes(condition, score, size = n, colour = condition)) +
    geom_point(alpha = 0.85, position = position_jitter(width = 0.15, height = 0, seed = 1)) +
    scale_colour_manual(values = c(A = mid, B = accent)) +
    scale_size_area(max_size = 12, guide = "none") +
    scale_y_continuous(breaks = 0:6, limits = c(-0.5, 6.5)) +
    labs(title = "re-analysed",
         subtitle = "biological n; rank-and-average within rep",
         x = NULL, y = NULL) +
    theme_talk() +
    theme(legend.position = "none")

  # banner annotations: 75/120 -> 42/120
  banner_left  <- annotate("text", x = 0, y = 0, label = "")
  hr <- p_left + p_right +
    plot_annotation(
      title = "75 / 120 significant   →   42 / 120 significant",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 17, colour = mid))
    )

  save_png(hr, "fig_hr_reanalysis_mock.png", width = 12, height = 5.5)
}

# Slide 19a — many overlapping curves ------------------------------------

{
  dx <- seq(-3, 8, length.out = 400)
  set.seed(11)
  curves <- lapply(1:25, function(i) {
    mu <- runif(1, -1, 5)
    sd <- runif(1, 0.6, 1.6)
    data.frame(x = dx, y = dnorm(dx, mu, sd), id = i)
  }) |> bind_rows()

  p <- ggplot(curves, aes(x, y, group = id)) +
    geom_line(colour = mid, alpha = 0.25, linewidth = 0.6) +
    annotate("text", x = 2.5, y = 0.7,
             label = "thousands of simultaneous comparisons",
             colour = "black", size = 5) +
    coord_cartesian(ylim = c(0, 0.75)) +
    theme_anatomy()

  save_png(p, "fig_many_curves.png", width = 6, height = 4)
}

# Slide 19b — missingness tile matrix ------------------------------------

{
  set.seed(3)
  n_pep <- 30
  n_sample <- 12
  # peptide-level NA probability, MNAR-shaped (high abundance = low NA)
  mu <- runif(n_pep, 8, 18)
  pna <- pmin(0.85, pmax(0, 1 - (mu - 7) / 12))
  obs <- vapply(seq_len(n_pep), function(i) {
    runif(n_sample) > pna[i]
  }, logical(n_sample))
  df <- expand.grid(sample = seq_len(n_sample), peptide = seq_len(n_pep))
  df$present <- as.vector(obs)

  # order peptides by NA rate (most missing at top)
  na_rate <- df |>
    group_by(peptide) |>
    summarise(na = mean(!present)) |>
    arrange(desc(na))
  df$peptide <- factor(df$peptide, levels = na_rate$peptide)

  p <- ggplot(df, aes(sample, peptide, fill = present)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = c(`TRUE` = mid, `FALSE` = "grey85"),
                      labels = c(`TRUE` = "observed", `FALSE` = "missing"),
                      name = NULL) +
    annotate("text", x = 6.5, y = -1.5,
             label = "non-random missingness",
             colour = "black", size = 4.5) +
    coord_cartesian(clip = "off") +
    theme_anatomy() +
    theme(legend.position = "top",
          plot.margin = margin(20, 20, 50, 20))

  save_png(p, "fig_missingness_tiles.png", width = 6, height = 4.5)
}

# Simulation recipe flowchart --------------------------------------------

{
  # --- top row: 4 step boxes ---
  top_y <- 0.74
  top_h <- 0.42
  top_w <- 0.21
  top_x <- seq(0.10, 0.90, length.out = 4)

  icon_y  <- top_y + 0.05
  label_y <- top_y - 0.13

  top_boxes <- data.frame(x = top_x, y = top_y)

  top_labels <- data.frame(
    x = top_x,
    y = label_y,
    label = c(
      "1. Fit distribution\nto pilot data",
      "2. Draw n control,\nn treatment\n(shifted by δ)",
      "3. Optional:\nknock out rᵢ\n(MNAR)",
      "4. Run test\nis p ≤ α?"
    )
  )

  # forward arrows between adjacent top boxes
  fwd <- data.frame(
    x    = head(top_x, -1) + top_w / 2 + 0.005,
    xend = tail(top_x, -1) - top_w / 2 - 0.005,
    y    = top_y,
    yend = top_y
  )

  # loop arrow from box 4 back to box 2
  loop_y <- top_y - top_h / 2 - 0.07
  loop_x_start <- top_x[4]
  loop_x_end   <- top_x[2]

  # --- icons ---

  # Icon 1: large bell curve filling most of the box
  bell_dx <- seq(-2.6, 2.6, length.out = 120)
  curve_w_half <- top_w * 0.40
  curve_h_max  <- 0.16
  icon1 <- data.frame(
    x = top_x[1] + bell_dx * (curve_w_half / 2.6),
    y = icon_y - curve_h_max * 0.45 + dnorm(bell_dx) * (curve_h_max / dnorm(0))
  )

  # Icon 2: two coloured groups of dots (control grey left, treatment red right)
  set.seed(7)
  ctrl_dx <- runif(5, -0.07, -0.018)
  ctrl_dy <- runif(5, -0.05, 0.05)
  trt_dx  <- runif(5,  0.018, 0.07)
  trt_dy  <- runif(5, -0.05, 0.05)
  icon2 <- rbind(
    data.frame(x = top_x[2] + ctrl_dx, y = icon_y + ctrl_dy, g = "ctrl"),
    data.frame(x = top_x[2] + trt_dx,  y = icon_y + trt_dy,  g = "trt")
  )

  # Icon 3: same dots shifted to box 3, with X through three of them
  icon3 <- icon2
  icon3$x <- icon3$x + (top_x[3] - top_x[2])
  knocked_idx <- c(2, 5, 8)
  icon3$knocked <- seq_len(nrow(icon3)) %in% knocked_idx
  knocked_pts <- icon3[icon3$knocked, ]
  x_half <- 0.014
  icon3_x1 <- data.frame(
    x = knocked_pts$x - x_half, xend = knocked_pts$x + x_half,
    y = knocked_pts$y - x_half, yend = knocked_pts$y + x_half
  )
  icon3_x2 <- data.frame(
    x = knocked_pts$x - x_half, xend = knocked_pts$x + x_half,
    y = knocked_pts$y + x_half, yend = knocked_pts$y - x_half
  )

  # Icon 4: decision diamond with '?'
  icon4 <- data.frame(
    x = top_x[4] + c(0, 0.055, 0, -0.055, 0),
    y = icon_y   + c(0.07, 0, -0.07, 0, 0.07)
  )

  # --- bottom row: box 5 (wide, with formula + explanation) ---
  bot_x <- 0.5
  bot_y <- 0.17
  bot_w <- 0.88
  bot_h <- 0.32

  # connector arrow from below top row to top of bottom box
  conn_x <- 0.5
  conn_y_start <- loop_y - 0.07
  conn_y_end   <- bot_y + bot_h / 2 + 0.005

  p <- ggplot() +
    # top boxes
    geom_tile(data = top_boxes, aes(x = x, y = y, width = top_w, height = top_h),
              fill = "white", colour = mid, linewidth = 0.6) +
    # forward arrows
    geom_segment(data = fwd, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = unit(0.18, "in"), type = "closed"),
                 colour = mid, linewidth = 0.8) +
    # loop arrow (three segments)
    annotate("segment", x = loop_x_start, xend = loop_x_start,
             y = top_y - top_h / 2, yend = loop_y,
             colour = accent, linewidth = 0.9) +
    annotate("segment", x = loop_x_start, xend = loop_x_end,
             y = loop_y, yend = loop_y,
             colour = accent, linewidth = 0.9) +
    annotate("segment", x = loop_x_end, xend = loop_x_end,
             y = loop_y, yend = top_y - top_h / 2,
             arrow = arrow(length = unit(0.2, "in"), type = "closed"),
             colour = accent, linewidth = 0.9) +
    annotate("text", x = (loop_x_start + loop_x_end) / 2,
             y = loop_y - 0.035, label = "repeat n_sim times",
             colour = accent, size = 6.5) +
    # connector arrow from top region down to bottom box
    annotate("segment", x = conn_x, xend = conn_x,
             y = conn_y_start, yend = conn_y_end,
             arrow = arrow(length = unit(0.22, "in"), type = "closed"),
             colour = mid, linewidth = 0.8) +
    annotate("text", x = conn_x + 0.02, y = (conn_y_start + conn_y_end) / 2,
             label = "after n_sim repeats", hjust = 0,
             colour = "black", size = 7) +
    # icons
    geom_path(data = icon1, aes(x, y), colour = accent, linewidth = 1.4) +
    geom_point(data = icon2, aes(x, y, colour = g), size = 4.5) +
    geom_point(data = icon3, aes(x, y, colour = g), size = 4.5) +
    geom_segment(data = icon3_x1, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "black", linewidth = 1.3) +
    geom_segment(data = icon3_x2, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "black", linewidth = 1.3) +
    scale_colour_manual(values = c(ctrl = mid, trt = accent), guide = "none") +
    geom_path(data = icon4, aes(x, y), colour = accent, linewidth = 1.1) +
    annotate("text", x = top_x[4], y = icon_y, label = "?",
             colour = accent, size = 10, fontface = "bold") +
    # text labels under top icons
    geom_text(data = top_labels, aes(x = x, y = y, label = label),
              size = 5.8, colour = "black", lineheight = 1.05) +
    # bottom box
    annotate("tile", x = bot_x, y = bot_y, width = bot_w, height = bot_h,
             fill = "white", colour = mid, linewidth = 0.6) +
    annotate("text", x = bot_x, y = bot_y + bot_h * 0.30,
             label = "5. Monte Carlo power estimate",
             colour = "black", size = 6.5, fontface = "bold") +
    annotate("text", x = bot_x, y = bot_y - bot_h * 0.02,
             label = "Power̂ᵢ(n, δ)  =  (1 / n_sim)  Σⱼ  𝟙{p⁽ʲ⁾ ≤ α}",
             colour = "black", size = 7) +
    annotate("text", x = bot_x, y = bot_y - bot_h * 0.34,
             label = "the fraction of the n_sim repeats in which the test rejected H₀",
             colour = "black", size = 9, fontface = "bold.italic") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_anatomy()

  save_png(p, "fig_simulation_recipe.png", width = 14, height = 8)
}

# Three tests icon column ------------------------------------------------

{
  # three rows stacked vertically (top → bottom)
  row_y <- c(0.85, 0.50, 0.15)
  icon_dy <-  0.05    # icon center slightly above row center
  name_dy <- -0.10    # name center below icon

  # 1. Wilcoxon — three ranked dots with rank labels
  rank_x  <- 0.58
  rank_dy <- c(-0.075, 0, 0.075)
  rank_dots <- data.frame(
    x = rank_x,
    y = row_y[1] + icon_dy + rank_dy
  )
  rank_labels <- data.frame(
    x = rank_x - 0.16,
    y = row_y[1] + icon_dy + rank_dy,
    label = c("1", "2", "3")
  )

  # 2. Bootstrap — circular arrow
  cx2 <- 0.5
  cy2 <- row_y[2] + icon_dy
  rr  <- 0.06
  t <- seq(0.35 * pi, 2.05 * pi, length.out = 80)
  loop_path <- data.frame(
    x = cx2 + rr * cos(t),
    y = cy2 + rr * sin(t)
  )
  ti <- length(t)
  arrow_seg <- data.frame(
    x    = cx2 + rr * cos(t[ti - 3]),
    xend = cx2 + rr * cos(t[ti]),
    y    = cy2 + rr * sin(t[ti - 3]),
    yend = cy2 + rr * sin(t[ti])
  )

  # 3. Bayes-factor — fraction visual P(D | H1) / P(D | H0)
  bay_x <- 0.5
  bay_y <- row_y[3] + icon_dy
  bay_top <- data.frame(x = bay_x, y = bay_y + 0.035,
                        label = "P(D ∣ H₁)")
  bay_bar <- data.frame(x = bay_x - 0.15, xend = bay_x + 0.15,
                        y = bay_y, yend = bay_y)
  bay_bot <- data.frame(x = bay_x, y = bay_y - 0.035,
                        label = "P(D ∣ H₀)")

  names_df <- data.frame(
    x = 0.5,
    y = row_y + name_dy,
    label = c("Wilcoxon rank-sum", "Bootstrap-t", "Bayes-factor t")
  )

  p <- ggplot() +
    # icon 1
    geom_point(data = rank_dots, aes(x, y), colour = mid, size = 12) +
    geom_text(data = rank_labels, aes(x, y, label = label),
              colour = accent, fontface = "bold", size = 14) +
    # icon 2
    geom_path(data = loop_path, aes(x, y),
              colour = accent, linewidth = 2.2) +
    geom_segment(data = arrow_seg, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = unit(0.3, "in"), type = "closed"),
                 colour = accent, linewidth = 2.2) +
    # icon 3
    geom_text(data = bay_top, aes(x, y, label = label),
              size = 11, colour = "black") +
    geom_segment(data = bay_bar, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "black", linewidth = 1.6) +
    geom_text(data = bay_bot, aes(x, y, label = label),
              size = 11, colour = "black") +
    # names underneath each icon
    geom_text(data = names_df, aes(x, y, label = label),
              size = 8, fontface = "bold") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_anatomy()

  save_png(p, "fig_three_tests.png", width = 6, height = 8)
}

# pepdiff workflow flow diagram ------------------------------------------
# CSV -> read_pepdiff() -> pepdiff_data -> compare() -> pepdiff_results -> plot()
# Rendered as four "state" nodes joined by three accent-red function arrows.

{
  # aspect correction: figure is wide (width:height ~ 14:4.8), so a horizontal
  # data offset must be scaled by kx to read as visually square as a vertical one.
  fig_w <- 14; fig_h <- 4.8
  kx <- fig_h / fig_w
  vx <- function(v) v * kx            # y-scale extent -> x-scale extent

  node_x <- c(0.11, 0.37, 0.63, 0.89)
  node_y <- 0.30
  box_w  <- 0.215
  box_h  <- 0.20
  icon_y <- 0.66

  # --- node boxes, names, and role sublabels ---
  boxes <- data.frame(x = node_x, y = node_y)
  names_df <- data.frame(
    x = node_x, y = node_y + 0.022,
    label = c("CSV", "pepdiff_data", "pepdiff_results", "figure")
  )
  roles_df <- data.frame(
    x = node_x, y = node_y - 0.045,
    label = c("input file", "S3 object", "S3 object", "ggplot")
  )

  # --- function arrows between adjacent nodes ---
  arr <- data.frame(
    x    = head(node_x, -1) + box_w / 2 + 0.006,
    xend = tail(node_x, -1) - box_w / 2 - 0.006,
    y    = node_y, yend = node_y
  )
  fn_df <- data.frame(
    x = (head(node_x, -1) + tail(node_x, -1)) / 2,
    y = node_y + 0.135,
    label = c("read_pepdiff()", "compare()", "plot()")
  )

  # --- icons ---------------------------------------------------------------

  # Icon 1: CSV page with folded corner + text lines
  c1 <- node_x[1]
  pg_hx <- vx(0.085); pg_hy <- 0.11; fold <- vx(0.045)
  page <- data.frame(
    x = c(c1 - pg_hx, c1 + pg_hx - fold, c1 + pg_hx, c1 + pg_hx, c1 - pg_hx),
    y = c(icon_y + pg_hy, icon_y + pg_hy, icon_y + pg_hy - 0.06,
          icon_y - pg_hy, icon_y - pg_hy)
  )
  page_fold <- data.frame(
    x = c(c1 + pg_hx - fold, c1 + pg_hx, c1 + pg_hx - fold),
    y = c(icon_y + pg_hy, icon_y + pg_hy - 0.06, icon_y + pg_hy - 0.06)
  )
  page_lines <- data.frame(
    x = c1 - pg_hx * 0.55, xend = c1 + pg_hx * 0.55,
    y = icon_y + c(0.03, 0, -0.03), yend = icon_y + c(0.03, 0, -0.03)
  )

  # Icon 2: tidy data table (3x3 grid, accent header row)
  c2 <- node_x[2]
  tb_hx <- vx(0.10); tb_hy <- 0.085
  tbl2_bg   <- data.frame(x = c2, y = icon_y)
  tbl2_head <- data.frame(x = c2, y = icon_y + tb_hy - 0.028)  # accent header band
  tbl2_v <- data.frame(x = c2 + c(-tb_hx/3, tb_hx/3), xend = c2 + c(-tb_hx/3, tb_hx/3),
                       y = icon_y - tb_hy, yend = icon_y + tb_hy)
  tbl2_h <- data.frame(x = c2 - tb_hx, xend = c2 + tb_hx,
                       y = icon_y + c(tb_hy - 0.057, tb_hy - 0.114),
                       yend = icon_y + c(tb_hy - 0.057, tb_hy - 0.114))

  # Icon 3: results table with one highlighted (significant) row + star
  c3 <- node_x[3]
  r3_hx <- vx(0.10); r3_hy <- 0.085
  tbl3_bg  <- data.frame(x = c3, y = icon_y)
  tbl3_hit <- data.frame(x = c3, y = icon_y)            # accent middle row = a hit
  tbl3_h <- data.frame(x = c3 - r3_hx, xend = c3 + r3_hx,
                       y = icon_y + c(tb_hy - 0.057, tb_hy - 0.114),
                       yend = icon_y + c(tb_hy - 0.057, tb_hy - 0.114))
  star <- data.frame(
    angle = seq(90, 90 + 360, length.out = 11)[-11] * pi / 180
  )
  star$r <- rep(c(0.04, 0.017), length.out = 10)
  star_x <- c3 + r3_hx + vx(0.075) + vx(star$r) * cos(star$angle)
  star_y <- icon_y + star$r * sin(star$angle)
  star_df <- data.frame(x = star_x, y = star_y)

  # Icon 4: volcano-style scatter (L axes + points, a few accent at top)
  c4 <- node_x[4]
  ax_hx <- vx(0.10); ax_hy <- 0.10
  axes <- data.frame(
    x    = c(c4 - ax_hx, c4 - ax_hx),
    xend = c(c4 - ax_hx, c4 + ax_hx),
    y    = c(icon_y - ax_hy, icon_y - ax_hy),
    yend = c(icon_y + ax_hy, icon_y - ax_hy)
  )
  set.seed(11)
  np <- 22
  px <- runif(np, -1, 1)
  py <- 1 - px^2 + rnorm(np, 0, 0.18)            # volcano cloud
  sig <- py > 0.85
  pts4 <- data.frame(
    x = c4 + vx(0.085) * px,
    y = icon_y - ax_hy * 0.55 + (ax_hy * 1.05) * pmin(pmax(py, -0.2), 1.2),
    sig = sig
  )

  p <- ggplot() +
    # arrows + function names
    geom_segment(data = arr, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = unit(0.20, "in"), type = "closed"),
                 colour = accent, linewidth = 1.0) +
    geom_text(data = fn_df, aes(x, y, label = label),
              colour = accent, family = "mono", fontface = "bold", size = 9.5) +
    # node boxes
    geom_tile(data = boxes, aes(x = x, y = y, width = box_w, height = box_h),
              fill = "white", colour = mid, linewidth = 0.7) +
    geom_text(data = names_df, aes(x, y, label = label),
              family = "mono", fontface = "bold", size = 6.6, colour = "black") +
    geom_text(data = roles_df, aes(x, y, label = label),
              fontface = "italic", size = 4.6, colour = mid) +
    # icon 1: CSV page
    geom_polygon(data = page, aes(x, y), fill = "white", colour = mid, linewidth = 0.7) +
    geom_polygon(data = page_fold, aes(x, y), fill = soft, colour = mid, linewidth = 0.5) +
    geom_segment(data = page_lines, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = mid, linewidth = 0.7) +
    # icon 2: data table
    geom_tile(data = tbl2_bg, aes(x, y, width = 2 * tb_hx, height = 2 * tb_hy),
              fill = "white", colour = mid, linewidth = 0.7) +
    geom_tile(data = tbl2_head, aes(x, y, width = 2 * tb_hx, height = 0.057),
              fill = accent, colour = NA, alpha = 0.85) +
    geom_segment(data = tbl2_v, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = mid, linewidth = 0.5) +
    geom_segment(data = tbl2_h, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = mid, linewidth = 0.5) +
    # icon 3: results table with a hit row + star
    geom_tile(data = tbl3_bg, aes(x, y, width = 2 * r3_hx, height = 2 * r3_hy),
              fill = "white", colour = mid, linewidth = 0.7) +
    geom_tile(data = tbl3_hit, aes(x, y, width = 2 * r3_hx, height = 0.057),
              fill = accent, colour = NA, alpha = 0.85) +
    geom_segment(data = tbl3_h, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = mid, linewidth = 0.5) +
    geom_polygon(data = star_df, aes(x, y), fill = accent, colour = NA) +
    # icon 4: volcano scatter
    geom_segment(data = axes, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = mid, linewidth = 0.7) +
    geom_point(data = pts4, aes(x, y, colour = sig), size = 2.6) +
    scale_colour_manual(values = c(`FALSE` = soft, `TRUE` = accent), guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.05, 0.92)) +
    theme_anatomy()

  save_png(p, "fig_pepdiff_workflow.png", width = fig_w, height = fig_h)
}

# per-peptide gamma GLM concept -----------------------------------------
# Two skewed gamma densities (ctrl/trt) spanning the full frame; the
# treatment shifts the mean. No border/cards/title — the curves breathe.

{
  mu_c <- 3.0; mu_t <- 5.4              # group means, data scale (1.8x fold change)
  shp  <- 2.2
  base_y <- 0.16
  curve_h <- 0.70
  mapX <- function(v) 0.06 + (0.97 - 0.06) * v / 10

  xseq <- seq(0.001, 10, length.out = 240)
  d_c <- dgamma(xseq, shape = shp, scale = mu_c / shp)
  d_t <- dgamma(xseq, shape = shp, scale = mu_t / shp)
  dmax <- max(d_c, d_t)
  to_y <- function(d) base_y + curve_h * d / dmax
  area_c <- data.frame(x = c(mapX(xseq[1]), mapX(xseq), mapX(xseq[length(xseq)])),
                       y = c(base_y, to_y(d_c), base_y))
  area_t <- data.frame(x = c(mapX(xseq[1]), mapX(xseq), mapX(xseq[length(xseq)])),
                       y = c(base_y, to_y(d_t), base_y))
  meanL <- data.frame(
    x = c(mapX(mu_c), mapX(mu_t)),
    yend = c(to_y(d_c[which.min(abs(xseq - mu_c))]),
             to_y(d_t[which.min(abs(xseq - mu_t))])),
    g = c("c", "t")
  )

  p_data <- ggplot() +
    geom_polygon(data = area_c, aes(x, y), fill = mid,    alpha = 0.20) +
    geom_polygon(data = area_t, aes(x, y), fill = accent, alpha = 0.20) +
    geom_path(data = area_c, aes(x, y), colour = mid,    linewidth = 1.3) +
    geom_path(data = area_t, aes(x, y), colour = accent, linewidth = 1.3) +
    geom_segment(data = meanL, aes(x = x, xend = x, y = base_y, yend = yend, colour = g),
                 linetype = "32", linewidth = 0.8) +
    annotate("segment", x = 0.05, xend = 0.98, y = base_y, yend = base_y,
             colour = mid, linewidth = 0.6) +
    annotate("text", x = 0.5, y = base_y - 0.08,
             label = "control vs treatment abundance", colour = mid, size = 6.5) +
    annotate("text", x = mapX(mu_c) - 0.008, y = base_y + 0.14,
             label = "mean\nμ_ctrl", colour = mid, fontface = "italic",
             size = 5.5, hjust = 1, lineheight = 0.9) +
    annotate("text", x = mapX(mu_t) + 0.008, y = base_y + 0.14,
             label = "mean\nμ_trt", colour = accent, fontface = "italic",
             size = 5.5, hjust = 0, lineheight = 0.9) +
    scale_colour_manual(values = c(c = mid, t = accent), guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.03, 0.92)) +
    theme_anatomy()

  save_png(p_data, "fig_pepdiff_glm_data.png", width = 12, height = 7)
}

# why fit a model — raw abundances -> effect size + significance ---------
# Motivation schematic: noisy per-peptide data through the GLM to a call.

{
  set.seed(3)
  # --- left: one peptide's noisy abundances, ctrl vs trt ---
  cx <- 0.135; tx <- 0.255
  ctrl_pts <- data.frame(x = cx + runif(6, -0.022, 0.022),
                         y = 0.40 + rnorm(6, 0, 0.055), g = "c")
  trt_pts  <- data.frame(x = tx + runif(6, -0.022, 0.022),
                         y = 0.60 + rnorm(6, 0, 0.065), g = "t")
  pts <- rbind(ctrl_pts, trt_pts)
  mean_seg <- data.frame(
    x = c(cx - 0.035, tx - 0.035), xend = c(cx + 0.035, tx + 0.035),
    y = c(0.40, 0.60), yend = c(0.40, 0.60), g = c("c", "t")
  )
  grp_lab <- data.frame(x = c(cx, tx), y = 0.24, label = c("ctrl", "trt"),
                        g = c("c", "t"))

  # --- middle: GLM box + flanking arrows ---
  glm_box <- data.frame(x = 0.525, y = 0.50)
  a1 <- data.frame(x = 0.33, xend = 0.435, y = 0.50, yend = 0.50)
  a2 <- data.frame(x = 0.615, xend = 0.70, y = 0.50, yend = 0.50)

  # --- right: result card ---
  res_card <- data.frame(x = c(0.71, 0.985, 0.985, 0.71),
                         y = c(0.30, 0.30, 0.70, 0.70))

  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.92,
             label = "raw abundances  →  effect size  +  significance, per peptide",
             colour = mid, size = 6.5) +
    # left data panel
    annotate("segment", x = 0.07, xend = 0.07, y = 0.22, yend = 0.78,
             colour = mid, linewidth = 0.6) +
    annotate("text", x = 0.065, y = 0.83, label = "abundance",
             colour = mid, size = 5, hjust = 0) +
    geom_segment(data = mean_seg, aes(x = x, xend = xend, y = y, yend = yend, colour = g),
                 linewidth = 1.6) +
    geom_point(data = pts, aes(x, y, colour = g), size = 5, alpha = 0.85) +
    geom_text(data = grp_lab, aes(x, y, label = label, colour = g),
              fontface = "bold", size = 6) +
    annotate("text", x = 0.195, y = 0.13, label = "one peptide, a few noisy values",
             colour = mid, size = 5) +
    # arrows
    geom_segment(data = a1, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = unit(0.18, "in"), type = "closed"),
                 colour = accent, linewidth = 1.1) +
    geom_segment(data = a2, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = unit(0.18, "in"), type = "closed"),
                 colour = accent, linewidth = 1.1) +
    # GLM box
    geom_tile(data = glm_box, aes(x, y, width = 0.18, height = 0.20),
              fill = "white", colour = accent, linewidth = 1.1) +
    annotate("text", x = 0.525, y = 0.535, label = "Gamma GLM",
             colour = accent, fontface = "bold", size = 7) +
    annotate("text", x = 0.525, y = 0.465, label = "one per peptide",
             colour = mid, fontface = "italic", size = 5) +
    # result card
    geom_polygon(data = res_card, aes(x, y), fill = "white", colour = mid, linewidth = 0.8) +
    annotate("text", x = 0.725, y = 0.615, label = "fold change",
             colour = mid, size = 5, hjust = 0) +
    annotate("text", x = 0.965, y = 0.615, label = "1.8×",
             colour = "black", fontface = "bold", size = 6.5, hjust = 1) +
    annotate("text", x = 0.725, y = 0.50, label = "p-value",
             colour = mid, size = 5, hjust = 0) +
    annotate("text", x = 0.965, y = 0.50, label = "0.02",
             colour = "black", fontface = "bold", size = 6.5, hjust = 1) +
    annotate("text", x = 0.725, y = 0.385, label = "FDR call",
             colour = mid, size = 5, hjust = 0) +
    annotate("text", x = 0.965, y = 0.385, label = "✓ significant",
             colour = accent, fontface = "bold", size = 6, hjust = 1) +
    scale_colour_manual(values = c(c = mid, t = accent), guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.05, 0.97)) +
    theme_anatomy()

  save_png(p, "fig_pepdiff_why.png", width = 13, height = 5.4)
}

# regression intuition — slope is the effect, scatter is the confidence ---
# Three panels: flat (no change); clear slope + tight (significant);
# same slope + loose (not significant).

{
  set.seed(8)
  panels <- list(
    list(cx = 0.18, cm = 0.55, tm = 0.55, sd = 0.030, line = soft,
         cap = "flat line — no change"),
    list(cx = 0.50, cm = 0.44, tm = 0.66, sd = 0.022, line = accent,
         cap = "clear slope, tight fit — significant"),
    list(cx = 0.82, cm = 0.44, tm = 0.66, sd = 0.075, line = mid,
         cap = "same slope, loose fit — not significant")
  )

  pts_all <- data.frame(); axis_all <- data.frame()
  glab_all <- data.frame(); cap_all <- data.frame()
  lines <- list()

  for (j in seq_along(panels)) {
    pn <- panels[[j]]
    ctrlx <- pn$cx - 0.06; trtx <- pn$cx + 0.06
    cy <- pmin(pmax(rnorm(5, pn$cm, pn$sd), 0.32), 0.78)
    ty <- pmin(pmax(rnorm(5, pn$tm, pn$sd), 0.32), 0.78)
    pts_all <- rbind(pts_all,
      data.frame(x = ctrlx + runif(5, -0.013, 0.013), y = cy, g = "c"),
      data.frame(x = trtx  + runif(5, -0.013, 0.013), y = ty, g = "t"))
    # axis (vertical) + baseline (horizontal)
    axis_all <- rbind(axis_all,
      data.frame(x = pn$cx - 0.115, xend = pn$cx - 0.115, y = 0.30, yend = 0.80),
      data.frame(x = pn$cx - 0.115, xend = pn$cx + 0.12,  y = 0.30, yend = 0.30))
    # regression line through the two group means, slightly extended
    slope <- (pn$tm - pn$cm) / (trtx - ctrlx)
    lines[[j]] <- list(
      x = ctrlx - 0.03, xend = trtx + 0.03,
      y = pn$cm - slope * 0.03, yend = pn$tm + slope * 0.03, col = pn$line)
    glab_all <- rbind(glab_all,
      data.frame(x = ctrlx, y = 0.255, label = "ctrl", g = "c"),
      data.frame(x = trtx,  y = 0.255, label = "trt",  g = "t"))
    cap_all <- rbind(cap_all,
      data.frame(x = pn$cx, y = 0.15, label = pn$cap))
  }

  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.93,
             label = "Fit a line through the groups: the slope is the change, the scatter is the confidence",
             colour = "black", fontface = "bold", size = 6.2) +
    annotate("text", x = 0.018, y = 0.56, label = "abundance",
             colour = mid, size = 5, angle = 90) +
    geom_segment(data = axis_all, aes(x = x, xend = xend, y = y, yend = yend),
                 colour = mid, linewidth = 0.6) +
    # regression lines (one per panel, colour encodes the verdict)
    annotate("segment", x = lines[[1]]$x, xend = lines[[1]]$xend,
             y = lines[[1]]$y, yend = lines[[1]]$yend,
             colour = lines[[1]]$col, linewidth = 2.0) +
    annotate("segment", x = lines[[2]]$x, xend = lines[[2]]$xend,
             y = lines[[2]]$y, yend = lines[[2]]$yend,
             colour = lines[[2]]$col, linewidth = 2.0) +
    annotate("segment", x = lines[[3]]$x, xend = lines[[3]]$xend,
             y = lines[[3]]$y, yend = lines[[3]]$yend,
             colour = lines[[3]]$col, linewidth = 2.0) +
    geom_point(data = pts_all, aes(x, y, colour = g), size = 4.2, alpha = 0.9) +
    geom_text(data = glab_all, aes(x, y, label = label, colour = g),
              fontface = "bold", size = 5.5) +
    geom_text(data = cap_all, aes(x, y, label = label),
              colour = "black", size = 5.2) +
    scale_colour_manual(values = c(c = mid, t = accent), guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.08, 0.98)) +
    theme_anatomy()

  save_png(p, "fig_pepdiff_regression.png", width = 14, height = 5.6)
}

message("done.")
