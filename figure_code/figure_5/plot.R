if (!requireNamespace("patchwork", quietly = TRUE)) webr::install("patchwork")

library(dplyr)
library(ggplot2)
library(patchwork)

JITTER_SEED <- 20250308L
MIN_PATIENTS_PER_GROUP <- 10L
MAX_COLLECT_BYTES <- 10 * 1024^2
NO_VIRUS <- "No virus detected"
OTHER_GROUP <- "Other"
LABEL_ANGLE <- 55
TITLE_OFFSET_PT <- 45

VIRUS_LEVELS <- c(
  "HPV", "EBV", "Merkel cell polyomavirus", "HBV",
  "HHV-8/KSHV", "CMV", "Other virus", "Multiple viruses"
)
FUSION_LEVELS <- c("0–20", "21–50", ">50")

# LaCroixColoR "Apricot" has six colours; interpolated to one per virus
# category, as lacroix_palette("Apricot", n, type = "continuous") does.
APRICOT <- c("#D72000", "#EE6100", "#FFAD0A", "#1BB6AF", "#9093A2", "#132157")
VIRUS_COLORS <- setNames(
  grDevices::colorRampPalette(APRICOT)(length(VIRUS_LEVELS)),
  VIRUS_LEVELS
)

# wesanderson "AsteroidCity3", first colours as wes_palette("AsteroidCity3", n).
ASTEROID_CITY_3 <- c("#FBA72A", "#D3D4D8", "#CB7A5C", "#5785C1")
FUSION_COLORS <- setNames(ASTEROID_CITY_3[seq_along(FUSION_LEVELS)], FUSION_LEVELS)


# ---- data -------------------------------------------------------------------

raw_data <- result |>
  select(
    datasource = `_datasource`, sampleId, donorId, tumor_category,
    t_cell_score, germline_telomere_length, somatic_telomere_length,
    virus_category, rna_fusion_count
  ) |>
  collect(max_bytes = MAX_COLLECT_BYTES)

if (as.numeric(object.size(raw_data)) > 0.95 * MAX_COLLECT_BYTES) {
  warning(
    "Collected data is close to MAX_COLLECT_BYTES; the cohort may be truncated.",
    call. = FALSE
  )
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

collapse_virus_categories <- function(x) {
  detected <- unique(x[!is.na(x) & x != NO_VIRUS])
  if (length(detected) == 0L) return(NO_VIRUS)
  if (length(detected) == 1L) return(detected)
  "Multiple viruses"
}

rare_locations <- raw_data |>
  distinct(datasource, donorId, tumor_category) |>
  count(tumor_category, name = "n_patients") |>
  filter(n_patients < MIN_PATIENTS_PER_GROUP) |>
  arrange(desc(n_patients), tumor_category)

raw_data <- raw_data |>
  mutate(
    tumor_group = if_else(
      tumor_category %in% rare_locations$tumor_category,
      OTHER_GROUP,
      tumor_category
    )
  )

patient_data <- raw_data |>
  group_by(datasource, donorId, tumor_group) |>
  summarise(
    t_cell_score = safe_mean(t_cell_score),
    germline_telomere_length = safe_mean(germline_telomere_length),
    somatic_telomere_length = safe_mean(somatic_telomere_length),
    virus = collapse_virus_categories(virus_category),
    .groups = "drop"
  )


# ---- tumour location labels -------------------------------------------------

group_labels <- patient_data |>
  count(tumor_group, name = "n_patients") |>
  arrange(tumor_group == OTHER_GROUP, desc(n_patients), tumor_group) |>
  mutate(
    tumor_label = sprintf(
      "%s (n=%d)",
      if_else(tumor_group == OTHER_GROUP, "Other*", tumor_group),
      n_patients
    )
  )

ordered_labels <- group_labels$tumor_label

add_tumor_label <- function(data) {
  data |>
    left_join(select(group_labels, tumor_group, tumor_label), by = "tumor_group") |>
    mutate(tumor_label = factor(tumor_label, levels = ordered_labels))
}

scale_x_tumor <- function(...) {
  scale_x_discrete(drop = FALSE, limits = ordered_labels, ...)
}

patient_data <- add_tumor_label(patient_data)

rna_sample_data <- raw_data |>
  filter(!is.na(rna_fusion_count)) |>
  add_tumor_label()


# ---- panel data -------------------------------------------------------------

paired_telomere_data <- patient_data |>
  filter(germline_telomere_length > 0, somatic_telomere_length > 0) |>
  mutate(log2_ratio = log2(somatic_telomere_length / germline_telomere_length))

telomere_means <- paired_telomere_data |>
  group_by(tumor_label) |>
  summarise(mean_log2_ratio = mean(log2_ratio), .groups = "drop")

virus_summary <- patient_data |>
  group_by(tumor_label) |>
  summarise(
    n_patients = n(),
    n_positive = sum(virus != NO_VIRUS),
    .groups = "drop"
  ) |>
  mutate(percent_label = sprintf("%.1f%%", 100 * n_positive / n_patients))

virus_composition <- patient_data |>
  filter(virus != NO_VIRUS) |>
  count(tumor_label, virus, name = "n_positive_patients") |>
  mutate(virus = factor(virus, levels = VIRUS_LEVELS))

rna_fusion_counts <- rna_sample_data |>
  mutate(
    fusion_category = case_when(
      rna_fusion_count <= 20 ~ "0–20",
      rna_fusion_count <= 50 ~ "21–50",
      .default = ">50"
    ),
    fusion_category = factor(fusion_category, levels = FUSION_LEVELS)
  ) |>
  count(tumor_label, fusion_category, name = "n_samples")

rna_sample_totals <- rna_sample_data |>
  group_by(tumor_label) |>
  summarise(n_rna_samples = n_distinct(sampleId), .groups = "drop") |>
  mutate(label = sprintf("(n=%d)", n_rna_samples))

# Annotation rows sit just above the tallest bar, so they track the data.
virus_label_y <- 1.04 * max(virus_summary$n_positive, 0)
fusion_label_y <- 1.04 * max(rna_sample_totals$n_rna_samples, 0)


# ---- theme ------------------------------------------------------------------

base_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 9), lineheight = 0.9),
    plot.caption.position = "plot",
    plot.caption = element_text(
      face = "bold",
      size = 10,
      hjust = 0,
      margin = margin(t = 5, l = TITLE_OFFSET_PT)
    ),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5),
    legend.key.height = grid::unit(3.5, "mm"),
    legend.key.width = grid::unit(4.5, "mm"),
    legend.box.margin = margin(0, 4, 0, 8),
    plot.margin = margin(3, 8, 3, 14)
  )

hide_x_axis <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

count_axis <- function(top_expansion) {
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 6),
    labels = scales::label_number(accuracy = 1),
    expand = expansion(mult = c(0, top_expansion))
  )
}


# ---- panels -----------------------------------------------------------------

panel_a <- patient_data |>
  filter(!is.na(t_cell_score)) |>
  ggplot(aes(x = tumor_label, y = t_cell_score)) +
  geom_boxplot(outlier.alpha = 0.22, width = 0.65, fill = "#4C78A8") +
  scale_x_tumor(position = "top") +
  scale_y_continuous(breaks = seq(0, 2, by = 0.5)) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(
    y = "PASS TRA/TRD + TRB\nassigned reads /\ntumor median coverage",
    caption = "A  Coverage-adjusted T-cell score"
  ) +
  base_theme +
  theme(
    axis.text.x.top = element_text(angle = LABEL_ANGLE, hjust = 0, vjust = 0, size = 7),
    plot.margin = margin(3, 8, 16, 14)
  )

panel_b <- ggplot(paired_telomere_data, aes(x = tumor_label, y = log2_ratio)) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_boxplot(outlier.shape = NA, width = 0.68, fill = "#D9E6F2", color = "#333333") +
  geom_point(
    position = position_jitter(width = 0.17, height = 0, seed = JITTER_SEED),
    alpha = 0.1,
    size = 0.38
  ) +
  geom_point(
    data = telomere_means,
    aes(y = mean_log2_ratio),
    shape = 23,
    size = 2.3,
    fill = "#D73027",
    color = "black",
    stroke = 0.3
  ) +
  scale_x_tumor() +
  scale_y_continuous(breaks = -5:5) +
  coord_cartesian(ylim = c(-5, 5)) +
  labs(
    y = "log2(somatic /\ngermline telomere length)",
    caption = "B  Paired tumor–normal telomere-length change — diamonds indicate group means"
  ) +
  base_theme +
  hide_x_axis +
  theme(plot.margin = margin(16, 8, 16, 14))

panel_c <- ggplot(
  virus_composition,
  aes(x = tumor_label, y = n_positive_patients, fill = virus)
) +
  geom_col(width = 0.68, show.legend = TRUE) +
  geom_text(
    data = virus_summary,
    aes(x = tumor_label, y = virus_label_y, label = percent_label),
    inherit.aes = FALSE,
    vjust = 0,
    size = 2.25
  ) +
  scale_x_tumor() +
  count_axis(top_expansion = 0.14) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(name = "Virus category", values = VIRUS_COLORS, drop = FALSE) +
  labs(
    y = "Virus-positive\npatients (n)",
    caption = "C  Virus presence — absolute counts; labels show prevalence"
  ) +
  base_theme +
  hide_x_axis +
  theme(plot.margin = margin(16, 8, 18, 14))

panel_d <- ggplot(
  rna_fusion_counts,
  aes(x = tumor_label, y = n_samples, fill = fusion_category)
) +
  geom_col(width = 0.7, show.legend = TRUE) +
  geom_text(
    data = rna_sample_totals,
    aes(x = tumor_label, y = fusion_label_y, label = label),
    inherit.aes = FALSE,
    angle = LABEL_ANGLE,
    hjust = 0,
    vjust = 0,
    size = 2.05
  ) +
  scale_x_tumor() +
  count_axis(top_expansion = 0.22) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(name = "RNA fusions per sample", values = FUSION_COLORS, drop = FALSE) +
  labs(
    y = "RNA-seq samples (n)",
    caption = "D  Expressed gene fusions — sample counts by fusion-count bin"
  ) +
  base_theme +
  hide_x_axis +
  theme(plot.margin = margin(18, 8, 3, 14))


# ---- figure -----------------------------------------------------------------

other_composition <- if (nrow(rare_locations) > 0L) {
  paste0(
    rare_locations$tumor_category, " (n=", rare_locations$n_patients, ")",
    collapse = "; "
  )
} else {
  "None"
}

figure_5 <- panel_a / panel_b / panel_c / panel_d +
  plot_layout(heights = c(1.25, 1.3, 1.2, 1.25), guides = "keep") +
  plot_annotation(
    title = "Figure 5. Immune, telomere, viral and RNA features by tumor location",
    caption = paste0("*Other composition: ", other_composition),
    theme = theme(
      plot.title = element_text(face = "bold", margin = margin(l = TITLE_OFFSET_PT))
    )
  )

print(figure_5)
