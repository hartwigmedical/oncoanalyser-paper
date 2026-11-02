library(dplyr)
library(ggplot2)

if (!requireNamespace("patchwork", quietly = TRUE)) {
  webr::install("patchwork")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Required package 'patchwork' could not be installed.")
}
library(patchwork)

JITTER_SEED <- 20250308L
MIN_PATIENTS_PER_GROUP <- 10L
NO_VIRUS <- "No virus detected"
OTHER_GROUP <- "Other"

raw_data <- result |>
  select(
    donorId,
    tumor_category,
    `t-cell score`,
    germline_telomere_length,
    somatic_telomere_length,
    virus_category
  ) |>
  rename(t_cell_score = `t-cell score`) |>
  collect(max_bytes = 8 * 1024^2)

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

collapse_virus_categories <- function(x) {
  detected <- unique(as.character(x[!is.na(x) & x != NO_VIRUS]))

  if (length(detected) == 0L) {
    NO_VIRUS
  } else if (length(detected) == 1L) {
    detected
  } else {
    "Multiple viruses"
  }
}

percent_labels <- function(x) paste0(x, "%")

# Combine tumor locations represented by fewer than 10 patients.
location_counts <- raw_data |>
  distinct(donorId, tumor_category) |>
  count(tumor_category, name = "n_patients") |>
  arrange(desc(n_patients), tumor_category)

rare_locations <- location_counts |>
  filter(n_patients < MIN_PATIENTS_PER_GROUP)

patient_data <- raw_data |>
  mutate(
    tumor_group = if_else(
      tumor_category %in% rare_locations$tumor_category,
      OTHER_GROUP,
      tumor_category
    )
  ) |>
  group_by(donorId, tumor_group) |>
  summarise(
    t_cell_score = safe_mean(t_cell_score),
    germline_telomere_length = safe_mean(germline_telomere_length),
    somatic_telomere_length = safe_mean(somatic_telomere_length),
    virus = collapse_virus_categories(virus_category),
    .groups = "drop"
  )

# Order tumor groups by decreasing patient count, with Other last.
group_counts <- patient_data |>
  distinct(donorId, tumor_group) |>
  count(tumor_group, name = "n_patients")

ordered_groups <- group_counts |>
  filter(tumor_group != OTHER_GROUP) |>
  arrange(desc(n_patients), tumor_group) |>
  pull(tumor_group)

if (OTHER_GROUP %in% group_counts$tumor_group) {
  ordered_groups <- c(ordered_groups, OTHER_GROUP)
}

group_labels <- group_counts |>
  mutate(
    tumor_label = if_else(
      tumor_group == OTHER_GROUP,
      sprintf("Other* (n=%d)", n_patients),
      sprintf("%s (n=%d)", tumor_group, n_patients)
    )
  )

ordered_labels <- group_labels$tumor_label[
  match(ordered_groups, group_labels$tumor_group)
]

patient_data <- patient_data |>
  left_join(group_labels, by = "tumor_group") |>
  mutate(tumor_label = factor(tumor_label, levels = ordered_labels))

# Paired tumor-normal telomere metrics.
paired_telomere_data <- patient_data |>
  filter(
    !is.na(germline_telomere_length),
    !is.na(somatic_telomere_length),
    germline_telomere_length > 0,
    somatic_telomere_length > 0
  ) |>
  mutate(
    log2_ratio = log2(somatic_telomere_length / germline_telomere_length),
    shorter_in_tumor = somatic_telomere_length < germline_telomere_length
  )

telomere_summary <- paired_telomere_data |>
  group_by(tumor_label) |>
  summarise(
    median_log2_ratio = median(log2_ratio),
    percent_shorter = 100 * mean(shorter_in_tumor),
    .groups = "drop"
  )

paired_telomere_data <- paired_telomere_data |>
  left_join(
    select(telomere_summary, tumor_label, median_log2_ratio),
    by = "tumor_label"
  )

# Virus bars use absolute counts; labels show prevalence percentages.
virus_summary <- patient_data |>
  mutate(is_positive = virus != NO_VIRUS) |>
  group_by(tumor_label) |>
  summarise(
    n_patients = n_distinct(donorId),
    n_positive = n_distinct(donorId[is_positive]),
    percent_positive = 100 * n_positive / n_patients,
    .groups = "drop"
  ) |>
  mutate(percent_label = sprintf("%.1f%%", percent_positive))

virus_composition <- patient_data |>
  filter(virus != NO_VIRUS) |>
  count(tumor_label, virus, name = "n_positive_patients")

virus_levels <- c(
  "HPV", "EBV", "Merkel cell polyomavirus", "HBV",
  "HHV-8/KSHV", "CMV", "Other virus", "Multiple viruses"
)

virus_colors <- c(
  "HPV" = "#F28E2B",
  "EBV" = "#E15759",
  "Merkel cell polyomavirus" = "#B07AA1",
  "HBV" = "#59A14F",
  "HHV-8/KSHV" = "#9C755F",
  "CMV" = "#76B7B2",
  "Other virus" = "#EDC948",
  "Multiple viruses" = "#4E79A7"
)

virus_composition <- virus_composition |>
  mutate(virus = factor(virus, levels = virus_levels))

# Legends are kept with their own panels. This prevents guides belonging to
# panels B and D from being merged into, or overlapping within, one guide box.
base_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 10.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 9), lineheight = 0.9),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7.5),
    legend.key.height = grid::unit(3.5, "mm"),
    legend.key.width = grid::unit(4.5, "mm"),
    legend.box.margin = margin(0, 4, 0, 8),
    plot.margin = margin(3, 8, 3, 14)
  )

hide_x_axis <- theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)

panel_a <- ggplot(
  filter(patient_data, !is.na(t_cell_score)),
  aes(x = tumor_label, y = t_cell_score)
) +
  geom_boxplot(outlier.alpha = 0.22, width = 0.65, fill = "#4C78A8") +
  scale_x_discrete(position = "top", drop = FALSE, limits = ordered_labels) +
  scale_y_continuous(breaks = seq(0, 2, by = 0.5)) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(
    title = "A  Coverage-adjusted T-cell score",
    y = "PASS TRA/TRD + TRB\nassigned reads /\ntumor median coverage"
  ) +
  base_theme +
  theme(
    axis.text.x.top = element_text(
      angle = 55, hjust = 0, vjust = 0, size = 7
    )
  )

panel_b1 <- ggplot(
  paired_telomere_data,
  aes(x = tumor_label, y = log2_ratio)
) +
  geom_hline(yintercept = 0) +
  geom_boxplot(
    aes(fill = median_log2_ratio),
    outlier.shape = NA,
    width = 0.68
  ) +
  geom_point(
    position = position_jitter(width = 0.17, height = 0, seed = JITTER_SEED),
    alpha = 0.1,
    size = 0.38
  ) +
  scale_x_discrete(drop = FALSE, limits = ordered_labels) +
  scale_y_continuous(breaks = -5:5) +
  coord_cartesian(ylim = c(-5, 5)) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Median log2 ratio"
  ) +
  labs(
    title = "B1  Patient-level tumor-versus-normal change",
    y = "log2(somatic /\ngermline\ntelomere\nlength)"
  ) +
  base_theme +
  hide_x_axis

panel_b2 <- ggplot(
  telomere_summary,
  aes(x = tumor_label, y = percent_shorter)
) +
  geom_col(width = 0.68, fill = "#4C78A8") +
  geom_text(
    aes(label = sprintf("%.1f", percent_shorter)),
    vjust = -0.2,
    size = 2.25
  ) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  scale_x_discrete(drop = FALSE, limits = ordered_labels) +
  scale_y_continuous(
    limits = c(0, 108),
    breaks = c(0, 25, 50, 75, 100),
    labels = percent_labels,
    expand = c(0, 0)
  ) +
  labs(
    title = "B2  Patients with shorter tumor telomeres",
    y = "Patients with\nsomatic < germline\ntelomere\nlength"
  ) +
  base_theme +
  hide_x_axis

# Panel D uses a standard linear y-axis.
panel_d <- ggplot(
  virus_composition,
  aes(x = tumor_label, y = n_positive_patients, fill = virus)
) +
  geom_col(width = 0.68) +
  geom_text(
    data = virus_summary,
    aes(x = tumor_label, y = n_positive, label = percent_label),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 2.25,
    fontface = "bold"
  ) +
  scale_x_discrete(drop = FALSE, limits = ordered_labels) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 5),
    labels = scales::label_number(accuracy = 1),
    expand = expansion(mult = c(0, 0.12))
  ) +
  scale_fill_manual(
    name = "Virus category",
    values = virus_colors,
    drop = TRUE
  ) +
  labs(
    title = "D  Virus-positive patients: absolute counts",
    y = "Virus-positive\npatients (n)"
  ) +
  base_theme +
  hide_x_axis

# Keep the B1 legend inside the B section rather than collecting it globally.
panel_b <- (
  panel_b1 / panel_b2 +
    plot_layout(heights = c(1.35, 0.8), guides = "keep")
) +
  plot_annotation(title = "B  Paired tumor–normal telomere-length change")

other_composition <- if (nrow(rare_locations) > 0L) {
  paste0(
    rare_locations$tumor_category,
    " (n=", rare_locations$n_patients, ")",
    collapse = "; "
  )
} else {
  "None"
}

# guides = "keep" leaves each legend in the row of its corresponding panel,
# preventing the telomere and virus legends from overlapping.
final_figure <- (
  panel_a / panel_b / panel_d +
    plot_layout(
      heights = c(1.3, 2.2, 1.2),
      guides = "keep"
    ) +
    plot_annotation(
      title = paste(
        "Coverage-adjusted T-cell score, paired telomere change,",
        "and viruses by tumor location"
      ),
      subtitle = paste(
        "Virus bars show absolute positive-patient counts on a linear scale;",
        "labels show prevalence among all patients."
      ),
      caption = paste0("*Other composition: ", other_composition)
    )
) &
  theme(plot.title = element_text(face = "bold"))

print(final_figure)

