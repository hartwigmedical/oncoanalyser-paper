library(dplyr)
library(ggplot2)

if (!requireNamespace("patchwork", quietly = TRUE)) webr::install("patchwork")
if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required")
library(patchwork)

NO_VIRUS <- "No virus detected"
OTHER_GROUP <- "Other"
TITLE_OFFSET_PT <- 45
LABEL_ANGLE <- 55

x <- result |>
  select(
    sampleId, donorId, tumor_category, t_cell_score,
    germline_telomere_length, somatic_telomere_length,
    virus_category, rna_fusion_count
  ) |>
  collect(max_bytes = 12 * 1024^2)

safe_mean <- function(z) if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
collapse_virus <- function(z) {
  detected <- unique(as.character(z[!is.na(z) & z != NO_VIRUS]))
  if (!length(detected)) NO_VIRUS else if (length(detected) == 1L) detected else "Multiple viruses"
}

location_counts <- x |>
  distinct(donorId, tumor_category) |>
  count(tumor_category, name = "n_patients") |>
  arrange(desc(n_patients), tumor_category)

rare_locations <- location_counts |>
  filter(n_patients < 10)

x <- x |>
  mutate(tumor_group = if_else(
    tumor_category %in% rare_locations$tumor_category,
    OTHER_GROUP,
    tumor_category
  ))

patient_data <- x |>
  group_by(donorId, tumor_group) |>
  summarise(
    t_cell_score = safe_mean(t_cell_score),
    germline_telomere_length = safe_mean(germline_telomere_length),
    somatic_telomere_length = safe_mean(somatic_telomere_length),
    virus = collapse_virus(virus_category),
    .groups = "drop"
  )

group_counts <- patient_data |>
  count(tumor_group, name = "n_patients")

ordered_groups <- group_counts |>
  filter(tumor_group != OTHER_GROUP) |>
  arrange(desc(n_patients), tumor_group) |>
  pull(tumor_group)

if (OTHER_GROUP %in% group_counts$tumor_group) ordered_groups <- c(ordered_groups, OTHER_GROUP)

group_labels <- group_counts |>
  mutate(tumor_label = if_else(
    tumor_group == OTHER_GROUP,
    sprintf("Other* (n=%d)", n_patients),
    sprintf("%s (n=%d)", tumor_group, n_patients)
  ))

ordered_labels <- group_labels$tumor_label[match(ordered_groups, group_labels$tumor_group)]

patient_data <- patient_data |>
  left_join(group_labels, by = "tumor_group") |>
  mutate(tumor_label = factor(tumor_label, levels = ordered_labels))

rna_data <- x |>
  filter(!is.na(rna_fusion_count)) |>
  left_join(group_labels, by = "tumor_group") |>
  mutate(tumor_label = factor(tumor_label, levels = ordered_labels))

paired_telomere_data <- patient_data |>
  filter(
    !is.na(germline_telomere_length), !is.na(somatic_telomere_length),
    germline_telomere_length > 0, somatic_telomere_length > 0
  ) |>
  mutate(log2_ratio = log2(somatic_telomere_length / germline_telomere_length))

telomere_means <- paired_telomere_data |>
  group_by(tumor_label) |>
  summarise(mean_log2_ratio = mean(log2_ratio), .groups = "drop")

virus_summary <- patient_data |>
  mutate(is_positive = virus != NO_VIRUS) |>
  group_by(tumor_label) |>
  summarise(
    n_patients = n_distinct(donorId),
    n_positive = n_distinct(donorId[is_positive]),
    percent_positive = 100 * n_positive / n_patients,
    .groups = "drop"
  ) |>
  mutate(prevalence_label = sprintf("%.1f", percent_positive))

virus_data <- patient_data |>
  filter(virus != NO_VIRUS) |>
  count(tumor_label, virus, name = "n_positive") |>
  left_join(select(virus_summary, tumor_label, n_patients), by = "tumor_label") |>
  mutate(percent = 100 * n_positive / n_patients)

virus_levels <- c(
  "HPV", "EBV", "Merkel cell polyomavirus", "HBV",
  "HHV-8/KSHV", "CMV", "Other virus", "Multiple viruses"
)
virus_data <- virus_data |>
  mutate(virus = factor(virus, levels = virus_levels))

virus_colors <- c(
  "HPV" = "#F28E2B", "EBV" = "#E15759", "Merkel cell polyomavirus" = "#B07AA1",
  "HBV" = "#59A14F", "HHV-8/KSHV" = "#9C755F", "CMV" = "#76B7B2",
  "Other virus" = "#EDC948", "Multiple viruses" = "#4E79A7"
)

fusion_counts <- rna_data |>
  mutate(
    fusion_bin = case_when(
      rna_fusion_count <= 20 ~ "0–20",
      rna_fusion_count <= 50 ~ "21–50",
      TRUE ~ ">50"
    ),
    fusion_bin = factor(fusion_bin, levels = c("0–20", "21–50", ">50"))
  ) |>
  count(tumor_label, fusion_bin, name = "n_samples")

rna_totals <- rna_data |>
  group_by(tumor_label) |>
  summarise(n_rna_samples = n_distinct(sampleId), .groups = "drop") |>
  mutate(label = sprintf("(n=%d)", n_rna_samples))

base_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "plain", margin = margin(r = 9)),
    axis.text = element_text(face = "plain"),
    plot.caption.position = "plot",
    plot.caption = element_text(
      face = "bold", size = 10, hjust = 0,
      margin = margin(t = 5, l = TITLE_OFFSET_PT)
    ),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "plain", size = 8),
    legend.text = element_text(face = "plain", size = 7.5),
    legend.key.height = grid::unit(3.5, "mm"),
    legend.key.width = grid::unit(4.5, "mm")
  )

hide_x_axis <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# A small sigma makes this pseudo-log scale behave logarithmically above 0.01,
# expanding low scores and compressing high scores while still showing zero.
panel_a <- ggplot(
  filter(patient_data, !is.na(t_cell_score)),
  aes(tumor_label, t_cell_score)
) +
  geom_boxplot(width = 0.65, fill = "#4C78A8", outlier.alpha = 0.22) +
  scale_x_discrete(position = "top", limits = ordered_labels, drop = FALSE) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 0.01),
    breaks = c(0, 0.01, 0.03, 0.1, 0.3, 1, 2),
    labels = c("0", "0.01", "0.03", "0.1", "0.3", "1", "2")
  ) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(
    y = "PASS TRA/TRD + TRB\nassigned reads /\ntumor median coverage\n(pseudo-log scale)",
    caption = "A  Coverage-adjusted T-cell score"
  ) +
  base_theme +
  theme(
    axis.text.x.top = element_text(
      angle = LABEL_ANGLE, hjust = 0, vjust = 0, size = 7, face = "plain"
    ),
    plot.margin = margin(3, 8, 16, 14)
  )

panel_b <- ggplot(paired_telomere_data, aes(tumor_label, log2_ratio)) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_boxplot(width = 0.68, fill = "#D9E6F2", color = "#333333", outlier.alpha = 0.15) +
  geom_point(
    data = telomere_means,
    aes(tumor_label, mean_log2_ratio),
    inherit.aes = FALSE,
    shape = 23, size = 2.3, fill = "#D73027", color = "black", stroke = 0.3
  ) +
  scale_x_discrete(limits = ordered_labels, drop = FALSE) +
  scale_y_continuous(breaks = -5:5) +
  coord_cartesian(ylim = c(-5, 5)) +
  labs(
    y = "log2(somatic /\ngermline telomere length)",
    caption = "B  Paired tumor–normal telomere-length change — diamonds indicate group means"
  ) +
  base_theme + hide_x_axis +
  theme(plot.margin = margin(16, 8, 16, 14))

panel_c <- ggplot(virus_data, aes(tumor_label, percent, fill = virus)) +
  geom_col(width = 0.68) +
  geom_text(
    data = virus_summary,
    aes(tumor_label, y = 100, label = prevalence_label),
    inherit.aes = FALSE,
    vjust = -0.35, size = 2.25, fontface = "plain"
  ) +
  scale_x_discrete(limits = ordered_labels, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
    labels = function(z) paste0(z, "%"), expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  scale_fill_manual(name = "Virus category", values = virus_colors, drop = TRUE) +
  labs(
    y = "Virus-positive\npatients (%)",
    caption = "C  Virus presence — percentage of patients; values show prevalence"
  ) +
  base_theme + hide_x_axis +
  theme(plot.margin = margin(16, 8, 18, 14))

panel_d <- ggplot(fusion_counts, aes(tumor_label, n_samples, fill = fusion_bin)) +
  geom_col(width = 0.7) +
  geom_text(
    data = rna_totals,
    aes(tumor_label, y = 800, label = label),
    inherit.aes = FALSE,
    angle = LABEL_ANGLE, hjust = 0, vjust = 0, size = 2.05, fontface = "plain"
  ) +
  scale_x_discrete(limits = ordered_labels, drop = FALSE) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 6),
    labels = scales::label_number(accuracy = 1),
    expand = expansion(mult = c(0, 0.22))
  ) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(
    name = "RNA fusions per sample",
    values = c("0–20" = "#FEE8C8", "21–50" = "#FDBB84", ">50" = "#E34A33"),
    drop = FALSE
  ) +
  labs(
    y = "RNA-seq samples (n)",
    caption = "D  RNA-detected gene fusions — sample counts by fusion-count bin"
  ) +
  base_theme + hide_x_axis +
  theme(plot.margin = margin(18, 8, 3, 14))

other_composition <- if (nrow(rare_locations) > 0L) {
  paste0(
    rare_locations$tumor_category,
    " (n=", rare_locations$n_patients, ")",
    collapse = "; "
  )
} else {
  "None"
}

figure_5 <- (
  panel_a / panel_b / panel_c / panel_d +
    plot_layout(heights = c(1.25, 1.3, 1.2, 1.25), guides = "keep") +
    plot_annotation(
      title = "Figure 5. Immune, telomere, viral and RNA features by tumor location",
      caption = paste0("*Other composition: ", other_composition),
      theme = theme(
        plot.title = element_text(face = "bold", margin = margin(l = TITLE_OFFSET_PT))
      )
    )
)

print(figure_5)
