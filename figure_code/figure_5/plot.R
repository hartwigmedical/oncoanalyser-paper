library(dplyr)
library(ggplot2)
if (!requireNamespace("patchwork", quietly = TRUE)) webr::install("patchwork")
if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required")
library(patchwork)

NO_VIRUS <- "No virus detected"
OTHER_GROUP <- "Other"
TITLE_OFFSET <- 45
ANGLE <- 55

x <- result |>
  select(
    donorId, tumor_category, t_cell_score,
    germline_telomere_length, somatic_telomere_length,
    virus_category, dna_fusion_count, rna_fusion_count
  ) |>
  collect(max_bytes = 12 * 1024^2)

location_counts <- x |>
  distinct(donorId, tumor_category) |>
  count(tumor_category, name = "n_patients") |>
  arrange(desc(n_patients), tumor_category)
rare <- filter(location_counts, n_patients < 10)

x <- x |>
  mutate(tumor_group = if_else(
    tumor_category %in% rare$tumor_category,
    OTHER_GROUP,
    tumor_category
  ))

group_counts <- x |>
  distinct(donorId, tumor_group) |>
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

x <- x |>
  left_join(group_labels, by = "tumor_group") |>
  mutate(tumor_label = factor(tumor_label, levels = ordered_labels))

paired <- x |>
  filter(
    !is.na(germline_telomere_length), !is.na(somatic_telomere_length),
    germline_telomere_length > 0, somatic_telomere_length > 0
  ) |>
  mutate(log2_ratio = log2(somatic_telomere_length / germline_telomere_length))

# Median indicators for Panels A and B.
t_cell_medians <- x |>
  filter(!is.na(t_cell_score)) |>
  group_by(tumor_label) |>
  summarise(median_score = median(t_cell_score), .groups = "drop")

telomere_medians <- paired |>
  group_by(tumor_label) |>
  summarise(median_ratio = median(log2_ratio), .groups = "drop")

virus_summary <- x |>
  mutate(positive = virus_category != NO_VIRUS) |>
  group_by(tumor_label) |>
  summarise(
    total = n_distinct(donorId),
    positive = n_distinct(donorId[positive]),
    prevalence = 100 * positive / total,
    .groups = "drop"
  ) |>
  mutate(label = sprintf("%.1f", prevalence))

virus_data <- x |>
  filter(virus_category != NO_VIRUS) |>
  count(tumor_label, virus_category, name = "n") |>
  left_join(select(virus_summary, tumor_label, total), by = "tumor_label") |>
  mutate(percent = 100 * n / total)

virus_colors <- c(
  "HPV" = "#F28E2B", "EBV" = "#E15759", "Merkel cell polyomavirus" = "#B07AA1",
  "HBV" = "#59A14F", "HHV-8/KSHV" = "#9C755F", "CMV" = "#76B7B2",
  "Other virus" = "#EDC948", "Multiple viruses" = "#4E79A7"
)

fusion_donors <- x |>
  filter(!is.na(rna_fusion_count))

fusion_long <- bind_rows(
  transmute(
    fusion_donors,
    donorId,
    tumor_label,
    mechanism = "DNA-all Linx",
    count = dna_fusion_count
  ),
  transmute(
    fusion_donors,
    donorId,
    tumor_label,
    mechanism = "RNA-detected",
    count = rna_fusion_count
  )
) |>
  mutate(
    bin = case_when(
      count <= 20 ~ "0–20",
      count <= 50 ~ "21–50",
      TRUE ~ ">50"
    ),
    bin = factor(bin, levels = c("0–20", "21–50", ">50")),
    mechanism = factor(
      mechanism,
      levels = c("DNA-all Linx", "RNA-detected")
    )
  )

fusion_dist <- fusion_long |>
  count(tumor_label, mechanism, bin, name = "n_donors") |>
  group_by(tumor_label, mechanism) |>
  mutate(percent = 100 * n_donors / sum(n_donors)) |>
  ungroup() |>
  mutate(
    tumor_index = as.numeric(tumor_label),
    x_pos = tumor_index + if_else(mechanism == "DNA-all Linx", -0.2, 0.2)
  )

fusion_totals <- fusion_donors |>
  group_by(tumor_label) |>
  summarise(n = n_distinct(donorId), .groups = "drop") |>
  mutate(
    tumor_index = as.numeric(tumor_label),
    label = sprintf("(n=%d)", n)
  )

base_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(face = "plain"),
    plot.caption.position = "plot",
    plot.caption = element_text(
      face = "bold", size = 10, hjust = 0,
      margin = margin(t = 5, l = TITLE_OFFSET)
    ),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "plain", size = 8),
    legend.text = element_text(face = "plain", size = 7.5)
  )
hide_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Bmlunge boxplot styling: opaque pale fill, turquoise median diamonds,
# and semi-transparent blue outliers.
p_a <- ggplot(filter(x, !is.na(t_cell_score)), aes(tumor_label, t_cell_score)) +
  geom_boxplot(
    width = 0.65,
    fill = "#DAEAEC",
    alpha = 1,
    outlier.colour = "#59A7B4",
    outlier.alpha = 0.5
  ) +
  geom_point(
    data = t_cell_medians,
    aes(tumor_label, median_score),
    inherit.aes = FALSE,
    shape = 18,
    size = 2.2,
    color = "#76CFD8"
  ) +
  scale_x_discrete(position = "top", limits = ordered_labels, drop = FALSE) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 0.01),
    breaks = c(0, 0.01, 0.03, 0.1, 0.3, 1, 2),
    labels = c("0", "0.01", "0.03", "0.1", "0.3", "1", "2")
  ) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(
    y = "PASS TRA/TRD + TRB\nassigned reads /\ntumor median coverage\n(pseudo-log scale)",
    caption = "A  Coverage-adjusted T-cell score — diamonds indicate medians"
  ) +
  base_theme +
  theme(
    axis.text.x.top = element_text(angle = ANGLE, hjust = 0, vjust = 0, size = 7),
    plot.margin = margin(3, 8, 16, 14)
  )

p_b <- ggplot(paired, aes(tumor_label, log2_ratio)) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_boxplot(
    width = 0.68,
    fill = "#DAEAEC",
    alpha = 1,
    outlier.colour = "#59A7B4",
    outlier.alpha = 0.5
  ) +
  geom_point(
    data = telomere_medians,
    aes(tumor_label, median_ratio),
    inherit.aes = FALSE,
    shape = 18,
    size = 2.2,
    color = "#76CFD8"
  ) +
  scale_x_discrete(limits = ordered_labels, drop = FALSE) +
  coord_cartesian(ylim = c(-5, 5)) +
  labs(
    y = "log2(somatic /\ngermline telomere length)",
    caption = "B  Paired tumor–normal telomere-length change — diamonds indicate medians"
  ) +
  base_theme + hide_x +
  theme(plot.margin = margin(16, 8, 16, 14))

# Labels sit immediately above each stacked prevalence bar.
p_c <- ggplot(virus_data, aes(tumor_label, percent, fill = virus_category)) +
  geom_col(width = 0.68) +
  geom_text(
    data = virus_summary,
    aes(tumor_label, y = prevalence, label = label),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 2.25
  ) +
  scale_x_discrete(limits = ordered_labels, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(z) paste0(z, "%"),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  scale_fill_manual(name = "Virus category", values = virus_colors) +
  labs(
    y = "Virus-positive\ndonors (%)",
    caption = "C  Virus presence — percentage of donors; values show prevalence"
  ) +
  base_theme + hide_x +
  theme(plot.margin = margin(16, 8, 18, 14))

# Bmlunge palette with subtle mechanism outlines.
p_d <- ggplot(
  fusion_dist,
  aes(x_pos, percent, fill = bin, color = mechanism)
) +
  geom_col(width = 0.34, linewidth = 0.3) +
  geom_text(
    data = fusion_totals,
    aes(tumor_index, y = 100, label = label),
    inherit.aes = FALSE,
    angle = ANGLE,
    hjust = 0,
    vjust = 0,
    size = 2.05
  ) +
  scale_x_continuous(
    limits = c(0.5, length(ordered_labels) + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(z) paste0(z, "%"),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  scale_fill_manual(
    name = "Fusions per donor",
    values = c(
      "0–20" = "#DAEAEC",
      "21–50" = "#76CFD8",
      ">50" = "#59A7B4"
    ),
    drop = FALSE
  ) +
  scale_color_manual(
    name = "Detection mechanism",
    values = c(
      "DNA-all Linx" = "#333333",
      "RNA-detected" = "#FFFFFF"
    ),
    drop = FALSE
  ) +
  labs(
    y = "Donors with RNA data (%)",
    caption = "D  Donor-level DNA versus RNA fusion-count distributions"
  ) +
  base_theme + hide_x +
  theme(plot.margin = margin(18, 8, 3, 14))

other <- if (nrow(rare)) {
  paste0(rare$tumor_category, " (n=", rare$n_patients, ")", collapse = "; ")
} else {
  "None"
}

fig <- p_a / p_b / p_c / p_d +
  plot_layout(heights = c(1.25, 1.3, 1.2, 1.25), guides = "keep") +
  plot_annotation(
    title = "Figure 5. Immune, telomere, viral and RNA features by tumor location",
    caption = paste0("*Other composition: ", other),
    theme = theme(
      plot.title = element_text(face = "bold", margin = margin(l = TITLE_OFFSET))
    )
  )

print(fig)
