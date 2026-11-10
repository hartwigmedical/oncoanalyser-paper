if (!requireNamespace("patchwork", quietly = TRUE)) webr::install("patchwork")

library(dplyr)
library(ggplot2)
library(patchwork)

MIN_DONORS_PER_GROUP <- 10L
NO_VIRUS <- "No virus detected"
OTHER_GROUP <- "Other"
TITLE_OFFSET <- 45
LABEL_ANGLE <- 55

BOX_FILL <- "#DAEAEC"
BOX_OUTLIER <- "#59A7B4"
MEDIAN_COLOR <- "#76CFD8"

VIRUS_LEVELS <- c(
  "HPV", "EBV", "Merkel cell polyomavirus", "HBV",
  "HHV-8/KSHV", "CMV", "Other virus", "Multiple viruses"
)

# tvthemes "Aquamarine" has seven colours; interpolated to one per virus
# category, as stevenUniverse_pal("Aquamarine", n, type = "continuous") does.
AQUAMARINE <- c(
  "#7ff9f7", "#0b2483", "#0f8cde", "#4c4cff", "#040509", "#23b4dc", "#3c7084"
)
VIRUS_COLORS <- setNames(
  grDevices::colorRampPalette(AQUAMARINE)(length(VIRUS_LEVELS)),
  VIRUS_LEVELS
)

FUSION_BINS <- c("0–20", "21–50", ">50")
FUSION_COLORS <- setNames(c("#DAEAEC", "#76CFD8", "#59A7B4"), FUSION_BINS)

DNA_FUSIONS <- "DNA-detected"
RNA_FUSIONS <- "RNA-detected"
MECHANISM_COLORS <- setNames(c("#333333", "#FFFFFF"), c(DNA_FUSIONS, RNA_FUSIONS))


# ---- donor-level data -------------------------------------------------------

donors <- result |>
  select(
    datasource = `_datasource`, donorId, tumor_category,
    t_cell_score, germline_telomere_length, somatic_telomere_length,
    virus_category, dna_fusion_count, rna_fusion_count
  ) |>
  collect(max_bytes = 12 * 1024^2)

if (anyDuplicated(select(donors, datasource, donorId)) > 0) {
  stop("selection.sql should return one row per donor.", call. = FALSE)
}


# ---- tumour location groups -------------------------------------------------

rare_locations <- donors |>
  count(tumor_category, name = "n_donors") |>
  filter(n_donors < MIN_DONORS_PER_GROUP) |>
  arrange(desc(n_donors), tumor_category)

donors <- donors |>
  mutate(
    tumor_group = if_else(
      tumor_category %in% rare_locations$tumor_category,
      OTHER_GROUP,
      tumor_category
    )
  )

group_labels <- donors |>
  count(tumor_group, name = "n_donors") |>
  arrange(tumor_group == OTHER_GROUP, desc(n_donors), tumor_group) |>
  mutate(
    tumor_label = sprintf(
      "%s (n=%d)",
      if_else(tumor_group == OTHER_GROUP, "Other*", tumor_group),
      n_donors
    )
  )

ordered_labels <- group_labels$tumor_label
n_groups <- length(ordered_labels)

donors <- donors |>
  left_join(select(group_labels, tumor_group, tumor_label), by = "tumor_group") |>
  mutate(tumor_label = factor(tumor_label, levels = ordered_labels))


# ---- panel data -------------------------------------------------------------

group_medians <- function(data, value) {
  data |>
    group_by(tumor_label) |>
    summarise(median = median({{ value }}), .groups = "drop")
}

t_cell <- filter(donors, !is.na(t_cell_score))

telomere <- donors |>
  filter(germline_telomere_length > 0, somatic_telomere_length > 0) |>
  mutate(log2_ratio = log2(somatic_telomere_length / germline_telomere_length))

virus_prevalence <- donors |>
  group_by(tumor_label) |>
  summarise(prevalence = 100 * mean(virus_category != NO_VIRUS), .groups = "drop") |>
  mutate(label = sprintf("%.1f", prevalence))

virus_composition <- donors |>
  add_count(tumor_label, name = "n_donors") |>
  filter(virus_category != NO_VIRUS) |>
  count(tumor_label, virus_category, n_donors, name = "n") |>
  mutate(
    percent = 100 * n / n_donors,
    virus_category = factor(virus_category, levels = VIRUS_LEVELS)
  )

fusion_donors <- filter(donors, !is.na(rna_fusion_count))

# DNA and RNA bars sit either side of each tumour location's position.
fusion_bins <- bind_rows(
  transmute(fusion_donors, tumor_label, mechanism = DNA_FUSIONS, count = dna_fusion_count),
  transmute(fusion_donors, tumor_label, mechanism = RNA_FUSIONS, count = rna_fusion_count)
) |>
  mutate(
    bin = cut(count, breaks = c(-Inf, 20, 50, Inf), labels = FUSION_BINS),
    mechanism = factor(mechanism, levels = names(MECHANISM_COLORS))
  ) |>
  count(tumor_label, mechanism, bin, name = "n") |>
  group_by(tumor_label, mechanism) |>
  mutate(percent = 100 * n / sum(n)) |>
  ungroup() |>
  mutate(x = as.numeric(tumor_label) + if_else(mechanism == DNA_FUSIONS, -0.2, 0.2))

fusion_totals <- fusion_donors |>
  count(tumor_label, name = "n") |>
  mutate(x = as.numeric(tumor_label), label = sprintf("(n=%d)", n))


# ---- shared styling ---------------------------------------------------------

base_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    plot.caption.position = "plot",
    plot.caption = element_text(
      face = "bold", size = 10, hjust = 0,
      margin = margin(t = 5, l = TITLE_OFFSET)
    ),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5)
  )

hide_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

scale_x_tumor <- function(...) {
  scale_x_discrete(limits = ordered_labels, drop = FALSE, ...)
}

scale_y_percent <- scale_y_continuous(
  limits = c(0, 100),
  breaks = seq(0, 100, 20),
  labels = function(z) paste0(z, "%"),
  expand = c(0, 0)
)

boxplot_with_medians <- function(medians, width) {
  list(
    geom_boxplot(
      width = width,
      fill = BOX_FILL,
      outlier.colour = BOX_OUTLIER,
      outlier.alpha = 0.5
    ),
    geom_point(
      data = medians,
      aes(y = median),
      shape = 18,
      size = 2.2,
      colour = MEDIAN_COLOR
    )
  )
}


# ---- panels -----------------------------------------------------------------

panel_a <- ggplot(t_cell, aes(tumor_label, t_cell_score)) +
  boxplot_with_medians(group_medians(t_cell, t_cell_score), width = 0.65) +
  scale_x_tumor(position = "top") +
  scale_y_continuous(
    transform = scales::transform_pseudo_log(sigma = 0.01, base = 10),
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
    axis.text.x.top = element_text(angle = LABEL_ANGLE, hjust = 0, vjust = 0, size = 7),
    plot.margin = margin(3, 8, 16, 14)
  )

panel_b <- ggplot(telomere, aes(tumor_label, log2_ratio)) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  boxplot_with_medians(group_medians(telomere, log2_ratio), width = 0.68) +
  scale_x_tumor() +
  coord_cartesian(ylim = c(-5, 5)) +
  labs(
    y = "log2(somatic /\ngermline telomere length)",
    caption = "B  Paired tumor–normal telomere-length change — diamonds indicate medians"
  ) +
  base_theme +
  hide_x +
  theme(plot.margin = margin(16, 8, 16, 14))

panel_c <- ggplot(virus_composition, aes(tumor_label, percent, fill = virus_category)) +
  geom_col(width = 0.68, show.legend = TRUE) +
  geom_text(
    data = virus_prevalence,
    aes(tumor_label, prevalence, label = label),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 2.25
  ) +
  scale_x_tumor() +
  scale_y_percent +
  coord_cartesian(clip = "off") +
  scale_fill_manual(name = "Virus category", values = VIRUS_COLORS, drop = FALSE) +
  labs(
    y = "Virus-positive\ndonors (%)",
    caption = "C  Virus presence — percentage of donors; values show prevalence"
  ) +
  base_theme +
  hide_x +
  theme(plot.margin = margin(16, 8, 18, 14))

# Continuous x so DNA and RNA bars can sit side by side; the expansion matches
# the discrete axis of panels A-C so the locations line up vertically.
panel_d <- ggplot(fusion_bins, aes(x, percent, fill = bin, colour = mechanism)) +
  geom_col(width = 0.34, linewidth = 0.3, show.legend = TRUE) +
  geom_text(
    data = fusion_totals,
    aes(x, 100, label = label),
    inherit.aes = FALSE,
    angle = LABEL_ANGLE,
    hjust = 0,
    vjust = 0,
    size = 2.05
  ) +
  scale_x_continuous(expand = expansion(add = 0.6)) +
  scale_y_percent +
  coord_cartesian(xlim = c(1, n_groups), clip = "off") +
  scale_fill_manual(name = "Fusions per donor", values = FUSION_COLORS, drop = FALSE) +
  scale_colour_manual(name = "Detection mechanism", values = MECHANISM_COLORS, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(fill = FUSION_COLORS[["21–50"]]))) +
  labs(
    y = "Donors with RNA data (%)",
    caption = "D  Donor-level DNA versus RNA fusion-count distributions"
  ) +
  base_theme +
  hide_x +
  theme(plot.margin = margin(18, 8, 3, 14))


# ---- figure -----------------------------------------------------------------

other_composition <- if (nrow(rare_locations) > 0L) {
  paste0(rare_locations$tumor_category, " (n=", rare_locations$n_donors, ")", collapse = "; ")
} else {
  "None"
}

figure_5 <- panel_a / panel_b / panel_c / panel_d +
  plot_layout(heights = c(1.25, 1.3, 1.2, 1.25), guides = "keep") +
  plot_annotation(
    title = "Figure 5. Immune, telomere, viral and RNA features by tumor location",
    caption = paste0("*Other composition: ", other_composition),
    theme = theme(
      plot.title = element_text(face = "bold", margin = margin(l = TITLE_OFFSET))
    )
  )

print(figure_5)
