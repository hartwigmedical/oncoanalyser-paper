library(dplyr)
library(ggplot2)

OTHER_GROUP <- "Other"
TITLE_OFFSET <- 45

hla_data <- result |>
  select(donorId, tumor_category, has_hla_loh) |>
  collect(max_bytes = 6 * 1024^2)

location_counts <- hla_data |>
  distinct(donorId, tumor_category) |>
  count(tumor_category, name = "n_donors") |>
  arrange(desc(n_donors), tumor_category)

rare_locations <- location_counts |>
  filter(n_donors < 10)

hla_data <- hla_data |>
  mutate(
    tumor_group = if_else(
      tumor_category %in% rare_locations$tumor_category,
      OTHER_GROUP,
      tumor_category
    )
  )

group_counts <- hla_data |>
  distinct(donorId, tumor_group) |>
  count(tumor_group, name = "n_donors")

ordered_groups <- group_counts |>
  filter(tumor_group != OTHER_GROUP) |>
  arrange(desc(n_donors), tumor_group) |>
  pull(tumor_group)

if (OTHER_GROUP %in% group_counts$tumor_group) {
  ordered_groups <- c(ordered_groups, OTHER_GROUP)
}

group_labels <- group_counts |>
  mutate(
    tumor_label = if_else(
      tumor_group == OTHER_GROUP,
      sprintf("Other* (n=%d)", n_donors),
      sprintf("%s (n=%d)", tumor_group, n_donors)
    )
  )

ordered_labels <- group_labels$tumor_label[
  match(ordered_groups, group_labels$tumor_group)
]

hla_summary <- hla_data |>
  group_by(tumor_group) |>
  summarise(
    n_donors = n_distinct(donorId),
    n_hla_loh = n_distinct(donorId[has_hla_loh == 1]),
    percent_hla_loh = 100 * n_hla_loh / n_donors,
    .groups = "drop"
  ) |>
  left_join(group_labels, by = "tumor_group") |>
  mutate(
    tumor_label = factor(tumor_label, levels = ordered_labels),
    label = sprintf("%.1f", percent_hla_loh)
  )

panel_e <- ggplot(
  hla_summary,
  aes(tumor_label, percent_hla_loh)
) +
  geom_col(width = 0.68, fill = "#76CFD8") +
  geom_text(
    aes(label = label),
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
  labs(
    x = NULL,
    y = "Donors with\nHLA LOH (%)",
    caption = "E  HLA allele loss of heterozygosity (gieHlaLoh)"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.caption.position = "plot",
    plot.caption = element_text(
      face = "bold",
      size = 10,
      hjust = 0,
      margin = margin(t = 5, l = TITLE_OFFSET)
    )
  )

print(panel_e)
