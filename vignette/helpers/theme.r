# Convert a grob to a ggplot via ggdraw — clip is off by default so annotations
# above y=1 (logo, website in the top margin) remain visible.
as_ggplot <- function(grob) {
  ggdraw() + draw_grob(grob, 0, 0, 1, 1)
}
go <- function(grob) as_ggplot(grob) + add_padding()

set_title <- function(hjust = 0.5, vjust = 0, size = 14, color = "#7F3D17", face = "bold") {
  theme(plot.title = element_text(hjust = hjust, vjust = vjust,
                                  size = size, colour = color, face = face))
}

add_padding  <- function(k = .2) theme(plot.margin = margin(k, k, k, k, "cm"))
add_legend   <- function(i = "none") theme(legend.position = i)

no_y_label   <- theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank())
no_x_label   <- theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank())
no_y_title   <- theme(axis.title.y = element_blank())
no_x_title   <- theme(axis.title.x = element_blank())
no_y         <- no_y_label + no_y_title
no_x         <- no_x_label + no_x_title
no_ticks     <- theme(axis.ticks   = element_blank())
no_y_line    <- theme(axis.line.y  = element_blank())
no_x_line    <- theme(axis.line.x  = element_blank())

bottom_legend    <- theme(legend.position = "bottom", legend.direction = "horizontal",
                          legend.text = element_text(size = 10))
top_legend       <- theme(legend.position = "top", legend.direction = "horizontal",
                          legend.text = element_text(size = 10))
no_legend        <- theme(legend.position = "none")
no_legend_title  <- theme(legend.title    = element_blank())
vertical_x       <- theme(axis.text.x     = element_text(angle = 90),
                           axis.title.x    = element_blank())

base_theme   <- theme_minimal(base_size = 13) + theme(panel.grid = element_blank())

vamos <- (
  theme(strip.placement      = "outside",
        strip.background     = element_blank(),
        legend.box.spacing   = unit(0, "pt"),
        legend.title         = element_blank())
  + add_padding()
  + add_legend()
  + set_title(color = "black")
)

panel_theme  <- vamos + set_title(hjust = 0, vjust = 1, size = 17) + add_padding()
common       <- base_theme + no_y_title + no_ticks + no_x_label + add_padding()
onco_theme   <- common + top_legend + no_y + no_legend_title
label_theme  <- common + no_legend + no_x_title + theme(axis.text.x = element_text(angle = 0))
title_set    <- set_title(size = 14, face = "plain", color = "black")
# tool_icon() / db_icon() return <img> tags for use inside subtitle strings.
# Requires ggtext::element_markdown() — blue_subtitle already sets this.
# Uses PNG (not SVG) — ggtext only supports raster images in <img> tags.
tool_icon <- function(size = 8) {
  sprintf("<img src='%s' width='%d' height='%d'>",
          file.path(REF_REPO, "ref/icons/tool.png"), size, size)
}
db_icon <- function(size = 10) {
  sprintf("<img src='%s' width='%d' height='%d'>",
          file.path(REF_REPO, "ref/icons/database.png"), size, size)
}

blue_subtitle <- theme(plot.subtitle = element_markdown(hjust = 0, size = 11,
                                                        color = "#304b99", face = "bold",
                                                        margin = margin(b = 4)))
