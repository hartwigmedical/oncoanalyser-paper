options(warn = -1)
shh <- suppressPackageStartupMessages

# User R library — extend path so packages installed outside conda are found
.libPaths(c("/home/jusset/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

shh(library(tidyverse))
shh(library(data.table))
shh(library(dtplyr))
shh(library(ggtext))
shh(library(cowplot))
shh(library(gridExtra))
shh(library(grid))
shh(library(png))

options(dplyr.summarise.inform = F)
options(repr.matrix.max.cols = 150, repr.matrix.max.rows = 200)

get_query_data <- function(query = "execute_sql_on_prod 'select * from hmfpatients.clinical'") {
    system_call <- system(query, intern = TRUE)
    as.data.frame(fread(text = system_call))
}

runQuery <- function(tableName, cols = "*", schema = "vignettes_hg38", where = NULL, sampleIds = NULL) {
    cols_str <- if (identical(cols, "*")) "*" else paste(cols, collapse = ", ")
    conditions <- c()
    if (!is.null(where))     conditions <- c(conditions, where)
    if (!is.null(sampleIds)) conditions <- c(conditions, paste0("sampleId in (\"", paste(sampleIds, collapse = "\",\""), "\")"))
    where_str <- if (length(conditions) == 0) "" else paste0(" where ", paste(conditions, collapse = " and "))
    get_query_data(paste0("execute_sql_on_prod 'select ", cols_str, " from ", schema, ".", tableName, where_str, "'"))
}

# Validate that a data frame contains required columns; stop with a clear message if not.
validate_columns <- function(df, required, context = deparse(substitute(df))) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0)
    stop(sprintf("[%s] Missing required columns: %s", context, paste(missing, collapse = ", ")))
  invisible(TRUE)
}

clean_filename <- function(x) {
    x <- gsub("[[:space:]]+", "_", x)       # replace spaces/tabs with underscores
    x <- gsub("[^a-zA-Z0-9_\\.-]", "", x)  # remove special characters except _, . and -
    x <- gsub("_+", "_", x)                # collapse multiple underscores
    x <- gsub("^_|_$", "", x)              # remove leading/trailing underscores
    return(x)
}

ggit <- function(png_file) {
  require(grid)
  require(ggplot2)
  img <- tryCatch({
    require(png)
    readPNG(png_file)
  }, error = function(e) {
    # Fallback: use magick to handle 16-bit or non-standard PNGs
    require(magick)
    tmp <- tempfile(fileext = ".png")
    on.exit(unlink(tmp))
    image_write(image_read(png_file), path = tmp, format = "png", depth = 8)
    readPNG(tmp)
  })
  g <- rasterGrob(img, interpolate = TRUE)
  ggplot() + theme_void() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}
