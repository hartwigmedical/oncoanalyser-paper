# ── Shared oncoplot utilities ─────────────────────────────────────────────────
# Used by: 4_drivers.ipynb, 5_actionable.ipynb
#
# Parameters:
#   entity_col  : column name for the entity axis ("gene" or "event")
#   status_col  : column used to test presence of an event
#   none_val    : value in status_col that means "no event"
#   pval_thresh : p-value cutoff for the direction arrow (default 0.1)

get_entity_order <- function(df, entity_col = "gene") {
  rev(df %>% select(all_of(c(entity_col, "rk"))) %>%
        drop_na(all_of(entity_col)) %>%
        unique() %>% arrange(rk) %>% pull(all_of(entity_col)))
}

get_sample_order <- function(data, entity_col = "gene",
                              status_col = "event", none_val = "none") {
  entities <- get_entity_order(data, entity_col)
  idx_col  <- paste0(entity_col, "_idx")
  data %>%
    distinct(sampleId, .data[[entity_col]], .keep_all = TRUE) %>%
    mutate(!!idx_col := if_else(.data[[status_col]] != none_val, 1L, 2L)) %>%
    select(sampleId, all_of(entity_col), all_of(idx_col)) %>%
    pivot_wider(names_from  = all_of(entity_col),
                values_from = all_of(idx_col), values_fill = 2L) %>%
    select(sampleId, all_of(rev(entities))) %>%
    arrange(across(-sampleId)) %>%
    pull(sampleId)
}

order_oncoplot <- function(data, entity_col = "gene",
                            status_col = "event", none_val = "none") {
  data[[entity_col]] <- factor(data[[entity_col]],
                                levels = get_entity_order(data, entity_col))
  data$sampleId      <- factor(data$sampleId,
                                levels = get_sample_order(data, entity_col,
                                                          status_col, none_val))
  data
}

get_fisher_highlights <- function(df, entity_col = "gene", pval_thresh = 0.1) {
  bind_rows(lapply(seq_len(nrow(df)), function(i) {
    mat <- matrix(c(df$events_PLT[i], df$events_REF[i],
                    df$non_events_PLT[i], df$non_events_REF[i]),
                  nrow = 2, byrow = FALSE)
    res <- fisher.test(mat)
    data.frame(entity = df[[entity_col]][i], pval = res$p.value, odds = res$estimate)
  })) %>%
  rename(!!entity_col := entity) %>%
  transmute(all_of(entity_col), field = "Direction", val0 = NA,
            val = case_when(pval <= pval_thresh & odds > 1 ~ "▲",
                            pval <= pval_thresh & odds < 1 ~ "▼",
                            TRUE ~ " "))
}

# get_label_pcts: compute per-entity percentage of non-none events across cohorts.
#   entity_col  : column that identifies the entity axis (e.g. "gene" or "event")
#   none_val    : value in status_col that means "no event"
#   status_col  : column whose value is compared to none_val (e.g. "event" or "top_level")
get_label_pcts <- function(df, entity_col, none_val, status_col) {
  df %>%
    group_by(cohort, .data[[entity_col]]) %>%
    summarise(pct = sum(.data[[status_col]] != none_val) / n(), .groups = "drop") %>%
    rename(field = cohort, val = pct) %>%
    mutate(val0 = val,
           val  = ifelse(100 * val <= 10,
                         as.character(round(100 * val, 1)),
                         as.character(round(100 * val, 0)))) %>%
    complete(field = c(plt_cohort, ref_cohort), !!sym(entity_col),
             fill = list(val = "0", val0 = 0))
}

# get_highlights: run Fisher test for each entity and return direction arrows.
#   entity_col  : column that identifies the entity axis (e.g. "gene" or "event")
#   none_val    : value in status_col that means "no event" (case-insensitive)
#   status_col  : column whose value is compared to none_val (e.g. "event" or "top_level")
get_highlights <- function(df, entity_col, none_val, status_col) {
  tmp <- df %>%
    mutate(cohort = ifelse(cohort == plt_cohort, "PLT", "REF")) %>%
    group_by(cohort, .data[[entity_col]]) %>%
    summarise(tot    = n(),
              events = sum(tolower(.data[[status_col]]) != tolower(none_val)),
              .groups = "drop") %>%
    drop_na()
  tmp2 <- tmp %>%
    transmute(cohort, !!sym(entity_col), events,
              non_events = as.double(tot - events)) %>%
    pivot_wider(names_from  = cohort,
                values_from = c(events, non_events), names_sep = "_",
                values_fill = 0)
  # Guard: REF or PLT columns may be absent when one cohort has no events at all
  if (!"events_REF"     %in% names(tmp2)) tmp2$events_REF     <- 0L
  if (!"non_events_REF" %in% names(tmp2)) tmp2$non_events_REF <- as.double(tmp2$non_events_PLT)
  if (!"events_PLT"     %in% names(tmp2)) tmp2$events_PLT     <- 0L
  if (!"non_events_PLT" %in% names(tmp2)) tmp2$non_events_PLT <- 0
  tmp2 <- tmp2 %>%
    mutate(across(starts_with("events_"), ~ replace_na(.x, 0)),
           non_events_REF = replace_na(non_events_REF,
                                       mean(.$non_events_PLT, na.rm = TRUE)))
  get_fisher_highlights(tmp2, entity_col)
}
