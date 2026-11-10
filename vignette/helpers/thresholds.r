# ── Clinical thresholds ────────────────────────────────────────────────────────
# Centralised constants used across notebooks and figure functions.
# Edit here to adjust thresholds globally.

# Mutational burden
TMB_HIGH_THRESHOLD  <- 10   # mut/Mb  — PURPLE/SAGE
MSI_HIGH_THRESHOLD  <- 4    # indels/Mb — PURPLE/SAGE

# Genome structure
PURITY_LOW_THRESHOLD <- 0.20  # samples below this are flagged as low purity

# ── Layout constants ───────────────────────────────────────────────────────────
VIGNETTE_WIDTH  <- 42   # cm
VIGNETTE_HEIGHT <- 26   # cm
VIGNETTE_DPI    <- 300

# ── CUPPA ─────────────────────────────────────────────────────────────────────
CUPPA_PROB_THRESHOLD          <- 0.5   # default min probability
CUPPA_DNA_COMBINED_THRESHOLD  <- 0.8   # stricter threshold for DNA Combined classifier
CUPPA_CLASSIFIER_THRESHOLD    <- 0.2   # threshold for GEN_POS, SNV96, EVENT classifiers

# ── Immune / SV thresholds ────────────────────────────────────────────────────
CDR3_INFILTRATED_THRESHOLD <- 15
TCELL_HIGH_THRESHOLD       <- 20
SV_HIGH_THRESHOLD          <- 300

# ── Oncoplot display limits ───────────────────────────────────────────────────
ONCO_DISPLAY_TOP_SOMATIC   <- 30
ONCO_DISPLAY_TOP_GERMLINE  <- 10

# ── Build layout ─────────────────────────────────────────────────────────────
BOTTOM_STRIP_LINES         <- 16
