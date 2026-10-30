#!/usr/bin/env bash
# Run all vignette notebooks sequentially, producing up to six vignette layouts.
#
# ── Usage ─────────────────────────────────────────────────────────────────────
#
#   bash run.sh cohorts.csv                          # standard run
#   bash run.sh cohorts.csv --exclude-test-from-ref  # drop test samples from reference
#   bash run.sh                                      # generate cohorts via 0_cohort.ipynb
#   bash run.sh 5_actionable.ipynb                   # resume from a specific notebook
#
# ── Options ───────────────────────────────────────────────────────────────────
#
#   --exclude-test-from-ref   Remove any reference-cohort rows whose sampleId
#                             also appears in the test cohort. Useful when the
#                             reference is a pan-cancer set that includes the
#                             test samples.
#
# ── Cohorts CSV ───────────────────────────────────────────────────────────────
#
#   Required columns: sampleId, cohort, ref, primaryTumorLocation, title, doids
#   ref = TRUE for reference cohort, FALSE for test cohort
#   Valid primaryTumorLocation values are defined in helpers/names.r (schema_map).
#
# ── Outputs ───────────────────────────────────────────────────────────────────
#
#   {cohort}_vignette.png   9_build  — stacked bar CUPPA, boxplot RNA, documentation box
#
# Requires: jupyter with IRkernel installed

set -euo pipefail

NOTEBOOKS=(
  0_cohort.ipynb
  1_clinical.ipynb
  2_signatures.ipynb
  3_markers.ipynb
  4_drivers.ipynb
  5_actionable.ipynb
  6_circos.ipynb
  7_coverage.ipynb
  8_schemas.ipynb
  9_build.ipynb
)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Parse flags ───────────────────────────────────────────────────────────────
EXCLUDE_TEST_FROM_REF=false
MULTI_SIZE=false
POSITIONAL=()
for arg in "$@"; do
  case "$arg" in
    --exclude-test-from-ref) EXCLUDE_TEST_FROM_REF=true ;;
    --sizes) MULTI_SIZE=true ;;
    *) POSITIONAL+=("$arg") ;;
  esac
done

ARG="${POSITIONAL[0]:-}"
START=""
SKIP=false

# ── Validate a cohorts CSV ────────────────────────────────────────────────────
validate_cohorts_csv() {
  local file="$1"

  # Check required columns in header
  local header
  header=$(head -1 "$file")
  for col in sampleId cohort ref primaryTumorLocation title; do
    if ! echo "$header" | grep -q "$col"; then
      echo "Error: cohorts file is missing required column '$col'." >&2
      echo "       Expected columns: sampleId, cohort, ref, primaryTumorLocation, title" >&2
      exit 1
    fi
  done

  # Check at least one ref=TRUE and one ref=FALSE row
  local ref_col
  ref_col=$(echo "$header" | tr ',' '\n' | grep -n "^ref$" | cut -d: -f1)
  if ! awk -F',' -v col="$ref_col" 'NR>1 {if ($col=="TRUE")  found=1} END {exit !found}' "$file"; then
    echo "Error: cohorts file has no rows with ref=TRUE (reference cohort)." >&2
    exit 1
  fi
  if ! awk -F',' -v col="$ref_col" 'NR>1 {if ($col=="FALSE") found=1} END {exit !found}' "$file"; then
    echo "Error: cohorts file has no rows with ref=FALSE (test cohort)." >&2
    exit 1
  fi

  echo "✓ Cohorts file validated: $file"
}

# ── Determine run mode ────────────────────────────────────────────────────────
if [[ -z "$ARG" ]]; then
  # No argument: generate cohorts via 0_cohort.ipynb then run pipeline
  echo "ℹ No cohorts file provided — running 0_cohort.ipynb to generate cohorts."
  echo "  Tip: 'bash run.sh cohorts.csv' to use a pre-built config instead."
  SKIP=false

elif [[ "$ARG" == *.csv ]]; then
  # CSV cohorts file: validate, copy, skip 0_cohort
  if [[ ! -f "$ARG" ]]; then
    echo "Error: cohorts file '$ARG' not found." >&2
    exit 1
  fi
  validate_cohorts_csv "$ARG"
  REF_DATA=$(grep 'REF_DATA' helpers/map.r | grep -oP '(?<=")[^"]+(?=")' | head -1)
  READY_DATA="${REF_DATA}util/data/ready/"
  mkdir -p "$READY_DATA"
  cp "$ARG" "${READY_DATA}cohorts.txt"
  # Normalise ref column to logical TRUE/FALSE (fread reads it as character from external CSVs)
  Rscript -e "
    library(data.table)
    f <- '${READY_DATA}cohorts.txt'
    dt <- fread(f, fill = TRUE)
    # Drop leading row-number column if present (no header, all numeric)
    if (names(dt)[1] == 'V1' && all(grepl('^[0-9]+$', as.character(dt[[1]])))) dt[, V1 := NULL]
    dt[, ref := as.logical(ref)]
    fwrite(dt, f)
  "
  echo "✓ Cohorts file copied to ${READY_DATA}cohorts.txt"
  # Optionally remove test samples from the reference cohort
  if $EXCLUDE_TEST_FROM_REF; then
    Rscript -e "
      library(data.table)
      f <- '${READY_DATA}cohorts.txt'
      dt <- fread(f)
      test_ids <- dt[ref == FALSE, sampleId]
      before   <- nrow(dt[ref == TRUE])
      dt       <- dt[!(ref == TRUE & sampleId %in% test_ids)]
      after    <- nrow(dt[ref == TRUE])
      message('--exclude-test-from-ref: removed ', before - after,
              ' sample(s) from reference (', after, ' remaining).')
      fwrite(dt, f)
    "
  fi
  # Write run mode so 9_build.ipynb can suffix the output filename
  echo "$( $EXCLUDE_TEST_FROM_REF && echo 'exclusive' || echo 'overlap' )" > "${READY_DATA}run_mode.txt"
  # Write sizes flag so 9_build.ipynb knows whether to produce thumb/preview exports
  echo "$( $MULTI_SIZE && echo 'true' || echo 'false' )" > "${READY_DATA}multi_size.txt"
  START="1_clinical.ipynb"
  SKIP=true

else
  # Notebook name: resume from that notebook
  if ! printf '%s\n' "${NOTEBOOKS[@]}" | grep -qx "$ARG"; then
    echo "Error: '$ARG' is not a valid notebook name or .csv file." >&2
    echo "Valid notebooks: ${NOTEBOOKS[*]}" >&2
    exit 1
  fi
  START="$ARG"
  SKIP=true
fi

# ── Execute notebooks ─────────────────────────────────────────────────────────
for nb in "${NOTEBOOKS[@]}"; do
  if $SKIP; then
    [[ "$nb" == "$START" ]] && SKIP=false || continue
  fi
  echo "▶ Running $nb ..."
  jupyter nbconvert \
    --to notebook \
    --execute \
    --inplace \
    --ExecutePreprocessor.timeout=3600 \
    --ExecutePreprocessor.kernel_name=ir \
    --log-level WARN \
    "$nb"
  echo "✓ $nb done"
done

echo ""
echo "All notebooks completed successfully."

# ── Print output paths ────────────────────────────────────────────────────────
REF_DATA=$(grep 'REF_DATA' helpers/map.r | grep -oP '(?<=")[^"]+(?=")' | head -1)
SHARE_FIG="${REF_DATA}util/figures/share/"
VIGNETTE=$(ls -t "${SHARE_FIG}"*_vignette*.png 2>/dev/null | head -1 || true)
[[ -n "$VIGNETTE" ]] && echo "📄 Vignette: $VIGNETTE"
