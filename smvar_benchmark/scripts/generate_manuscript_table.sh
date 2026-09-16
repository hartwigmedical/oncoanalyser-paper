#!/usr/bin/env bash
# Regenerates the "Preliminary results" manuscript tables from scratch, for
# BOTH the original (full-depth) and downsampled (~30x normal / ~90x tumor)
# query VCFs.
#
# Original: 4 pipeline/technology combinations (DRAGEN, Illumina SAGE/PURPLE,
# SBX SAGE/PURPLE, Ultima SAGE/PURPLE).
# Downsampled: 3 combinations (Illumina/SBX/Ultima SAGE/PURPLE) -- DRAGEN has
# no downsampled VCF, so it's omitted from the downsampled tables entirely.
#
# Both run against the same 5 confidence-region tiers:
#   1. widest         -- all.bed (no germline-interference exclusion, no HCR)
#   2. giabhcr         -- all.bed ∩ GIAB_hcr.bed
#   3. ultimahcr        -- all.bed ∩ ug_hcr_true.bed
#   4. giabhcr_ni50    -- nogermlinewithin50bp.bed ∩ GIAB_hcr.bed
#   5. ultimahcr_ni50   -- nogermlinewithin50bp.bed ∩ ug_hcr_true.bed
#
# All runs use the finalized methodology: full PON/PONGnomad rescue
# (rescue_filter="PONGnomad,PON;PONGnomad,PON") -- see README.md
# ("Relaxing PONGnomad") for why this replaced symmetric region masking.
#
# Output (all under results/):
#   Preliminary_results_wigits_manuscript.tsv (+_SNV/_Insertion/_Deletion)
#     -- original, 20 rows (4 pipelines x 5 tiers), ALL-variant-type +
#        per-variant-type GT-mode numbers.
#   Preliminary_results_wigits_manuscript_downsampled.tsv (+_SNV/_Insertion/_Deletion)
#     -- downsampled, 15 rows (3 pipelines x 5 tiers).
#
# Safe to re-run: uses fixed run_ids (manuscript_<pipeline>_<tier> /
# manuscript_ds_<pipeline>_<tier>); if results/runs/<run_id>/.complete exists,
# that run is skipped and its cached ledger rows are reused. The marker -- not
# the directory -- is the test: a run directory left behind by a failed or
# interrupted attempt has no marker, and run_benchmark.sh discards it and redoes
# the run, so an interrupted sweep resumes by simply re-invoking this script.
# Delete results/runs/<run_id>/ and its results.tsv rows to force a completed
# run to re start.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REFDATA_DIR="$ROOT_DIR/data/reference_data"
VCF_DIR="$ROOT_DIR/data"
RESULTS_DIR="$ROOT_DIR/results"
RESULTS_TSV="$RESULTS_DIR/results.tsv"

export PATH="/opt/homebrew/bin:$PATH"

mkdir -p "$RESULTS_DIR/logs"

RESCUE="PONGnomad,PON;PONGnomad,PON"

echo "[0/2] Building any missing confidence-region BED tiers..."
build_tier () {
  local out="$REFDATA_DIR/regions/HG008-T_somatic_smvar_benchmark_v0.3_$1.bed"
  local a="$REFDATA_DIR/regions/HG008-T_somatic_smvar_benchmark_v0.3_$2.bed"
  local b="$3"
  if [ ! -f "$out" ]; then
    echo "  building $1.bed = $2.bed ∩ $(basename "$b")"
    bedtools intersect -a "$a" -b "$b" > "$out"
  fi
}
build_tier "all_giabhcr" "all" "$REFDATA_DIR/regions/GIAB_hcr.bed"
build_tier "all_hcr" "all" "$REFDATA_DIR/regions/ug_hcr_true.bed"
build_tier "nogermlinewithin50bp_giabhcr" "nogermlinewithin50bp" "$REFDATA_DIR/regions/GIAB_hcr.bed"
build_tier "nogermlinewithin50bp_hcr" "nogermlinewithin50bp" "$REFDATA_DIR/regions/ug_hcr_true.bed"

# pipeline_key:tech_label:query_vcf:query_sample
#
# Paths follow the manuscript-wide convention (2026-08-27 restructure):
# data/<SEQ_Tech>/<sampleId>/<sampleId>.<stage>.somatic.vcf.gz (and a
# <SEQ_Tech>_ds variant for downsampled data), sampleId being the tumor
# sample name in the VCF -- matches how the rest of the manuscript's
# (many more) samples are organized.
#
# SBX: the original (full-depth) set uses the finalized sbx-r-f release
# (see README.md -- "SBX pipeline finalized to sbx-r-f"). The downsampled
# set is intentionally still on the preceding sbx-r release -- no
# downsampled sbx-r-f SAGE file exists yet; update DOWNSAMPLED_PIPELINES
# once one is produced. Superseded original-depth sbx-r files live under
# data/_legacy_sbx_versions/, outside this convention.
ORIGINAL_PIPELINES=(
  "dragen:illumina_dragen:$VCF_DIR/Dragen/HG008-T-mosaic/HG008-T-mosaic.hard-filtered.vcf.gz:HG008-T-mosaic"
  "ilmnpurple:illumina_purple:$VCF_DIR/Ilmn/HG008T/HG008T.purple.somatic.vcf.gz:HG008T"
  "sbxrf:sbx_purple:$VCF_DIR/SBX/HG008T/HG008T.purple.somatic.vcf.gz:HG008T"
  "ultimapurple:ultima_purple:$VCF_DIR/Ultima/HG008T/HG008T.purple.somatic.vcf.gz:HG008T"
)
# no DRAGEN downsampled VCF exists -- omitted from the downsampled set
DOWNSAMPLED_PIPELINES=(
  "ilmnpurple:illumina_purple_ds:$VCF_DIR/Ilmn_ds/HG008T/HG008T.purple.somatic.vcf.gz:HG008T"
  "sbxr:sbx_purple_ds:$VCF_DIR/SBX_ds/HG008T/HG008T.purple.somatic.vcf.gz:HG008T"
  "ultimapurple:ultima_purple_ds:$VCF_DIR/Ultima_ds/HG008T/HG008T.purple.somatic.vcf.gz:HG008T"
)

# tier_key:bed_tier
TIERS=(
  "widest:all"
  "giabhcr:all_giabhcr"
  "ultimahcr:all_hcr"
  "giabhcr_ni50:nogermlinewithin50bp_giabhcr"
  "ultimahcr_ni50:nogermlinewithin50bp_hcr"
)

run_all () {
  local run_prefix="$1"; shift
  for pipeline in "$@"; do
    IFS=: read -r pkey tech vcf sample <<< "$pipeline"
    for tier in "${TIERS[@]}"; do
      IFS=: read -r tkey bed_tier <<< "$tier"
      run_id="${run_prefix}_${pkey}_${tkey}"
      RUN_IDS+=("$run_id")
      # Only a run marked .complete by run_benchmark.sh counts as cached. A bare
      # directory can be an aardvark/ledger failure that wrote no ledger rows, and
      # treating that as a success silently dropped a pipeline x tier from the
      # tables. run_benchmark.sh discards such a directory and redoes the run.
      if [ -f "$RESULTS_DIR/runs/$run_id/.complete" ]; then
        echo "  $run_id already completed, skipping run (reusing cached result)"
        continue
      fi
      echo "  running $run_id ($tech, tier=$bed_tier)"
      "$SCRIPT_DIR/run_benchmark.sh" "$run_id" "$tech" "$vcf" "$sample" "$bed_tier" "$RESCUE" \
        > "$RESULTS_DIR/logs/$run_id.log" 2>&1
    done
  done
}

echo "[1/2] Running all 20 original pipeline x tier combinations..."
RUN_IDS=()
run_all "manuscript" "${ORIGINAL_PIPELINES[@]}"
ORIGINAL_RUN_IDS=("${RUN_IDS[@]}")

echo "[1/2] Running all 15 downsampled pipeline x tier combinations..."
RUN_IDS=()
run_all "manuscript_ds" "${DOWNSAMPLED_PIPELINES[@]}"
DOWNSAMPLED_RUN_IDS=("${RUN_IDS[@]}")

echo "[2/2] Building manuscript tables from the ledger..."
python3 "$SCRIPT_DIR/build_manuscript_tables.py" \
  --results-tsv "$RESULTS_TSV" \
  --out-tsv "$RESULTS_DIR/Preliminary_results_wigits_manuscript.tsv" \
  --original-run-ids "${ORIGINAL_RUN_IDS[@]}" \
  --downsampled-run-ids "${DOWNSAMPLED_RUN_IDS[@]}"

echo "Done."
