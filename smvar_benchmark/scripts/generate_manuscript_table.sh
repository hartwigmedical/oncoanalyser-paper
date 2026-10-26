#!/usr/bin/env bash
# Regenerates the "Preliminary results" manuscript tables from scratch, for
# BOTH the original (full-depth) and downsampled (~30x normal / ~90x tumor)
# query VCFs.
#
# Original: 4 pipeline/technology combinations (DRAGEN, Illumina SAGE/PURPLE,
# SBX SAGE/PURPLE, Ultima SAGE/PURPLE).
# Downsampled: the same 4 combinations. DRAGEN's downsampled VCF comes from a
# newer DRAGEN version than its original -- see DOWNSAMPLED_PIPELINES below.
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
#     -- downsampled, 20 rows (4 pipelines x 5 tiers).
#   Supplementary_PCR_vs_PCRfree.tsv (+_SNV/_Insertion/_Deletion)
#     -- 10 rows: Illumina SAGE/PURPLE PCR-free (the manuscript_ilmnpurple_*
#        runs) and PCR-amplified (supp_pcr_*), 5 tiers each. Skipped with a
#        warning if the PCR VCF hasn't been placed (see README.md, "Folder
#        layout" -- it isn't fetched by download_resources.sh).
#   Supplementary_REPC_le8.tsv (+_SNV/_Insertion/_Deletion)
#     -- 20 rows: the 4 original pipelines on the 5 tiers with long tandem
#        repeats (REP_C > 8) removed from truth and query alike (supp_repcle8_*).
#        See README.md ("Region tiers").
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

# REP_C <= 8 supplementary tiers: each manuscript tier minus every reference
# tandem repeat SAGE would annotate REP_C > 8 (unit 1-5bp, >= 9 copies, +-1bp).
REPC_MASK="$REFDATA_DIR/regions/repc_gt8_mask.bed"
if [ ! -f "$REPC_MASK" ]; then
  echo "  building $(basename "$REPC_MASK") (~18 min, one pass over the reference)"
  python3 "$SCRIPT_DIR/build_repc_mask.py" \
    "$REFDATA_DIR/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz" "$REPC_MASK.raw"
  sort -k1,1 -k2,2n "$REPC_MASK.raw" | bedtools merge > "$REPC_MASK.tmp"
  mv "$REPC_MASK.tmp" "$REPC_MASK"
  rm "$REPC_MASK.raw"
fi
for base in all all_giabhcr all_hcr nogermlinewithin50bp_giabhcr nogermlinewithin50bp_hcr; do
  out="$REFDATA_DIR/regions/HG008-T_somatic_smvar_benchmark_v0.3_${base}_repcle8.bed"
  if [ ! -f "$out" ]; then
    echo "  building $(basename "$out") = ${base}.bed minus $(basename "$REPC_MASK")"
    bedtools subtract -a "$REFDATA_DIR/regions/HG008-T_somatic_smvar_benchmark_v0.3_${base}.bed" \
      -b "$REPC_MASK" > "$out"
  fi
done

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
# DRAGEN downsampled: DRAGEN 4.5.4 re-aligned from the downsampled FASTQs, vs.
# 4.2.4 from BAMs for the original -- see README.md ("Notes on the downsampled
# VCFs"). Fetched by download_resources.sh.
DOWNSAMPLED_PIPELINES=(
  "dragen:illumina_dragen_ds:$VCF_DIR/Dragen_ds/HG008T/HG008T.hard-filtered.vcf.gz:HG008T"
  "ilmnpurple:illumina_purple_ds:$VCF_DIR/Ilmn_ds/HG008T/HG008T.purple.somatic.vcf.gz:HG008T"
  "sbxr:sbx_purple_ds:$VCF_DIR/SBX_ds/HG008T/HG008T.purple.somatic.vcf.gz:HG008T"
  "ultimapurple:ultima_purple_ds:$VCF_DIR/Ultima_ds/HG008T/HG008T.purple.somatic.vcf.gz:HG008T"
)
# Supplementary: the same HG008 tumor/normal on Illumina with a PCR-amplified
# library (NovaSeq X), same oncoanalyser pipeline. Placed by hand, not fetched
# by download_resources.sh -- see README.md ("Folder layout").
PCR_PIPELINE="ilmnpcr:illumina_purple_pcr:$VCF_DIR/Ilmn_pcr/HG008T/HG008T.purple.somatic.vcf.gz:HG008T"

# tier_key:bed_tier
TIERS=(
  "widest:all"
  "giabhcr:all_giabhcr"
  "ultimahcr:all_hcr"
  "giabhcr_ni50:nogermlinewithin50bp_giabhcr"
  "ultimahcr_ni50:nogermlinewithin50bp_hcr"
)
REPC_TIERS=()
for t in "${TIERS[@]}"; do REPC_TIERS+=("${t}_repcle8"); done

# run_all <run_prefix> <pipeline>... -- runs every pipeline over RUN_TIERS,
# appending each run_id to RUN_IDS
run_all () {
  local run_prefix="$1"; shift
  for pipeline in "$@"; do
    IFS=: read -r pkey tech vcf sample <<< "$pipeline"
    for tier in "${RUN_TIERS[@]}"; do
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

RUN_TIERS=("${TIERS[@]}")

echo "[1/2] Running all 20 original pipeline x tier combinations..."
RUN_IDS=()
run_all "manuscript" "${ORIGINAL_PIPELINES[@]}"
ORIGINAL_RUN_IDS=("${RUN_IDS[@]}")

echo "[1/2] Running all 20 downsampled pipeline x tier combinations..."
RUN_IDS=()
run_all "manuscript_ds" "${DOWNSAMPLED_PIPELINES[@]}"
DOWNSAMPLED_RUN_IDS=("${RUN_IDS[@]}")

EXTRA_TABLES=()
IFS=: read -r _ _ pcr_vcf _ <<< "$PCR_PIPELINE"
if [ -f "$pcr_vcf" ]; then
  echo "[1/2] Running all 5 supplementary Illumina PCR tier combinations..."
  RUN_IDS=()
  run_all "supp_pcr" "$PCR_PIPELINE"
  # PCR-free rows reuse the original Illumina runs; rows ordered PCR-free then PCR
  PCRFREE_RUN_IDS=()
  for rid in "${ORIGINAL_RUN_IDS[@]}"; do
    case "$rid" in manuscript_ilmnpurple_*) PCRFREE_RUN_IDS+=("$rid") ;; esac
  done
  EXTRA_TABLES+=(--extra-table "$RESULTS_DIR/Supplementary_PCR_vs_PCRfree.tsv" \
    "${PCRFREE_RUN_IDS[@]}" "${RUN_IDS[@]}")
else
  echo "WARNING: Illumina PCR VCF not found at $pcr_vcf -- skipping Supplementary_PCR_vs_PCRfree" >&2
fi

echo "[1/2] Running all 20 supplementary REP_C<=8 pipeline x tier combinations..."
RUN_TIERS=("${REPC_TIERS[@]}")
RUN_IDS=()
run_all "supp_repcle8" "${ORIGINAL_PIPELINES[@]}"
EXTRA_TABLES+=(--extra-table "$RESULTS_DIR/Supplementary_REPC_le8.tsv" "${RUN_IDS[@]}")

echo "[2/2] Building manuscript tables from the ledger..."
python3 "$SCRIPT_DIR/build_manuscript_tables.py" \
  --results-tsv "$RESULTS_TSV" \
  --out-tsv "$RESULTS_DIR/Preliminary_results_wigits_manuscript.tsv" \
  --original-run-ids "${ORIGINAL_RUN_IDS[@]}" \
  --downsampled-run-ids "${DOWNSAMPLED_RUN_IDS[@]}" \
  "${EXTRA_TABLES[@]}"

echo "Done."
