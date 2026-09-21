#!/usr/bin/env bash
# Reusable HG008-T somatic SNV/INDEL benchmark runner.
#
# Usage:
#   run_benchmark.sh <run_id> <tech_label> <query_vcf> <query_sample> [bed_tier] [rescue_filter]
#
#   run_id        unique label for this run, e.g. sbx_purple_20260724. Must not be
#                 empty or contain '/'. A finished run is marked by
#                 results/runs/<run_id>/.complete (written as the very last step);
#                 a run directory WITHOUT that marker is an interrupted/failed run
#                 and is discarded and redone rather than reused -- see
#                 "Completion marker" below.
#   tech_label    free-text tech/pipeline tag, e.g. sbx, illumina, ultima
#   query_vcf     path to the pipeline's somatic VCF (bgzipped)
#   query_sample  sample name of the TUMOR sample column in query_vcf
#                 (must be explicit -- aardvark silently uses the first
#                 sample in the VCF if this is wrong, which for a
#                 normal-first tumor-normal VCF gives a false all-zero result)
#   bed_tier      all | nogermlineinterference (default) | nogermlinewithin50bp |
#                 all_hcr | all_giabhcr | nogermlinewithin50bp_hcr |
#                 nogermlinewithin50bp_giabhcr
#                 The last 4 are derived tiers (bedtools intersect of a primary
#                 tier with an independent HCR track) and are built
#                 automatically on first use if missing -- see
#                 "ensure_bed_tier" below and README.md.
#   rescue_filter optional, comma-separated list, default none. Any query
#                 record whose FILTER column is EXACTLY one of the listed
#                 values (whole-string match, e.g. "PONGnomad" or
#                 "PON;PONGnomad" -- each entry is one exact FILTER string,
#                 not a substring/contains check) is rewritten to PASS before
#                 the PASS filter below. Used to relax a specific PURPLE
#                 filter on the query side independently of region masking.
#                 No-op for pipelines that never emit that exact FILTER value
#                 (e.g. DRAGEN, which has no PONGnomad-equivalent filter at
#                 all). See README.md for the methodology this supports.
#
# Query calls with tumor FORMAT/AF < 10% are dropped before comparison. The
# truth benchmark is truncal-somatic-only by construction; low-VAF query
# calls fall outside that scope by design, not because they're necessarily
# wrong -- see README.md.
#
# Query calls are restricted to the reference's own contigs (no ALT contigs)
# before atomizing, since --atomize needs reference sequence for every
# record it touches and some pipelines call against ALT-contig-aware
# references.
#
# Query calls are also atomized (bcftools norm --atomize): MNV/block
# substitutions become consecutive single-base SNVs, matching how the truth
# VCF always represents adjacent phased substitutions (as separate SNVs,
# never as MNV blocks). See README.md ("MNV atomization").
#
# See README.md in the repo root for the full methodology.

set -euo pipefail

if [ $# -lt 4 ]; then
  echo "Usage: $0 <run_id> <tech_label> <query_vcf> <query_sample> [bed_tier] [rescue_filter]" >&2
  exit 1
fi

RUN_ID="$1"
TECH_LABEL="$2"
QUERY_VCF="$3"
QUERY_SAMPLE="$4"
BED_TIER="${5:-nogermlineinterference}"
RESCUE_FILTER="${6:-}"

# run_id becomes a directory name under results/runs/ that this script may delete
# (see "Completion marker" below), so reject anything that could escape it.
case "$RUN_ID" in
  ""|.|..|*/*) echo "Invalid run_id '$RUN_ID': must be non-empty and contain no '/'." >&2; exit 1 ;;
esac

export PATH="/opt/homebrew/bin:$PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REFDATA_DIR="$ROOT_DIR/data/reference_data"
TRUTH_VCF="$REFDATA_DIR/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz"
BED="$REFDATA_DIR/regions/HG008-T_somatic_smvar_benchmark_v0.3_${BED_TIER}.bed"
REF="$REFDATA_DIR/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz"
RESULTS_DIR="$ROOT_DIR/results"
RUN_DIR="$RESULTS_DIR/runs/$RUN_ID"
RESULTS_TSV="$RESULTS_DIR/results.tsv"

# Derived confidence-region tiers (bedtools intersect of a primary NIST tier
# with an independently-sourced HCR track) are not shipped -- they're
# regenerated here on first use so the region-BED footprint of this repo
# stays small. Safe to re-run: skipped if the file already exists.
ensure_bed_tier () {
  local tier="$1"
  local out="$REFDATA_DIR/regions/HG008-T_somatic_smvar_benchmark_v0.3_${tier}.bed"
  [ -f "$out" ] && return 0
  local primary hcr
  case "$tier" in
    all_hcr)                        primary="all";                    hcr="ug_hcr_true.bed" ;;
    all_giabhcr)                    primary="all";                    hcr="GIAB_hcr.bed" ;;
    nogermlinewithin50bp_hcr)       primary="nogermlinewithin50bp";    hcr="ug_hcr_true.bed" ;;
    nogermlinewithin50bp_giabhcr)   primary="nogermlinewithin50bp";    hcr="GIAB_hcr.bed" ;;
    *) return 0 ;;  # not a known derived tier; let the missing-file check below fail with a clear message
  esac
  local primary_bed="$REFDATA_DIR/regions/HG008-T_somatic_smvar_benchmark_v0.3_${primary}.bed"
  local hcr_bed="$REFDATA_DIR/regions/$hcr"
  echo "Building missing derived tier: ${tier}.bed = ${primary}.bed ∩ $(basename "$hcr_bed")" >&2
  bedtools intersect -a "$primary_bed" -b "$hcr_bed" > "$out"
}
ensure_bed_tier "$BED_TIER"

for f in "$TRUTH_VCF" "$BED" "$REF" "$QUERY_VCF"; do
  [ -f "$f" ] || { echo "Missing required file: $f" >&2; exit 1; }
done

# Completion marker.
#
# A run is "done" only once $RUN_DIR/.complete exists.
# A directory with the marker is a real result and is never touched. One without
# it is an interrupted/failed attempt, so it is discarded and the run redone.
if [ -d "$RUN_DIR" ]; then
  if [ -f "$RUN_DIR/.complete" ]; then
    echo "run_id '$RUN_ID' already completed at $RUN_DIR -- pick a unique run_id (or delete it first)." >&2
    exit 1
  fi
  echo "run_id '$RUN_ID' exists at $RUN_DIR but has no .complete marker (interrupted or failed run):" >&2
  echo "  discarding it and re-running. If that attempt already appended rows for this run_id to" >&2
  echo "  $RESULTS_TSV, append_ledger_rows.py will refuse to append again -- delete those rows first." >&2
  rm -rf "$RUN_DIR"
fi
mkdir -p "$RUN_DIR"

# Clean-up
cleanup_incomplete_run () {
  if [ -f "$RUN_DIR/.complete" ]; then
    return 0
  fi
  if [ -s "$RUN_DIR/aardvark.log" ]; then
    mkdir -p "$RESULTS_DIR/logs"
    cp "$RUN_DIR/aardvark.log" "$RESULTS_DIR/logs/$RUN_ID.aardvark.failed.log"
    echo "Kept aardvark's log at $RESULTS_DIR/logs/$RUN_ID.aardvark.failed.log" >&2
  fi
  rm -rf "$RUN_DIR"
}
trap cleanup_incomplete_run EXIT

if [ -n "$RESCUE_FILTER" ]; then
  echo "[1/3] Filtering query to PASS (rescuing FILTER in {$RESCUE_FILTER} exactly) + tumor AF>=10%, splitting multiallelics + atomizing MNVs, sorting..."
else
  echo "[1/3] Filtering query to PASS + tumor AF>=10%, splitting multiallelics + atomizing MNVs, sorting..."
fi
if ! bcftools query -l "$QUERY_VCF" | grep -qx "$QUERY_SAMPLE"; then
  echo "Sample '$QUERY_SAMPLE' not found in $QUERY_VCF. Available samples:" >&2
  bcftools query -l "$QUERY_VCF" >&2
  exit 1   # the EXIT trap removes the unfinished run directory
fi
REF_CONTIGS="$(cut -f1 "${REF}.fai" | paste -sd, -)"
bcftools view -s "$QUERY_SAMPLE" "$QUERY_VCF" \
  | awk -F'\t' -v OFS='\t' -v rescue="$RESCUE_FILTER" '
      BEGIN { n = split(rescue, rlist, ",") }
      /^#/ { print; next }
      {
        for (i = 1; i <= n; i++) { if ($7 == rlist[i]) { $7 = "PASS"; break } }
        print
      }
    ' \
  | bcftools view -f PASS \
  | bcftools view -i 'FMT/AF>=0.10' \
  | bcftools view -t "$REF_CONTIGS" \
  | bcftools norm -m -any --atomize -f "$REF" 2>/dev/null \
  | bcftools sort -Oz -o "$RUN_DIR/query.pass.vcf.gz" 2>/dev/null
tabix -f -p vcf "$RUN_DIR/query.pass.vcf.gz"

echo "[2/3] Running aardvark compare (linux/amd64 via Docker)..."
docker run --rm --platform linux/amd64 \
  -v "$REFDATA_DIR:/refdata" \
  -v "$SCRIPT_DIR/tools:/tools" \
  -v "$RUN_DIR:/out" \
  debian:bookworm-slim \
  /tools/aardvark compare \
    --reference "/refdata/genome/$(basename "$REF")" \
    --truth-vcf "/refdata/truth/$(basename "$TRUTH_VCF")" \
    --query-vcf /out/query.pass.vcf.gz \
    --query-sample "$QUERY_SAMPLE" \
    --regions "/refdata/regions/$(basename "$BED")" \
    --output-dir /out/aardvark \
    --compare-label "$RUN_ID" \
    --threads 4 \
  > "$RUN_DIR/aardvark.log" 2>&1 || { echo "aardvark failed -- its log is preserved, see the path printed below." >&2; exit 1; }

echo "[3/3] Appending raw + EA-aware dosage-corrected rows to results ledger..."
python3 "$SCRIPT_DIR/append_ledger_rows.py" \
  "$RUN_DIR" "$RESULTS_TSV" "$RUN_ID" "$TECH_LABEL" "$QUERY_VCF" "$QUERY_SAMPLE" "$BED_TIER"

# Last step, and only reached if everything above succeeded: mark the run finished.
date -u +%Y-%m-%dT%H:%M:%SZ > "$RUN_DIR/.complete"

echo "Done. Per-run detail: $RUN_DIR/aardvark/ (summary.tsv, truth.vcf.gz, query.vcf.gz w/ BD/EA/OA tags)"
echo "Ledger rows appended to: $RESULTS_TSV (dosage_scoring=raw is aardvark's native GT/BASEPAIR scoring;"
echo "dosage_scoring=ea_corrected treats EA>=1 FPs -- truth already expected an allele at that exact"
echo "position -- as not false positives; see README.md for why)."
