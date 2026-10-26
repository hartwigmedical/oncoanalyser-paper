#!/usr/bin/env bash
# Downloads the large input resources (query VCFs, reference genome,
# truth VCF, region BEDs) into data/, since these aren't shipped in git
# (see README.md -- "Folder layout").
#
# Usage:
#   download_resources.sh                  # just the pipeline VCFs below
#   download_resources.sh <gs://bucket/prefix>   # also mirrors reference_data/
#                                                 # (and anything else) from a
#                                                 # bucket laid out like data/
#
# Part 1 (always runs): fetches SBX/Ilmn/Ultima's SAGE+PURPLE somatic VCFs
# and the DRAGEN downsampled VCF from gs://hartwig-oapaper-celllines (read
# access only) into data/<SEQ_Tech>/HG008T/ -- see README.md ("Folder
# layout") for the full <SEQ_Tech>/<sampleId>/... convention. The full-depth
# DRAGEN VCF, the downsampled SAGE/PURPLE sets, Ilmn_pcr/ and
# reference_data/ are NOT in that bucket and not covered by this part --
# populate them as you already do.
#
# Part 2 (only if a bucket prefix is given): the original generic mirror --
# `gsutil -m rsync -r <bucket> data/`, unchanged, for anything (e.g.
# reference_data/) that still lives in a bucket laid out exactly like data/
# locally, e.g. gs://bucket/prefix/reference_data/genome/...
#
# Requires: gsutil (part of the Google Cloud SDK), authenticated with read
# access (gcloud auth login, or application-default credentials already
# configured).
#
# Safe to re-run: `gsutil cp` overwrites in place; the Part 2 rsync only
# transfers new/changed files and deliberately does NOT pass -d (which
# would delete local files absent from the bucket) -- neither part ever
# removes anything you already have in data/.

set -euo pipefail

if ! command -v gsutil >/dev/null 2>&1; then
  echo "Error: gsutil not found. Install the Google Cloud SDK and run 'gcloud auth login' first." >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$ROOT_DIR/data"

mkdir -p "$DATA_DIR"

# Cell-line results bucket: <tech>/bwa/tn/oa3/ is the tumor/normal
# oncoanalyser run; sage/somatic/ and purple/ hold the two stages.
CELLLINE_RUNS="gs://hartwig-oapaper-celllines/HG008/wgs"

# seq_tech:bucket_tech_folder
# SBX: the bucket still holds a SAGE 5.0-beta.12 / PURPLE 4.4-beta.12 run, not
# the final sbx-r-f release (SAGE 5.0 / PURPLE 4.4) the manuscript tables
# use; it is being replaced. Check `##sageVersion=5.0` in the fetched
# HG008T.sage.somatic.vcf.gz header before running anything on SBX.
PIPELINE_SOURCES=(
  "SBX:roche_sbx"
  "Ilmn:illumina_novaseq6000"
  "Ultima:ultima_ug100"
)

# fetch <src> <dest>: copy one VCF plus its index (fetched if the bucket has
# one, otherwise built locally with bcftools)
fetch () {
  local src="$1" dest="$2"
  echo "  $src -> $dest"
  mkdir -p "$(dirname "$dest")"
  gsutil cp "$src" "$dest"
  if gsutil cp "${src}.tbi" "${dest}.tbi" 2>/dev/null; then
    :
  elif command -v bcftools >/dev/null 2>&1; then
    echo "    no remote .tbi -- indexing locally"
    bcftools index -t -f "$dest"
  else
    echo "    WARNING: no remote .tbi and bcftools not found -- $dest is unindexed" >&2
  fi
}

echo "[1/2] Fetching pipeline VCFs..."
for entry in "${PIPELINE_SOURCES[@]}"; do
  IFS=: read -r tech folder <<< "$entry"
  run="$CELLLINE_RUNS/$folder/bwa/tn/oa3"
  fetch "$run/sage/somatic/HG008T.sage.somatic.vcf.gz" "$DATA_DIR/$tech/HG008T/HG008T.sage.somatic.vcf.gz"
  fetch "$run/purple/HG008T.purple.somatic.vcf.gz" "$DATA_DIR/$tech/HG008T/HG008T.purple.somatic.vcf.gz"
done
# DRAGEN downsampled: DRAGEN 4.5.4 from the downsampled FASTQs (the bucket
# path doesn't say downsampled, but it is -- see README.md, "Notes on the
# downsampled VCFs")
fetch "$CELLLINE_RUNS/illumina_novaseq6000/dragen454/tn/linear/HG008R_HG008T.hard-filtered.vcf.gz" \
  "$DATA_DIR/Dragen_ds/HG008T/HG008T.hard-filtered.vcf.gz"

if [ $# -ge 1 ]; then
  BUCKET="${1%/}"  # strip any trailing slash for consistent concatenation/logging

  if [[ "$BUCKET" != gs://* ]]; then
    echo "Error: expected a gs:// URI, got: $BUCKET" >&2
    exit 1
  fi

  echo "[2/2] Syncing $BUCKET -> $DATA_DIR ..."
  gsutil -m rsync -r "$BUCKET" "$DATA_DIR"
else
  echo "[2/2] No bucket prefix given -- skipping the reference_data/misc mirror step."
fi

echo "Done. Contents of $DATA_DIR:"
find "$DATA_DIR" -maxdepth 3 -type d | sort
