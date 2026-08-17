#!/usr/bin/env bash
# Downloads the large input resources (query VCFs, reference genome,
# truth VCF, region BEDs) from a GCS bucket into data/, since these
# aren't shipped in git (see README.md -- "Folder layout").
#
# Usage:
#   download_resources.sh <gs://bucket/prefix>
#
# The bucket (at the given prefix) is expected to already mirror the
# exact layout data/ should have locally, e.g.:
#   gs://bucket/prefix/query_vcfs/original/...
#   gs://bucket/prefix/query_vcfs/downsampled/...
#   gs://bucket/prefix/reference_data/genome/...
#   gs://bucket/prefix/reference_data/truth/...
#   gs://bucket/prefix/reference_data/regions/...
# i.e. everything directly under the given prefix is copied as-is to
# the same relative path under data/.
#
# Requires: gsutil (part of the Google Cloud SDK), authenticated with
# read access to the bucket (gcloud auth login, or application-default
# credentials already configured).
#
# Safe to re-run: uses `gsutil -m rsync -r`, which only transfers
# new/changed files. Deliberately does NOT pass -d (which would delete
# local files absent from the bucket) -- this only ever adds/updates,
# never removes, so it won't touch anything you've placed in data/
# yourself (e.g. a manually-downloaded aardvark binary would live
# under scripts/tools/, not data/, so it's unaffected either way).

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <gs://bucket/prefix>" >&2
  exit 1
fi

BUCKET="${1%/}"  # strip any trailing slash for consistent concatenation/logging

if [[ "$BUCKET" != gs://* ]]; then
  echo "Error: expected a gs:// URI, got: $BUCKET" >&2
  exit 1
fi

if ! command -v gsutil >/dev/null 2>&1; then
  echo "Error: gsutil not found. Install the Google Cloud SDK and run 'gcloud auth login' first." >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$ROOT_DIR/data"

mkdir -p "$DATA_DIR"

echo "Syncing $BUCKET -> $DATA_DIR ..."
gsutil -m rsync -r "$BUCKET" "$DATA_DIR"

echo "Done. Contents of $DATA_DIR:"
find "$DATA_DIR" -maxdepth 3 -type d | sort
