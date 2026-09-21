# HG008-T smvar benchmark — self-contained handover package

Compares somatic small-variant (SNV/indel) VCFs from four
technology/pipeline combinations (DRAGEN, Illumina SAGE/PURPLE, SBX
SAGE/PURPLE, Ultima SAGE/PURPLE) against the NIST HG008-T draft somatic
SNV/INDEL benchmark (v0.3), using `aardvark`. Includes both the
original (full-depth) query VCFs and a downsampled set (~30x normal /
~90x tumor coverage) for Illumina, SBX, and Ultima — DRAGEN has no
downsampled VCF and is omitted from the downsampled outputs.

This package is self-contained and holds no intermediate results: run
the two scripts below and they build everything (region-tier BEDs,
per-run comparisons, manuscript tables, variant archives) from the
inputs already on disk. Everything is safe to re-run — both scripts
skip work that's already done and never overwrite existing output. A run
counts as "already done" only once it is marked complete, so an
interrupted or failed run is redone rather than reused (see "Resuming"
below).

## Folder layout

Query VCFs follow the manuscript-wide convention (2026-08-27
restructure): `data/<SEQ_Tech>/<sampleId>/<sampleId>.<stage>.somatic.vcf.gz`,
`sampleId` being the tumor sample name inside the VCF — matching how the
rest of the manuscript's (many more) samples are organized. `<SEQ_Tech>`
is one of `Dragen`/`Ilmn`/`SBX`/`Ultima`, plus a `<SEQ_Tech>_ds` variant
for downsampled data.

```
data/
  Dragen/HG008-T-mosaic/    HG008-T-mosaic.hard-filtered.vcf.gz(.tbi)
  Ilmn/HG008T/              HG008T.{sage,purple}.somatic.vcf.gz(.tbi)
  SBX/HG008T/               HG008T.{sage,purple}.somatic.vcf.gz(.tbi)     -- sbx-r-f
  Ultima/HG008T/            HG008T.{sage,purple}.somatic.vcf.gz(.tbi)
  Ilmn_ds/HG008T/           ~30x/90x downsampled Illumina
  SBX_ds/HG008T/            ~30x/90x downsampled SBX -- still sbx-r; no
                            downsampled sbx-r-f SAGE file exists yet, so
                            it wasn't swapped when the live/original set
                            was finalized to sbx-r-f (see
                            scripts/generate_variant_archive.py). DRAGEN
                            has no downsampled VCF, so there's no Dragen_ds/.
  Ultima_ds/HG008T/         ~30x/90x downsampled Ultima
  _legacy_sbx_versions/     superseded original-depth sbx-r files, outside
                            the SEQ_Tech convention
  reference_data/
    genome/           GRCh38 no-ALT reference (GCA_000001405.15), + .fai/.gzi
    truth/            NIST HG008-T v0.3 truth VCF (tumor-only, confirmed truncal somatic)
    regions/          Confidence-region BED tiers (see "Region tiers" below)
scripts/
  download_resources.sh          fetches data/ from a GCS bucket (see below)
  run_benchmark.sh               single pipeline x tier comparison (calls aardvark)
  append_ledger_rows.py          called by run_benchmark.sh; appends ledger rows
  generate_manuscript_table.sh   runs everything needed for the manuscript tables
  build_manuscript_tables.py     called by generate_manuscript_table.sh
  generate_variant_archive.py    builds the per-variant archive TSV
  tools/aardvark                 aardvark binary (Linux x86_64, run via Docker)
results/           NOT included -- created by the scripts (runs/, results.tsv,
                    manuscript TSVs, variant archive TSVs)
```

`data/` is intentionally the only top-level folder named that way — this
package is meant to live in a git repo whose `.gitignore` already
excludes any directory named `data` (or `figures`) anywhere in the
tree, so the large VCFs/reference genome under `data/` are never
accidentally committed; they're archived separately from git, in a
GCS bucket, and fetched with `download_resources.sh` (see "Quick
start" below). If you place this package somewhere without that
`.gitignore` convention already in place, add an equivalent ignore
rule for `data/` and `results/` yourself before running anything. A
local `.gitignore` in this folder also excludes `scripts/tools/aardvark`
(a third-party binary, not authored code — see "Prerequisites" below
for where to get it).

## Prerequisites

- `gsutil` (part of the Google Cloud SDK), authenticated with read
  access to the resource bucket — only needed for
  `download_resources.sh` (`gcloud auth login`, or
  application-default credentials already configured).
- `bcftools`, `tabix`, `bgzip` (htslib) — tested with 1.24
- `bedtools` — tested with v2.31.1
- `python3` — no third-party packages required (stdlib only)
- `docker`, with the daemon running — `aardvark` runs inside a
  `debian:bookworm-slim` container via `--platform linux/amd64`
  (Docker's own x86_64 emulation; no native macOS/ARM build of
  `aardvark` exists). Start Docker Desktop (or equivalent) before
  running either script.
- **`aardvark` binary** — not included in this repo (third-party tool,
  gitignored — see `.gitignore`). Download the Linux x86_64 build
  (tested with v0.10.5) from
  [github.com/PacificBiosciences/aardvark](https://github.com/PacificBiosciences/aardvark)
  releases, make it executable, and place it at
  `scripts/tools/aardvark`:
  ```bash
  mkdir -p scripts/tools
  # download the linux-x86_64 aardvark binary for the release you want, then:
  chmod +x scripts/tools/aardvark
  ```
- macOS users: `bash` must resolve to something in `PATH`; both
  scripts are written to run under the stock macOS bash (3.2), no
  Homebrew bash required. If `bcftools`/`tabix`/`bedtools` are
  installed via Homebrew and not yet on `PATH`, either add
  `/opt/homebrew/bin` to `PATH` yourself or note that
  `run_benchmark.sh` already does this internally.

## Quick start

```bash
cd scripts

# 0. Fetch data/ (query VCFs, reference genome, truth VCF, region BEDs)
#    from wherever it's archived. Safe to re-run -- only adds/updates,
#    never deletes.
./download_resources.sh gs://<bucket>/<prefix>

# 1. Manuscript tables (original + downsampled), ~35 aardvark runs total.
#    Takes a while (each run invokes Docker); safe to re-run/resume.
#    Must run before step 2 -- the archive imports its BD/EA columns
#    from these runs' aardvark output (specifically the *_widest runs).
./generate_manuscript_table.sh

# 2. Variant archive (original + downsampled). Runs no aardvark itself,
#    but depends on the *_widest runs from step 1 already existing.
./generate_variant_archive.py
./generate_variant_archive.py --downsampled
```

After both finish, `results/` contains:

```
results/
  runs/<run_id>/aardvark/          per-run detail: summary.tsv, truth.vcf.gz,
                                    query.vcf.gz (BD/EA/OA-tagged)
  runs/<run_id>/aardvark.log
  logs/<run_id>.log                run_benchmark.sh's own stdout/stderr per run
  results.tsv                      append-only ledger, one row per
                                    (run, comparison-mode, variant-type)
  Preliminary_results_wigits_manuscript.tsv              (+_SNV/_Insertion/_Deletion)
  Preliminary_results_wigits_manuscript_downsampled.tsv  (+_SNV/_Insertion/_Deletion)
  Variant_archive_truthset.tsv
  Variant_archive_truthset_downsampled.tsv
```

## Step-by-step: downloading resources

```bash
./download_resources.sh <gs://bucket/prefix>
```

Fetches `data/` from a GCS bucket via `gsutil -m rsync -r`. The bucket
is expected to already mirror the exact layout `data/` should have
locally — i.e. `<gs://bucket/prefix>/SBX/...` and
`<gs://bucket/prefix>/reference_data/...` map directly onto
`data/SBX/...` and `data/reference_data/...` (see "Folder layout"
above for the full `<SEQ_Tech>/<sampleId>/...` convention). Everything
under the given prefix is synced as-is; the script doesn't pick and
choose individual files.

Safe to re-run — `rsync` only transfers new/changed files, and the
script deliberately does not pass `-d` (which would delete local
files absent from the bucket), so it will never remove anything you
already have. If you're populating `data/` some other way (a shared
network mount, manual `gsutil cp`, etc.), this script is just a
convenience — nothing else in this repo depends on how `data/` gets
populated, only on its final layout matching "Folder layout" above.

## Step-by-step: manuscript tables

`generate_manuscript_table.sh` does 3 things:

1. **Builds any missing derived confidence-region BED tiers.** Only
   the primary/independent region tracks are shipped in
   `data/reference_data/regions/` (see "Region tiers" below); the 4 derived
   intersections (`all_hcr`, `all_giabhcr`, `nogermlinewithin50bp_hcr`,
   `nogermlinewithin50bp_giabhcr`) are built via `bedtools intersect`
   on first run and cached there afterward.
2. **Runs `run_benchmark.sh` for every pipeline x tier combination** —
   20 original (4 pipelines x 5 tiers) + 15 downsampled (3 pipelines x
   5 tiers, no DRAGEN). Each run is skipped if it is already marked
   complete — see "Resuming" below.
3. **Builds the manuscript TSVs from `results/results.tsv`**, one row
   per pipeline x tier, for `ALL` variants and each of `Snv` /
   `Insertion` / `Deletion` separately (8 files total: 4 original + 4
   downsampled).

### Table columns

`pipeline`, `bed_tier`, `truth_TP`, `query_TP`, `FP`, `FN`, `recall`,
`precision`, `f1` (the three metrics as percentages, or `NA` where the
ledger has no number for that pipeline x tier x variant type).

There are two TP columns because recall and precision are counted on
opposite sides of the comparison:

```
recall    = truth_TP / (truth_TP + FN)     # truth alleles aardvark matched
precision = query_TP / (query_TP + FP)     # query records aardvark matched
```

### Resuming

`run_benchmark.sh` writes `results/runs/<run_id>/.complete` as its very
last step, after the ledger rows are appended. 

* **marker present** → a finished run; `generate_manuscript_table.sh`
  skips it and reuses its ledger rows, and `run_benchmark.sh` refuses to
  overwrite it.
* **marker absent** → an interrupted or failed attempt (aardvark or the
  ledger step can fail after the directory is created).
  `run_benchmark.sh` discards the directory and redoes the run, so an
  interrupted sweep resumes by just re-invoking
  `generate_manuscript_table.sh`. A failed aardvark compare leaves its
  log behind at `results/logs/<run_id>.aardvark.failed.log`.

To force a *completed* run to redo, delete both
`results/runs/<run_id>/` and that run's rows in `results/results.tsv`;
`append_ledger_rows.py` refuses to append a second set of rows for a
`run_id` the ledger already has.

To run a single pipeline x tier combination directly (e.g. while
debugging), call `run_benchmark.sh` yourself:

```bash
./run_benchmark.sh <run_id> <tech_label> <query_vcf> <query_sample> [bed_tier] [rescue_filter]

# example: Illumina SAGE/PURPLE at the "Ultima HCR" tier, full PON/PONGnomad rescue
./run_benchmark.sh il_ultimahcr illumina_purple \
  ../data/Ilmn/HG008T/HG008T.purple.somatic.vcf.gz HG008T \
  all_hcr "PONGnomad,PON;PONGnomad,PON"
```

- `bed_tier` — one of `all`, `nogermlineinterference` (default),
  `nogermlinewithin50bp`, `all_hcr`, `all_giabhcr`,
  `nogermlinewithin50bp_hcr`, `nogermlinewithin50bp_giabhcr`.
- `rescue_filter` — comma-separated exact FILTER values to rewrite to
  `PASS` before scoring. The manuscript tables use
  `"PONGnomad,PON;PONGnomad,PON"` throughout (relaxes PURPLE's PON/gnomAD
  population filters on the query side, rather than masking those
  regions out of the comparison for every pipeline symmetrically).
- Query calls are filtered to `PASS` + tumor `FMT/AF >= 0.10` + reference
  contigs only, then multiallelic-split and MNV-atomized, before being
  compared against the truth VCF.
- Ledger rows are written with two `dosage_scoring` values:
  `raw` (aardvark's native GT/BASEPAIR numbers) and `ea_corrected`
  (reclassifies EA>=1 FPs — sites the truth side already expects an
  allele at, just with a different reported dosage/genotype — as TP,
  not FP; recall is unaffected since it depends only on the truth
  side). The manuscript tables use `ea_corrected` throughout.

## Step-by-step: variant archive

`generate_variant_archive.py` builds one row per truth-set allele,
plus one row per additional cross-technology "FP" candidate (FILTER
exactly `PASS`/`PON`/`PON;PONGnomad`/`PONGnomad` **and** tumour
`AF >= 0.10` at that pipeline's final stage.
REP_C diagnostic columns — no dependency on `results/results.tsv` for
those. It does, however, require `generate_manuscript_table.sh` to
have already run: the `<Pipeline>_BD`/`<Pipeline>_EA` columns (see
below) are imported directly from that pipeline's `*_widest` aardvark
run rather than re-derived from FILTER/AF, because aardvark does local
haplotype reconstruction that a simple FILTER/AF-based proxy cannot
replicate exactly (see the script's own docstring for a worked
example). Run step 1 in "Quick start" above first.

```bash
./generate_variant_archive.py               # original VCFs -> Variant_archive_truthset.tsv
./generate_variant_archive.py --downsampled # downsampled VCFs -> Variant_archive_truthset_downsampled.tsv
```

Columns: `CHROM/POS/REF/ALT/Variant_Type`, `in_truth_set`, region-membership
flags (`in_all_bed`, `in_GIAB_HCR`, `in_Ultima_HCR`,
`in_nogermlinewithin50bp`), `in_homopolymer` (true if any pipeline's
SAGE call carries a `REP_C` repeat-length annotation there), then per
pipeline: raw FILTER at each stage, tumor `AF`, `REP_C` (SAGE only,
homopolymer/repeat length at that position), and finally `<Pipeline>_BD`
/`<Pipeline>_EA` (aardvark's own TP/FN/FP verdict and dosage-match
indicator — see above; this is what any TP/FN/FP reconstruction from
this archive should use, not the FILTER/AF columns). The downsampled
archive has no `Dragen_*` columns (not just blank ones) since no
downsampled DRAGEN VCF exists.

## Region tiers

Shipped in `data/reference_data/regions/` (primary/independent inputs
only — derived tiers are built by the scripts, see above):

- `HG008-T_somatic_smvar_benchmark_v0.3_all.bed` — NIST's broadest
  confidence tier (full clonal truth set, no germline-interference
  exclusion, no HCR restriction). A position outside this bed is not
  scored by `aardvark` at all, regardless of any other tier.
- `HG008-T_somatic_smvar_benchmark_v0.3_nogermlineinterference.bed` —
  NIST's default/recommended tier, excludes regions where a somatic
  variant sits near/overlaps a germline variant in a way comparison
  tools can't reliably resolve.
- `HG008-T_somatic_smvar_benchmark_v0.3_nogermlinewithin50bp.bed` —
  NIST's strictest germline-interference exclusion.
- `ug_hcr_true.bed` — Ultima Genomics' published UG-HCR track (v3.1),
  independent of every pipeline compared here.
- `GIAB_hcr.bed` — GIAB's high-confidence-region track
  (`HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_..._v.3.3.2`),
  also independent of every pipeline compared here, but noted to be
  somewhat less technology-agnostic than `ug_hcr_true.bed` (retains
  some long-homopolymer territory Ultima's own HCR excludes).

The 5 tiers used throughout the manuscript tables:

| Tier key | BED tier | Built from |
|---|---|---|
| widest | `all` | shipped directly |
| giabhcr | `all_giabhcr` | `all` ∩ `GIAB_hcr.bed` |
| ultimahcr | `all_hcr` | `all` ∩ `ug_hcr_true.bed` |
| giabhcr_ni50 | `nogermlinewithin50bp_giabhcr` | `nogermlinewithin50bp` ∩ `GIAB_hcr.bed` |
| ultimahcr_ni50 | `nogermlinewithin50bp_hcr` | `nogermlinewithin50bp` ∩ `ug_hcr_true.bed` |

## Notes on the downsampled VCFs

The downsampled files (`*.ds.*`) target ~30x normal / ~90x tumor
coverage, closer to expected daily-practice sequencing depth
(cost-driven) than the original full-depth files. The *amount* of
actual depth reduction differs substantially by technology, since the
original files started at different depths (median normal/tumor DP):
Illumina 99x/121x, SBX 43x/119x, Ultima 35x/119x — so Illumina was cut
the most, Ultima barely at all. Bear this in mind before treating the
downsampled comparison as an apples-to-apples "all three technologies
at 30x/90x" result; it's better read as "cost of downsampling from
each technology's own original depth."
