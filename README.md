# Cross-platform cancer genomics pipeline
## oncoanalyser

Code and data used in the cross-platform publication of the
[oncoanalyser](https://nf-co.re/oncoanalyser/) pipeline, based on the
[WiGiTS](https://github.com/hartwigmedical/hmftools/#wigits) tools.

## Figure notebooks

Each manuscript figure is one self-contained notebook in the repository root. Notebooks are named
for their analysis rather than their figure number: `fig_` for a main figure, `ext_` for an
extended/supplementary figure.


| Notebook | Contents |
| --- | --- |
| `fig_redux_overview.ipynb` | REDUX overview — unmapping, deduplication/consensus, base-quality recalibration, microsatellite jitter (GIAB HG008) |
| `fig_cross_platform.ipynb` | REDUX read processing and tumour molecular profiles across Illumina, Roche and Ultima |
| `fig_ffpe_and_tumor_only.ipynb` | Fresh-frozen vs FFPE (patients) and tumour/normal vs tumour-only (cell lines) |
| `fig_wgs_vs_panel.ipynb` | WGS tumour/normal truth vs targeted panel |
| `ext_cross_platform_region_spectra.ipynb` | Per-region SNV & indel spectra and VAF across the three platforms |
| `ext_ffpe_qc_benchmark.ipynb` | Fresh-frozen vs FFPE, sequencing QC and GIAB high-confidence benchmark |
| `ext_ffpe_false_negatives.ipynb` | Why FFPE misses variants called in fresh-frozen |
| `ext_tumor_only_benchmark.ipynb` | Tumour/normal vs tumour-only, benchmarked on Illumina and compared across the three platforms |
| `ext_hg008_pcr_vs_pcr_free.ipynb` | GIAB HG008 PCR-free vs PCR-amplified library prep — REDUX, variant burden and overlap |
| `ext_compar_concordance.ipynb` | COMPAR concordance of reportable output — FF tumour/normal vs tumour-only, NovaSeq 6000 vs X |

### Notebook → manuscript figure

Each notebook also states its figure number directly under its title.

**Main figures**

| Figure | Content | Source |
| --- | --- | --- |
| 1 | Schematics plus REDUX principle | Panel D: `fig_redux_overview.ipynb`; other panels not made here |
| 2 | Multiplatform performance on cell lines, incl. HG008 benchmark | `fig_cross_platform.ipynb` |
| 3 | Multimodal: FFPE, tumour-only | `fig_ffpe_and_tumor_only.ipynb` |
| 4 | Application: panel, exome | `fig_wgs_vs_panel.ipynb` |
| 5 | TCR, telomere, viral | Not made here |
| 6 | Longitudinal analysis | Not made here |
| 7 | Vignette | Not made here |

**Extended figures**

| Figure | Content | Source |
| --- | --- | --- |
| 1 | SAGE variant call visualisation vs IGV | oncoanalyser pipeline output (SAGE visualisation) |
| 2 | SNV and indel mutational profiles, shared and private variants per technology | `ext_cross_platform_region_spectra.ipynb` |
| 3 | PCR-free vs PCR, HG008 Illumina | `ext_hg008_pcr_vs_pcr_free.ipynb` |
| 4 | Circos plots HG008 Illumina, Roche, Ultima | oncoanalyser pipeline output (Circos plots) |
| 5 | FF vs FFPE differences in detail | `ext_ffpe_qc_benchmark.ipynb` |
| 6 | FFPE false negatives plus SAGE visualisation | `ext_ffpe_false_negatives.ipynb` |
| 7 | Tumour-only multiplatform | `ext_tumor_only_benchmark.ipynb` |
| 8 | COMPAR results, Illumina T/N vs tumour-only, reportable only, 5 cell lines | `ext_compar_concordance.ipynb` |

Every notebook opens with a **Data access** section that lists exactly which files it reads,
how large they are, and the three ways to supply them. Read that first.

## Input data

All inputs are public oncoanalyser output in two Google Cloud Storage buckets:

| Bucket | Content | Size |
| --- | --- | --- |
| `gs://hartwig-oapaper-celllines` | 5 cell lines x 3 sequencing platforms, WGS and targeted panel | 36 GB |
| `gs://hartwig-oapaper-patients` | 2 consented patients, fresh-frozen and FFPE WGS, targeted panel | — |

### Folder scheme

Every pipeline run has one folder, built from the same fields in the same order —
**sample, modality, platform, aligner, run type, pipeline**:

```text
gs://hartwig-oapaper-celllines/<sample>/<modality>/<platform>/<aligner>/<run type>/<pipeline>/
gs://hartwig-oapaper-patients/<sample>/<modality>/<preservation>/<platform>/<aligner>/<run type>/<pipeline>/

<sample>        COLO829 | HCC1187 | HCC1395 | HCC1954 | HG008      cell lines
                H00005238 | H00006429                              patients
<modality>      wgs | panel                 panel = HMF oncopanel v3, Illumina only
<preservation>  ff | ffpe                   patients only
<platform>      illumina_novaseq6000 | roche_sbx | ultima_ug100
<aligner>       bwa
<run type>      tn = tumour/normal | to = tumour-only
<pipeline>      oa3 = oncoanalyser 3.x
```

So COLO829, Illumina WGS, tumour/normal is:

```text
gs://hartwig-oapaper-celllines/COLO829/wgs/illumina_novaseq6000/bwa/tn/oa3/
```

Inside `oa3/` are the oncoanalyser output directories exactly as the pipeline writes them —
`alignments/`, `bamtools/`, `purple/`, `sage/`, `linx/`, `sigs/`, `logs/`, ... — with **no**
`results/<subject>` level in between.

Sample names inside a run are not uniform: tumours are `<sample>T`, normals are `<sample>R`
for COLO829, HG008 and the patients but `<sample>BL` for HCC1187/HCC1395/HCC1954, and the
Illumina COLO829 WGS run uses `COLO829v012T` / `COLO829v012R`. Nothing hard-codes them — every
notebook reads the real names out of `alignments/` and `purple/`.

Browse it with:

```bash
gsutil ls gs://hartwig-oapaper-celllines/
gsutil ls gs://hartwig-oapaper-celllines/COLO829/wgs/illumina_novaseq6000/bwa/tn/oa3/
```

### Three ways to supply the data

1. **Read straight from the bucket** (the default). Leave the `gs://` values in each
   Configuration cell alone. Each notebook stages only the files it needs into `data/` —
   between 11 MB and about 1 GB depending on the figure — and re-uses that cache afterwards.
2. **Copy the data down first.** Mirror the run folders a notebook needs with
   `gcloud storage rsync -r`, then set `CELLLINES_ROOT` / `PATIENTS_ROOT` in its Configuration
   cell to the local copy. Each notebook's Data access section has the exact commands.
3. **Run it on your own oncoanalyser output.** Point the Configuration cell at any
   oncoanalyser 3.x run directory. Sample names are read from the output, so they do not have
   to match ours.

Local paths must be **absolute** — `~` is not expanded, so a value like `~/oapaper-data` fails
with `FileNotFoundError`.

Inputs are staged into a local cache under `data/`, and figures are written to `figures/`. Both
directories are git-ignored and created on first run, so a fresh clone contains no data.

### COMPAR output

`ext_compar_concordance.ipynb` is the exception to the buckets: it reads
[COMPAR](https://github.com/hartwigmedical/hmftools/tree/master/compar) output, which is not
published. Run COMPAR on two sets of oncoanalyser run folders yourself (the notebook's Data access
section has the command) and set `COMPAR_OUTPUT` in its Configuration cell to each run's
reportable output folder, local or `gs://`.

Its root-cause panel is based on **manual curation** of this paper's discordant records, written
into the notebook itself (`MANUAL_ANNOTATIONS`). On your own data that panel shows every
record as *Not annotated*: you have to redo the curation and adapt the root-cause code to get a
similar figure.

### Reference data

Two notebooks need a reference file that is not part of either bucket:

* `ext_ffpe_qc_benchmark.ipynb` and `ext_ffpe_false_negatives.ipynb` download the **GIAB
  high-confidence region BED** (GRCh38, chr-prefixed, NISTv3.3.2 HG001) straight from the NIST
  GIAB FTP site over https — 10 MB, no credentials.
* `fig_wgs_vs_panel.ipynb` needs the targeted panel's capture regions. By default it reads them from
  the panel runs' own `bamtools/<sample>/<sample>.bam_metric.target_regions.tsv`, which covers
  exactly the same bases as the HMF oncopanel v3 capture BED; set `PANEL_BED` only if you have
  your own capture BED.

## Prerequisites

* **Python >= 3.11** (the figures were generated with 3.14.7).
* **[Google Cloud SDK](https://cloud.google.com/sdk/docs/install)** — only when reading inputs
  from `gs://`. The notebooks shell out to `gsutil`, and `ext_ffpe_false_negatives` additionally
  calls `gcloud auth print-access-token` to hand htslib a bearer token. Authenticate once with
  `gcloud auth login`. Not needed for local inputs.
* **`samtools` and `tabix` on `PATH`** — only for `ext_ffpe_false_negatives`, which reads BAM/VCF
  files through their indexes (remote or local). Skip these if you do not run that notebook.
* Optional: the figures use **Arial**; without it matplotlib falls back to Helvetica and then
  DejaVu Sans, which changes the typeface but nothing else.

## Set up the environment

Either a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

or conda/mamba, which also brings `samtools` and `tabix`:

```bash
conda env create -f environment.yml   # or: mamba env create -f environment.yml
conda activate oncoanalyser-paper
```

## Register the Jupyter kernel

The notebooks are saved against a kernel named `oncoanalyser-paper`. Register the environment
under that name so the notebooks open without having to pick a kernel:

```bash
python -m ipykernel install --user --name oncoanalyser-paper \
    --display-name "Python (oncoanalyser-paper)"
```

Then run the notebooks with `jupyter lab`, or open them in VS Code and select the
**Python (oncoanalyser-paper)** kernel.
