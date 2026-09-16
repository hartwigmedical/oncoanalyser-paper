# Cross-platform cancer genomics pipeline
## oncoanalyser

Code and data used in the cross-platform publication of the
[oncoanalyser](https://nf-co.re/oncoanalyser/) pipeline, based on the
[WiGiTS](https://github.com/hartwigmedical/hmftools/#wigits) tools.

## Figure notebooks

Each manuscript figure is one self-contained notebook in the repository root. Extended Figures 1
and 3 are not generated here.


| Notebook | Contents |
| --- | --- |
| `figure_1_panel_d_redux_overview.ipynb` | Figure 1D — REDUX overview — unmapping, deduplication/consensus, base-quality recalibration, microsatellite jitter (GIAB HG008) |
| `figure_2.ipynb` | Figure 2 — REDUX read processing and tumour molecular profiles across Illumina, Roche and Ultima |
| `figure_3.ipynb` | Figure 3 — fresh-frozen vs FFPE (patients) and tumour/normal vs tumour-only (cell lines) |
| `figure_4.ipynb` | Figure 4 — WGS tumour/normal truth vs targeted panel |
| `extended_figure_2.ipynb` | Extended Figure 2 — per-region SNV & indel spectra and VAF across the three platforms |
| `extended_figure_4.ipynb` | Extended Figure 4 — fresh-frozen vs FFPE, sequencing QC and GIAB high-confidence benchmark |
| `extended_figure_5.ipynb` | Extended Figure 5 — why FFPE misses variants called in fresh-frozen |

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

### Reference data

Two notebooks need a reference file that is not part of either bucket:

* `extended_figure_4.ipynb` and `extended_figure_5.ipynb` download the **GIAB
  high-confidence region BED** (GRCh38, chr-prefixed, NISTv3.3.2 HG001) straight from the NIST
  GIAB FTP site over https — 10 MB, no credentials.
* `figure_4.ipynb` needs the targeted panel's capture regions. By default it reads them from
  the panel runs' own `bamtools/<sample>/<sample>.bam_metric.target_regions.tsv`, which covers
  exactly the same bases as the HMF oncopanel v3 capture BED; set `PANEL_BED` only if you have
  your own capture BED.

## Prerequisites

* **Python >= 3.11** (the figures were generated with 3.14.7).
* **[Google Cloud SDK](https://cloud.google.com/sdk/docs/install)** — only when reading inputs
  from `gs://`. The notebooks shell out to `gsutil`, and `extended_figure_5` additionally
  calls `gcloud auth print-access-token` to hand htslib a bearer token. Authenticate once with
  `gcloud auth login`. Not needed for local inputs.
* **`samtools` and `tabix` on `PATH`** — only for `extended_figure_5`, which reads BAM/VCF
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
