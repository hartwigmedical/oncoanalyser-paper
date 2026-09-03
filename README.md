# Cross platfrom cancer genomics pipeline
## oncoanalyser
Code and Data used in the crossplatform publication of the [oncoanalyser](https://nf-co.re/oncoanalyser/) pipeline, based on the [WiGiTS](https://github.com/hartwigmedical/hmftools/#wigits) tools.

## Figure notebooks

Each manuscript figure is one self-contained notebook in the repository root:

| Notebook | Figure |
| --- | --- |
| `figure_1_panel_d_redux_overview.ipynb` | Figure 1, panel D |
| `figure_2.ipynb` | Figure 2 |
| `figure_3.ipynb` | Figure 3 |
| `figure_4.ipynb` | Figure 4 |
| `supplementary_figure_1.ipynb` | Supplementary Figure 1 |
| `supplementary_figure_2.ipynb` | Supplementary Figure 2 |
| `supplementary_figure_3.ipynb` | Supplementary Figure 3 |

### Input data

Every notebook opens with a configuration cell that points at one or more oncoanalyser
run directories. Each of those may be **either a local path or a `gs://` URI** — the IO
helpers branch on the `gs://` prefix and otherwise fall back to plain `pathlib` and
`shutil`, so listing, reading and staging behave identically for both. The values
committed here are the `gs://` locations the published figures were produced from;
repoint them at a local run directory and no cloud access is needed.

Inputs are staged into a local cache under `data/`, and figures are written to
`figures/`. Both directories are git-ignored and created on first run, so a fresh clone
contains no data.

When configuring a local path, give an absolute one: `~` is not expanded, so a value
like `~/HG008` fails with `FileNotFoundError`.

### Prerequisites

* **Python >= 3.11** (the figures were generated with 3.14.7).
* **[Google Cloud SDK](https://cloud.google.com/sdk/docs/install)** — only when reading
  inputs from `gs://`. The notebooks shell out to `gsutil`, and `supplementary_figure_3`
  additionally calls `gcloud auth print-access-token` to hand htslib a bearer token.
  Authenticate once with `gcloud auth login`. Not needed for local inputs.
* **`samtools` and `tabix` on `PATH`** — only for `supplementary_figure_3`, which reads
  BAM/VCF files through their indexes (remote or local). Skip these if you do not run
  that notebook.
* Optional: the figures use **Arial**; without it matplotlib falls back to Helvetica and
  then DejaVu Sans, which changes the typeface but nothing else.

### Set up the environment

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

### Register the Jupyter kernel

The notebooks are saved against a kernel named `oncoanalyser-paper`. Register the
environment under that name so the notebooks open without having to pick a kernel:

```bash
python -m ipykernel install --user --name oncoanalyser-paper \
    --display-name "Python (oncoanalyser-paper)"
```

Then run the notebooks with `jupyter lab`, or open them in VS Code and select the
**Python (oncoanalyser-paper)** kernel.
