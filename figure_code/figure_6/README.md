# Figure 6 — Wisp ctDNA detection in a COLO829 titration (Illumina vs Ultima)

`fig_longitudinal.Rmd` regenerates Main Figure 6 from the public Wisp output in
`gs://hartwig-oapaper-celllines`. Knit it (RStudio *Knit* button, or
`rmarkdown::render("fig_longitudinal.Rmd")`).

## Experiment

A COLO829 titration series with expected tumor fractions of 0, 0.1, 1, 10 and 100% was
sequenced on two platforms (see the paper's Methods for sample preparation):

- **Illumina NovaSeq 6000**
- **Ultima Genomics UG 100**

Both series were analysed against one tissue reference: the Illumina whole-genome tumor/normal
run of COLO829 (`COLO829T` / `COLO829R`). The tumor SNVs Wisp tracks in plasma therefore come
from Illumina data for both platforms. For Ultima this is a cross-platform (Illumina-informed)
analysis.

## Tissue reference

Wisp is tumor-informed: it does not discover mutations in plasma, but looks up what it already
knows about the tumor from a tissue run, then measures that in each plasma sample.

The tissue reference is a standard oncoanalyser whole-genome tumor/normal run (`--mode wgts`)
on Illumina data of COLO829 (tumor `COLO829T`, normal `COLO829R`), stored at:

```text
gs://hartwig-oapaper-celllines/COLO829/ctdna/ultima_ug100/bwa/tissue_reference/oa3/COLO829/
```

The same tissue reference was used for both the Illumina and the Ultima plasma runs. What each
plasma run takes from it:

| From the tissue run | File | Used for |
|---|---|---|
| Tumor somatic variants | `purple/COLO829T.purple.somatic.vcf.gz` | The sites SAGE append measures in plasma, the basis of the SNV tumor fraction |
| Tumor purity, ploidy and copy number | `purple/` | Converting plasma VAFs into a tumor fraction (variant copy number per SNV), clonality correction, CNV method |
| Tumor LOH regions and B-allele frequencies | `purple/` copy number, `amber/COLO829T.amber.baf.tsv.gz` | The heterozygous sites that are lost in the tumor, which the LOH method tests in plasma |
| Normal sample alignment | `COLO829R` | The plasma AMBER and COBALT steps are run against the normal, as a germline baseline for the LOH and CNV methods |

In the Wisp summary, `TumorPurity` (0.99) and `TumorPloidy` (2.96) come from this tissue run, and
`TotalVariants` (49,036) is the number of tumor variants Wisp loaded from it. `CalcVariants` is
the number of SNVs left after Wisp's filters, which is what the SNV tumor fraction is based on.

## Sample mapping

The same titration point has a different sample id per platform. The figure's x-axis is the
expected tumor fraction:

| Expected tumor fraction (x-axis) | Illumina NovaSeq 6000 sample | Ultima UG 100 sample |
|---:|---|---|
| 0%   | `10ng0p-tr`   | `Colomimic0`   |
| 0.1% | `10ng01p-tr`  | `Colomimic01`  |
| 1%   | `10ng1p-tr`   | `Colomimic1`   |
| 10%  | `10ng10p-tr`  | `Colomimic10`  |
| 100% | `10ng100p-tr` | `Colomimic100` |

Naming: `<x>p` is the expected tumor percentage with the decimal point dropped
(`01` = 0.1%, `001` = 0.01%). Wisp prefixes its output with the patient id, e.g.
`COLO829_10ng1p-tr.wisp.summary.tsv`.

The 0.01% Illumina sample is in the bucket but not in the figure, because there is no Ultima
counterpart. The Rmd only reads the samples listed in its `SAMPLES` table.

## How the data were generated

| | |
|---|---|
| Pipeline | [nf-core/oncoanalyser](https://nf-co.re/oncoanalyser) v3.0.0 |
| Mode | `--mode purity_estimate --purity_estimate_mode wgts` (one run per platform) |
| Genome | GRCh38 (`GRCh38_hmf`), HMF reference data 3.0.0 |
| Tools | REDUX 2.0.5, SAGE 5.0.2 (append mode), AMBER 4.3, COBALT 3.0, **Wisp 1.3.1** |
| Wisp methods | `SOMATIC_VARIANT;AMBER_LOH;COPY_NUMBER` (SNV, LOH and CNV tumor fraction) |
| Plasma input, Illumina | FASTQ, aligned with BWA-MEM2 within oncoanalyser |
| Plasma input, Ultima | CRAM

Per plasma sample, oncoanalyser:
1. Runs REDUX (duplicate marking, BQR) on the plasma alignment.
2. Runs SAGE append to measure the tumor's somatic variants (`COLO829T.purple.somatic.vcf.gz` from
   the tissue reference) in the plasma sample.
3. Runs AMBER and COBALT on the plasma sample.
4. Runs Wisp with the tissue reference's PURPLE output (`purple/`: purity, ploidy, copy number), the plasma
   sample's SAGE append VCF, REDUX BQR, AMBER and COBALT output.

Tool versions and pipeline parameters are stored with each run in `software_versions.yml` and
`pipeline_info/`.

## Bucket layout

```text
gs://hartwig-oapaper-celllines/COLO829/ctdna/
├── illumina_novaseq6000/bwa/wisp/oa3/            # Illumina plasma run
│   ├── software_versions.yml
│   ├── pipeline_info/                              # params, execution report/trace
│   └── <sample>/                                   # 10ng0p-tr … 10ng100p-tr, 10ng001p-tr
│       ├── alignments/<sample>/                    # REDUX BQR, duplicate, jitter tables
│       ├── amber/
│       ├── cobalt/
│       ├── sage_append/<sample>/                   # tumor SNVs measured in plasma
│       └── wisp/                                   # Wisp output  ← used by the figure
└── ultima_ug100/bwa/
    ├── wisp/oa3/                                   # Ultima plasma run, same structure
    │   └── <sample>/                               # Colomimic0 … Colomimic100
    └── tissue_reference/oa3/COLO829/               # Illumina WGS tumor/normal tissue run,
        └── purple/  sage/  amber/  cobalt/  …      # used by BOTH plasma runs (see above)
```

The figure reads only `<sample>/wisp/COLO829_<sample>.wisp.summary.tsv`, 10 files of about 1.3 KB.
The other Wisp files (`somatic_variants`, `amber_loh`, `cn_segments`, `somatic_peak`,
`fragment_length`) contain per-variant and per-segment detail for further analysis.

## Key Wisp summary columns used

| Column | Meaning |
|---|---|
| `SNVPurity` | Reported SNV tumor fraction (0 when not significant) |
| `RawSNVPurity` | SNV tumor fraction before clonality correction |
| `SNVPurityLow` / `SNVPurityHigh` | 98% interval (Poisson 1% / 99% quantiles) |
| `SNVLod` | SNV limit of detection, **as a fraction** (×100 for %) |
| `SNV_MRD` | Detected: Poisson p < 0.01 against the background error rate |
| `CalcVariants` | Tumor SNVs used in the estimate |
| `WeightedAvgDepth` | Depth-weighted average plasma depth at those SNV sites |
| `LOHPurity`, `LOH_MRD` / `CopyNumberPurity`, `CNV_MRD` | LOH- and CNV-based tumor fraction and detection call |

See the [Wisp documentation](https://github.com/hartwigmedical/hmftools/tree/master/wisp) for details.

## Reproducing

Requirements:
- R with `rmarkdown`, `knitr`, `dplyr`, `tidyr`, `readr`, `stringr`, `ggplot2`, `patchwork`
- [Google Cloud SDK](https://cloud.google.com/sdk/docs/install) (`gsutil`) for reading the bucket
