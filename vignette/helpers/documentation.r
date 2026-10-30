# ── Vignette bottom panel documentation strings ───────────────────────────────
# Edit these strings to update the annotations printed below the vignette.
# Sourced by 9_build.ipynb.

doc_header <- "Panel annotations and abbreviations"

# ── Tool glossary — 5 columns × 2 rows, richtext HTML for blue bold names ──────
# Tool names match subtitle colour (#304b99, bold). Used with annotate("richtext").
.tool <- function(name, url_suffix, desc) {
  paste0(
    tool_icon(size = 8), "&ensp;",
    "<b><span style='color:#304b99'>", name, "</span></b>",
    "  <span style='color:#aaaaaa'>github.com/hartwigmedical/hmftools/", url_suffix, "</span><br>",
    desc
  )
}

.db <- function(name, url, desc) {
  paste0(
    db_icon(size = 8), "&ensp;",
    "<b><span style='color:#304b99'>", name, "</span></b>",
    "  <span style='color:#aaaaaa'>", url, "</span><br>",
    desc
  )
}

doc_tools_c1 <- paste0(
  .tool("PURPLE", "purple",
        "Tumor purity, ploidy, & allele-specific copy number estimator."),
  "<br><br>",
  .tool("SAGE", "sage",
        "High-sensitivity somatic variant calling (SNV, INDEL, MNV) across the whole genome.")
)

doc_tools_c2 <- paste0(
  .tool("LINX", "linx",
        "SV interpretation — gene fusions, disruptions and ecDNA amplification detection."),
  "<br><br>",
  .tool("CHORD", "chord",
        "Homologous recombination deficiency (HRD) prediction from WGS mutational patterns.")
)

doc_tools_c3 <- paste0(
  .tool("CUPPA", "cuppa",
        "Cancer tissue-of-origin prediction from combined DNA and RNA genomic features."),
  "<br><br>",
  .tool("LILAC", "lilac",
        "HLA allele typing and somatic LOH detection, enabling immune escape analysis.")
)

doc_tools_c4 <- paste0(
  .tool("CIDER", "cider",
        "CDR3 assembly for T-cell and B-cell receptor diversity characterisation."),
  "<br><br>",
  .tool("TEAL", "teal",
        "Telomere length estimation; somatic/germline ratio reflects maintenance activity.")
)

doc_tools_c5 <- paste0(
  .tool("ISOFOX", "isofox",
        "RNA expression quantification and gene set activity scores for pathway analysis."),
  "<br><br>",
  .tool("VirusInterpreter", "virus-interpreter",
        "Viral genome annotation post-processed from VirusBreakend tool.")
)

# ── Row 2 — additional tools and databases ─────────────────────────────────────
doc_row2_c1 <- .tool("SIGS", "sigs",
  "Fits somatic SNVs to COSMIC reference signatures to quantify mutational process contributions.")

doc_row2_c2 <- .tool("PAVE", "pave",
  "Somatic & germline variant effect annotation.")

doc_row2_c3 <- .tool("COBALT", "cobalt",
  "Read-depth ratios for copy number analysis.")

doc_row2_c4 <- .db("Clinical", "hartwigmedicalfoundation.nl/en/data/database/",
  "Hartwig clinical database: location, biopsy, and treatments.")

doc_row2_c5 <- .db("CKB", "ckbhome.com",
  "Clinical knowledge database: Genomic alteration–therapy evidence.")

doc_left <- paste0(
  "Cohort Metadata: WGS - whole genome sequencing, RNASeq - RNA sequencing, ",
  "Chemo - Chemotherapy, Targeted - Targeted therapy, Hormonal - Hormonal therapy, ",
  "Immuno - Immunotherapy, Radio - Radiotherapy.\n",
  "Treatments - some samples received multiple treatment types. ",
  "Pre and post-biopsy treatment available. ",
  "Responders - complete or partial response per RECIST criteria.\n",
  "Tumor Characteristics: TMB - tumor mutational burden, MSI - microsatellite instability, ",
  "WGD - whole genome doubling, GIE - genetic immune escape, ",
  "HLA - human leukocyte antigen, LOH - loss of heterozygosity, ",
  "HRD - homologous recombination deficiency.\n",
  "Genomic Landscape: SNV - single nucleotide variant, INDEL - insertion/deletion, ",
  "MNV - multi-nucleotide variant, SV - structural variant ",
  "(deletions, duplications, inversions, translocations, complex clusters, LINE insertions).\n",
  "CDR3 - complementarity-determining region 3 ",
  "(T-cell and B-cell receptor diversity via CIDER). ",
  "Telomere ratio - somatic/germline telomere length ratio (TEAL).\n",
  "Transcriptomic Landscape (if available): Pathway activity scores from ISOFOX. ",
  "EMT - epithelial-mesenchymal transition, IFN-γ - interferon gamma, ",
  "TGF-β - transforming growth factor beta."
)

doc_middle <- paste0(
  "Processes Underlying Mutations: Estimated from single base substitutions using COSMIC v2. ",
  "Signatures grouped by aetiology.\n",
  "Copy Number Alteration Profile (PURPLE): ",
  "Inner ring - homozygous deletion (dark blue), ",
  "LOH with significant loss <0.6× ploidy (blue), near copy-neutral LOH (light blue).\n",
  "Outer ring - high amplification >3× ploidy (dark red), moderate >2× (red), ",
  "low >1.4× (light red). ",
  "Driver amplifications (red) and homozygous deletions (blue) annotated."
)

doc_oncoplots <- paste0(
  "Cancer Driver Landscape (PURPLE · LINX · VirusInterpreter): ",
  "TSG - tumor suppressor gene, ONC - oncogene, ecDNA - extrachromosomal DNA amplification. ",
  "Summary percentages based on all samples; oncoplot shows samples with ≥1 event.\n",
  "Pathway abbreviations: DDR - DNA damage repair, CC - cell cycle, ",
  "EPI - epigenetic regulation, ICK - immune checkpoint/HLA, ",
  "RAS - RAS/MAPK, PI3K - PI3K pathway, WNT - WNT signalling.\n",
  "Germline Predisposition (PURPLE · LINX): Germline driver events annotated by pathway.\n",
  "Potentially Actionable Events (Source CKB): ",
  "FDA - FDA approved, Off - drug repurposing, Exp - experimental, Oth - other.\n",
  "TMB high >", TMB_HIGH_THRESHOLD, " mut/Mb, MSI high >", MSI_HIGH_THRESHOLD, " indels/Mb."
)

doc_right <- paste0(
  "WGS vs Panel Coverage Summary: In-silico coverage study comparing WGS to comprehensive ",
  "and targeted panels across\nvarious genomic biomarkers. ",
  "WGS actionable biomarkers identified by the Jackson Clinical Knowledgebase (CKB) used as\n",
  "a baseline, and panel coverage determined by rules below:\n",
  "   -Signatures: TMB High and MSI High captured by comprehensive panel (no HRD); ",
  "only MSI High covered by targeted panel.\n",
  "   -Small variants/splice sites: Coverage measured by comparing genomic coordinates ",
  "of events to panel BED input files.\n",
  "   -Copy number: Events for panel genes with min copy number < 0.5 or ",
  "max copy number > 6 assumed covered.\n",
  "   -Other events: Fusions, HMZ (Homozygous) Disruptions, Viral inserts, ",
  "HLA type assumed not captured by panels.\n",
  "** See documentation for further details on the WGS vs Panel coverage study."
)
