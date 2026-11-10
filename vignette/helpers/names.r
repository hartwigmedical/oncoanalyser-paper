# ── Tumor characteristic bar chart names ─────────────────────────────────────
namer <- function(i, map) {
  if (i %in% names(map)) map[[i]] else i
}

bar_names <- list(
  "purity"                 = "Purity",           # baseline
  "wholeGenomeDuplication" = "WGD",              # genome structure
  "lohGenomePct"           = "LOH",              # copy number (mean %)
  "meanClonality"          = "Clonal",           # tumour evolution
  "hrStatus"               = "HRD",             # DNA repair
  "svTmb"                  = "SV High",          # structural instability
  "svDoubleMinute"         = "ecDNA",
  "tmbStatus"              = "TMB High",         # immunotherapy biomarkers
  "msStatus"               = "MSI High",
  "gie"                    = "GIE",              # immune
  "gieHlaLoh"              = "HLA LOH",
  "cdr3Total"              = "Immune Infiltrated",
  "viralReported"          = "Viral reported"    # etiology
)

# ── Human-readable names for all marker columns ───────────────────────────────
# Usage:  marker_names[["tmbPerMb"]]              → "TMB (mut/Mb)"
#         df %>% rename(any_of(marker_names))     → renames matching cols
#         df %>% mutate(label = marker_names[col]) → in gathered data

marker_names <- c(
  # Mutational burden
  tmbPerMb                        = "TMB (mut/Mb)",
  tmbStatus                       = "TMB high",
  tml                             = "Tumor mutational load",
  msIndelsPerMb                   = "MSI (indels/Mb)",
  msStatus                        = "MSI high",
  indelsPerMb                     = "INDELs (per Mb)",
  snpPerMb                        = "SNVs (per Mb)",
  mnpPerMb                        = "MNVs (per Mb)",
  # Structural variants
  svTmb                           = "SV burden",
  svTotal                         = "SVs (total)",
  svDeletions                     = "SV deletions",
  svDuplications                  = "SV duplications",
  svInversions                    = "SV inversions",
  svTranslocations                = "SV translocations",
  svFoldBack                      = "SV fold-backs",
  svComplexClusters               = "SV complex clusters",
  svLine                          = "SV LINE insertions",
  svDoubleMinute                  = "SV double minutes",
  svSingleBreakend                = "SV single breakends",
  svInferred                      = "SV inferred",
  svLowVaf                        = "SV low VAF",
  # Genome structure
  purity                          = "Tumor purity",
  ploidy                          = "Mean ploidy",
  diploidProportion               = "Diploid proportion",
  wholeGenomeDuplication          = "Whole genome duplication",
  autosomesDoubled                = "Autosomes doubled",
  deletedGenes                    = "Deleted genes",
  copyNumberSegments              = "CN segments",
  unsupportedCopyNumberSegments   = "Unsupported CN segments",
  meanClonality                   = "Mean clonality",
  chrArmDeviations                = "Aneuploidy score",
  chrArmDeviationsPloidy          = "Aneuploidy score (ploidy adj.)",
  # Copy number
  cnaGenomePct                    = "CNA (%)",
  lohGenomePct                    = "LOH (%)",
  copyNeutralLohGenomePct         = "Copy-neutral LOH (%)",
  doublingGenomePct               = "Genome doubling (%)",
  biallelicDeletionsGenomePct     = "Biallelic deletions (%)",
  shallowDeletionsGenomePct       = "Shallow deletions (%)",
  highAmpGenomePct                = "High amplifications (%)",
  modAmpGenomePct                 = "Moderate amplifications (%)",
  # Signatures
  pctAging                        = "Aging (%)",
  pctApobec                       = "APOBEC (%)",
  pctChemicalExposure             = "Chemical exposure (%)",
  pctDnaDamageRepair              = "DNA damage repair (%)",
  pctTobacco                      = "Tobacco (%)",
  pctUvLight                      = "UV light (%)",
  allocationAging                 = "Aging (mut)",
  allocationApobec                = "APOBEC (mut)",
  allocationChemicalExposure      = "Chemical exposure (mut)",
  allocationDnaDamageRepair       = "DNA damage repair (mut)",
  allocationTobacco               = "Tobacco (mut)",
  allocationUvLight               = "UV light (mut)",
  # HRD
  hrd                             = "HRD score",
  hrStatus                        = "HRD",
  hrBRCA1                         = "HRD BRCA1",
  hrBRCA2                         = "HRD BRCA2",
  # Drivers
  driverTotal                     = "Total drivers",
  driverOnco                      = "Oncogene drivers",
  driverTsg                       = "TSG drivers",
  driverBiallelic                 = "Biallelic drivers",
  # Actionability
  protectFdaApproved              = "Actionable (FDA approved)",
  protectOffLabel                 = "Actionable (off-label)",
  protectExperimental             = "Actionable (experimental)",
  protectOther                    = "Actionable (other)",
  # Immune evasion
  gie                             = "Immune evasion",
  gieEvents                       = "Immune evasion events",
  gieHla                          = "HLA alleles",
  gieHlaLof                       = "HLA LOF",
  gieHlaLoh                       = "HLA LOH",
  gieApm                          = "Antigen processing (APM)",
  gieIfg                          = "IFN-γ pathway",
  giePdl1                         = "PD-L1 loss",
  gieCd58                         = "CD58 loss",
  gieSetdb1                       = "SETDB1 alteration",
  # Telomeres
  germlineTelomereLength          = "Telomere length (germline)",
  somaticTelomereLength           = "Telomere length (somatic)",
  somaticGermlineTelomereRatio    = "Telomere ratio (somatic/germline)",
  # Immune repertoire
  cdr3Total                       = "CDR3 sequences",
  cdr3Tcell                       = "CDR3 T-cell",
  cdr3Bcell                       = "CDR3 B-cell",
  cdr3TcellLocusTRB               = "CDR3 TRB",
  cdr3TcellLocusTRAD              = "CDR3 TRAD",
  cdr3TcellLocusTRG               = "CDR3 TRG",
  cdr3BcellLocusIGH               = "CDR3 IGH",
  cdr3BcellLocusIGK               = "CDR3 IGK",
  cdr3BcellLocusIGL               = "CDR3 IGL",
  # Viral
  viralReported                   = "Viral sequences (reported)",
  viralTotal                      = "Viral sequences (total)",
  # CUPPA
  cuppaCancertype_DNA_COMBINED    = "CUPPA type (DNA combined)",
  cuppaCancertype_GEN_POS         = "CUPPA type (genomic position)",
  cuppaCancertype_SNV96           = "CUPPA type (SBS96)",
  cuppaCancertype_EVENT           = "CUPPA type (events)",
  cuppaProb_DNA_COMBINED          = "CUPPA probability (DNA combined)",
  cuppaProb_GEN_POS               = "CUPPA probability (genomic position)",
  cuppaProb_SNV96                 = "CUPPA probability (SBS96)",
  cuppaProb_EVENT                 = "CUPPA probability (events)",
  # RNA signatures
  rnaAndrogen                     = "Androgen response",
  rnaEmt                          = "EMT",
  rnaEstrogenResponse             = "Estrogen response",
  rnaHrActivity                   = "HR activity",
  rnaHypoxia                      = "Hypoxia",
  rnaIfnGamma                     = "IFN-γ",
  rnaImmuneKilling                = "Immune killing",
  rnaMacrophageM1                 = "Macrophage M1",
  rnaMacrophageM2                 = "Macrophage M2",
  rnaMycTargets                   = "MYC targets",
  rnaNkCell                       = "NK cell",
  rnaP53Pathway                   = "p53 pathway",
  rnaPi3kAkt                      = "PI3K/AKT",
  rnaProliferation                = "Proliferation",
  rnaTgfBeta                      = "TGF-β",
  rnaWnt                          = "WNT"
)

mname <- function(col) if (col %in% names(marker_names)) marker_names[[col]] else col

marker_groups <- list(
  "Mutational burden"  = c("tmbPerMb","tmbStatus","tml","msIndelsPerMb","msStatus","indelsPerMb","snpPerMb","mnpPerMb"),
  "Structural variants"= c("svTmb","svTotal","svDeletions","svDuplications","svInversions","svTranslocations","svFoldBack","svComplexClusters","svLine","svDoubleMinute","svSingleBreakend"),
  "Genome structure"   = c("purity","ploidy","diploidProportion","wholeGenomeDuplication","autosomesDoubled","deletedGenes","copyNumberSegments","meanClonality","chrArmDeviations"),
  "Copy number"        = c("cnaGenomePct","lohGenomePct","copyNeutralLohGenomePct","doublingGenomePct","biallelicDeletionsGenomePct","shallowDeletionsGenomePct","highAmpGenomePct","modAmpGenomePct"),
  "Signatures"         = c("pctAging","pctApobec","pctChemicalExposure","pctDnaDamageRepair","pctTobacco","pctUvLight"),
  "HRD"                = c("hrd","hrStatus","hrBRCA1","hrBRCA2"),
  "Drivers"            = c("driverTotal","driverOnco","driverTsg","driverBiallelic"),
  "Actionability"      = c("protectFdaApproved","protectOffLabel","protectExperimental","protectOther"),
  "Immune evasion"     = c("gie","gieEvents","gieHla","gieHlaLof","gieHlaLoh","gieApm","gieIfg","giePdl1","gieCd58","gieSetdb1"),
  "Telomeres"          = c("germlineTelomereLength","somaticTelomereLength","somaticGermlineTelomereRatio"),
  "Immune repertoire"  = c("cdr3Total","cdr3Tcell","cdr3Bcell","cdr3TcellLocusTRB","cdr3TcellLocusTRAD","cdr3TcellLocusTRG","cdr3BcellLocusIGH","cdr3BcellLocusIGK","cdr3BcellLocusIGL"),
  "Viral"              = c("viralReported","viralTotal"),
  "CUPPA"              = c("cuppaCancertype_DNA_COMBINED","cuppaCancertype_GEN_POS","cuppaCancertype_SNV96","cuppaCancertype_EVENT","cuppaProb_DNA_COMBINED","cuppaProb_GEN_POS","cuppaProb_SNV96","cuppaProb_EVENT"),
  "RNA signatures"     = c("rnaAndrogen","rnaEmt","rnaEstrogenResponse","rnaHrActivity","rnaHypoxia","rnaIfnGamma","rnaImmuneKilling","rnaMacrophageM1","rnaMacrophageM2","rnaMycTargets","rnaNkCell","rnaP53Pathway","rnaPi3kAkt","rnaProliferation","rnaTgfBeta","rnaWnt")
)
    
# ── Primary tumor location → schema PNG filename ──────────────────────────────
schema_map <- c(
  "Breast"                               = "Breast",
  "Lung"                                 = "Lung",
  "Colon"                                = "Colon",
  "Prostate"                             = "Prostate",
  "Skin"                                 = "Skin",
  "Soft tissue"                          = "Soft_tissue",
  "Cancer of unknown primary (CUP)"      = "CUP",
  "Pancreas"                             = "Pancreas",
  "Ovary"                                = "Ovary",
  "Esophagus/gastroesophageal junction"  = "Esophagus",
  "Bladder"                              = "Bladder",
  "Kidney"                               = "Kidney",
  "Brain/central nervous system"         = "Brain_CNS",
  "Biliary tract"                        = "Biliary_gallbladder",
  "Rectum"                               = "Rectum",
  "Small bowel/intestine"                = "Small_intestine",
  "Endometrium/uterus"                   = "Endometrium",
  "Liver"                                = "Liver",
  "Stomach"                              = "Stomach",
  "Cervix"                               = "Cervix",
  "Ureter/renal pelvis"                  = "Urothelial",
  "Thyroid gland"                        = "Thyroid",
  "Lymphoid"                             = "Lymphoid",
  "Gallbladder"                          = "Gall_bladder",
  "Salivary gland"                       = "Head_and_Neck",
  "Anus"                                 = "Anus",
  "Pleura"                               = "Pleura",
  "Urethra"                              = "Urethra",
  "Vulva"                                = "Vulva",
  "Bone"                                 = "Bone",
  "Thymus"                               = "Thymus",
  "Oral cavity/tongue"                   = "Head_and_Neck",
  "Penis"                                = "Penis",
  "Nasal cavity"                         = "Head_and_Neck",
  "Oropharynx"                           = "Head_and_Neck",
  "Appendix"                             = "Appendix",
  "Adrenal gland"                        = "Adrenal_gland",
  "Nasopharynx"                          = "Head_and_Neck",
  "Vagina"                               = "Vulva_vagina",
  "Hypopharynx"                          = "Head_and_Neck",
  "Ampulla of Vater"                     = "Pancreaticobiliary",
  "Larynx"                               = "Head_and_Neck",
  "Paranasal sinuses"                    = "Head_and_Neck",
  "Eye"                                  = "All",
  "Testis"                               = "All",
  "Urachus"                              = "Urothelial",
  "Bone marrow/myeloid"                  = "All",
  "Peritoneum"                           = "All",
  "Unknown"                              = "All",
  "Other"                                = "All",
  "NULL"                                 = "All",
  "ALL"                                  = "All",
  # ── Batch file format (UPPERCASE_UNDERSCORE) ──────────────────────────────
  "ADRENAL_GLAND"                        = "Adrenal_gland",
  "LUNG"                                 = "Lung",
  "SOFT_TISSUE"                          = "Soft_tissue",
  "BOWEL"                                = "Colorectum",
  "BRAIN"                                = "Brain_CNS",
  "SKIN"                                 = "Skin",
  "BILIARY_TRACT"                        = "Biliary_gallbladder",
  "BLADDER"                              = "Bladder",
  "BONE"                                 = "Bone",
  "BREAST"                               = "Breast",
  "OVARY"                                = "Ovary",
  "CERVIX"                               = "Cervix",
  "STOMACH"                              = "Stomach",
  "UTERUS"                               = "Uterus",
  "HEAD_NECK"                            = "Head_and_Neck",
  "KIDNEY"                               = "Kidney",
  "LIVER"                                = "Liver",
  "LYMPH"                                = "Lymphoid",
  "PANCREAS"                             = "Pancreas",
  "PENIS"                                = "Penis",
  "PLEURA"                               = "Pleura",
  "PROSTATE"                             = "Prostate",
  "TESTIS"                               = "All",
  "THYMUS"                               = "Thymus",
  "THYROID"                              = "Thyroid",
  "VULVA"                                = "Vulva"
)

# ── Per-schema label positions ────────────────────────────────────────────────
# Override the default annotate() positions (x_left=0.16, y=2.1, x_right=0.82)
# when the inset or anatomy overlaps with a label.
# Format: list(x_left, y_left, x_right, y_right)
schema_label_positions <- list(
  "Biliary_gallbladder" = list(x_left = 0.05, y_left = 1.1, x_right = 0.82, y_right = 2.1)
)

get_schema_label_pos <- function(schema_file) {
  if (schema_file %in% names(schema_label_positions))
    schema_label_positions[[schema_file]]
  else
    list(x_left = 0.16, y_left = 2.1, x_right = 0.82, y_right = 2.1)
}

# ── Title → schema PNG override ──────────────────────────────────────────────
# Takes priority over schema_map when the cohort title matches exactly.
# Values are PNG filenames without extension (same convention as schema_map).
title_schema_map <- c(
  "Colorectal Adenocarcinoma"                                        = "Colorectum",
  "Colon Adenocarcinoma"                                             = "Colon",
  "Rectal Adenocarcinoma"                                            = "Rectum",
  "Gastrointestinal Neuroendocrine Tumors"                           = "Gastrointestinal",
  "Small Bowel Well-Differentiated Neuroendocrine Tumor"             = "Small_intestine",
  "Mucinous Adenocarcinoma of the Colon and Rectum"                  = "Colorectum",
  "Anal Squamous Cell Carcinoma"                                     = "Anal",
  "Appendiceal Adenocarcinoma"                                       = "Appendix",
  "Small Bowel Cancer"                                               = "Small_intestine",
  "High-Grade Neuroendocrine Carcinoma of the Colon and Rectum"      = "Colorectum",
  "Bowel"                                                            = "Colorectum"
)

# ── Vignette title → display title ───────────────────────────────────────────
# Location-only titles get "Cancer" appended; all others pass through unchanged.
# Usage: title_display(raw_title)
title_cancer_map <- c(
  "Adrenal Gland"        = "Adrenal Gland Cancer",
  "Biliary Tract"        = "Biliary Tract Cancer",
  "Bladder/Urinary Tract"= "Bladder/Urinary Tract Cancer",
  "Bone"                 = "Bone Cancer",
  "Bowel"                = "Bowel Cancer",
  "Breast"               = "Breast Cancer",
  "Cervix"               = "Cervix Cancer",
  "CNS/Brain"            = "CNS/Brain Cancer",
  "Esophagus/Stomach"    = "Esophagus/Stomach Cancer",
  "Head and Neck"        = "Head and Neck Cancer",
  "Kidney"               = "Kidney Cancer",
  "Liver"                = "Liver Cancer",
  "Lung"                 = "Lung Cancer",
  "Lymphoid"             = "Lymphoid Cancer",
  "Ovary/Fallopian Tube" = "Ovary/Fallopian Tube Cancer",
  "Pancreas"             = "Pancreatic Cancer",
  "Penis"                = "Penile Cancer",
  "Pleura"               = "Pleural Cancer",
  "Prostate"             = "Prostate Cancer",
  "Skin"                 = "Skin Cancer",
  "Soft Tissue"          = "Soft Tissue Cancer",
  "Testis"               = "Testicular Cancer",
  "Thymus"               = "Thymic Cancer",
  "Thyroid"              = "Thyroid Cancer",
  "Uterus"               = "Uterine Cancer",
  "Vulva/Vagina"         = "Vulvar/Vaginal Cancer"
)

title_display <- function(t) {
  if (t %in% names(title_cancer_map)) title_cancer_map[[t]] else t
}

cuppa_name_map <- c(
  "HPB: Bile duct/Gallbladder"                 = "CCA",
  "HPB: Pancreas"                              = "PDAC",
  "HPB: Liver"                                 = "HCC",
  "Esophagus/Stomach"                          = "ESCA/STAD",
  "Colorectum/Small intestine/Appendix"        = "CRC",
  "NET: Colorectum/Small intestine"            = "NET-GI",
  "NET: Pancreas"                              = "pNET",
  "NET: Lung"                                  = "lNET",
  "Lung: Non-small cell: LUAD"                 = "LUAD",
  "Lung: Non-small cell: LUSC"                 = "LUSC",
  "Lung: Small cell"                           = "SCLC",
  "Breast: Triple negative"                    = "TNBC",
  "Breast: Other"                              = "BRCA",
  "Gynecologic: Ovary/Fallopian tube"          = "OV",
  "Gynecologic: Endometrium"                   = "UCEC",
  "Anogenital"                                 = "CESC",
  "Prostate"                                   = "PRAD",
  "Urothelial tract"                           = "BLCA",
  "Kidney: Chromophobe"                        = "chRCC",
  "Kidney: Other"                              = "RCC",
  "Thyroid gland"                              = "THCA",
  "Head and neck: Salivary gland"              = "SGC",
  "Head and neck: Adenoid cystic"              = "ACC",
  "Head and neck: Other"                       = "HNSC",
  "Skin: Melanoma"                             = "SKCM",
  "Skin: Other"                                = "SKIN",
  "CNS: Glioma"                                = "GBM",
  "CNS: Medulloblastoma"                       = "MB",
  "CNS: Pilocytic astrocytoma"                 = "PA",
  "Bone/Soft tissue: GIST"                     = "GIST",
  "Bone/Soft tissue: Leiomyosarcoma"           = "LMS",
  "Bone/Soft tissue: Liposarcoma"              = "LPS",
  "Bone/Soft tissue: Osteosarcoma"             = "OS",
  "Bone/Soft tissue: Undiff. sarcoma"          = "UPS",
  "Bone/Soft tissue: Other"                    = "STS",
  "Myeloid: Acute myeloid leukemia"            = "AML",
  "Myeloid: Myeloproliferative neoplasm"       = "MPN",
  "Lymphoid tissue"                            = "LYM",
  "Mesothelium"                                = "MESO"
)

# ── Driver pathway abbreviations (≤4 chars) ───────────────────────────────────
pathway_map <- c(
  "DNA_damage_repair"      = "DDR",
  "Cell_cycle"             = "CC",
  "WNT"                    = "WNT",
  "Epigenetic_regulation"  = "EPI",
  "PI3K"                   = "PI3K",
  "RTK"                    = "RTK",
  "Apoptosis"              = "APO",
  "Immune_checkpoint_HLA"  = "ICK",
  "p53"                    = "p53",
  "Notch"                  = "NCH",
  "MYC"                    = "MYC",
  "RAS_MAPK"               = "RAS",
  "TGFB"                   = "TGFB",
  "NFKB"                   = "NFKB",
  "Hippo"                  = "HIP",
  "Hedgehog"               = "HH",
  "NRF2"                   = "NRF2",
  "Angiogenesis"           = "ANGI",
  "Telomere"               = "TEL",
  "Ferroptosis"            = "FER",
  "Autophagy"              = "ATG",
  "JAK_STAT"               = "JAK",
  "NULL"                   = "OTH",
  "Other"                  = "OTH"
)
