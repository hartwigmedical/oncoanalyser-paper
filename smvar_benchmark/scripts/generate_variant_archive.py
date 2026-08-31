#!/usr/bin/env python3
"""Archive every truth-set variant, PLUS every candidate "FP" variant
(FILTER exactly PASS/PON/PON;PONGnomad/PONGnomad at the final pipeline
stage -- i.e. what the full-rescue methodology would count as PASS --
in ANY of the included techs, but absent from the truth set), with
region membership and each pipeline's raw call status, into a single
convenient TSV.

Depends on `results/runs/manuscript_<pipeline>_widest/` (or
`manuscript_ds_<pipeline>_widest/` for --downsampled) already existing
-- run generate_manuscript_table.sh first. See "<Pipeline>_BD"/
"<Pipeline>_EA" columns below.

Usage:
  ./generate_variant_archive.py                 # original (full-depth) VCFs, 4 pipelines incl. DRAGEN
  ./generate_variant_archive.py --downsampled    # ~30x/90x downsampled VCFs, 3 pipelines (no DRAGEN)

Output (under results/):
  Variant_archive_truthset.tsv               (--downsampled omitted)
  Variant_archive_truthset_downsampled.tsv   (--downsampled passed)
Each: one row per truth-set allele, plus one row per additional
cross-tech "FP" candidate not in the truth set. An `in_truth_set`
column (1/0) distinguishes the two kinds of row.

Per-pipeline `<Pipeline>_BD`/`<Pipeline>_EA` columns carry aardvark's
own verdict (TP/FN for truth rows, FP/TP + dosage EA for FP-only rows),
imported directly from that pipeline's `manuscript[_ds]_<pipeline>_widest`
aardvark run rather than re-derived from FILTER/AF here. This is
necessary, not cosmetic: aardvark does local haplotype reconstruction,
matching truth alleles against a differently-represented query allele
nearby when they describe the same net sequence change (e.g. truth's
`CA>C` + `A>T` + `C>T` atomic calls vs. a caller's own `TCAA>T` +
`C>CTT` block substitution covering the same bases) -- a case no
exact-`(chrom,pos,ref,alt)`-key lookup can ever see, however the
FILTER/AF-based "does this look PASS-quality" proxy is computed. Using
these imported columns (combined with the region-membership columns,
which use full reference-span containment -- see build_bed_keyset)
reproduces the manuscript tables' TP/FN/FP counts exactly, rather than
approximately; the FILTER/AF/REP_C columns remain for diagnostic use
(homopolymer context, raw AF, etc.), just no longer for TP/FN/FP
determination. The widest tier is used as the single source for all 5
tiers because aardvark's BD/EA decision at a given key is tier-
invariant (region restriction only gates what's in scope, not how a
variant is matched) -- verified by diffing widest-tier vs. Ultima-HCR-
tier tags for every shared key: 0 mismatches.
"""
import argparse
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
REFDATA_DIR = ROOT_DIR / "data" / "reference_data"
VCF_DIR = ROOT_DIR / "data"
RESULTS_DIR = ROOT_DIR / "results"

REF = REFDATA_DIR / "genome" / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz"
TRUTH_VCF = REFDATA_DIR / "truth" / "HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz"

BEDS = {
    # NIST's broadest confidence-region tier -- the base region every tier
    # used in this report is built from (directly or via further
    # restriction). A position outside this isn't scored by aardvark at
    # all in any run in this report, regardless of what the other three
    # (independent, non-NIST-truth-derived) tracks below say.
    "in_all_bed": REFDATA_DIR / "regions" / "HG008-T_somatic_smvar_benchmark_v0.3_all.bed",
    "in_GIAB_HCR": REFDATA_DIR / "regions" / "GIAB_hcr.bed",
    "in_Ultima_HCR": REFDATA_DIR / "regions" / "ug_hcr_true.bed",
    "in_nogermlinewithin50bp": REFDATA_DIR / "regions" / "HG008-T_somatic_smvar_benchmark_v0.3_nogermlinewithin50bp.bed",
}

# pipeline -> {stage_label: (vcf_path, tumor_sample)}
# Paths follow the manuscript-wide convention (2026-08-27 restructure):
# data/<SEQ_Tech>/<sampleId>/<sampleId>.<stage>.somatic.vcf.gz (and a
# <SEQ_Tech>_ds variant for downsampled data). SBX (original/full-depth)
# uses the finalized sbx-r-f release, matching the manuscript-table
# methodology (see README.md -- "SBX pipeline finalized to sbx-r-f") --
# not sbx-r or the earlier plain sbx files, which live under
# data/_legacy_sbx_versions/, outside this convention. The downsampled
# set below is intentionally still on sbx-r -- no downsampled sbx-r-f
# SAGE file exists yet.
PIPELINES_ORIGINAL = {
    "Dragen": {
        "": (VCF_DIR / "Dragen/HG008-T-mosaic/HG008-T-mosaic.hard-filtered.vcf.gz", "HG008-T-mosaic"),
    },
    "Ilmn": {
        "SAGE": (VCF_DIR / "Ilmn/HG008T/HG008T.sage.somatic.vcf.gz", "HG008T"),
        "PURPLE": (VCF_DIR / "Ilmn/HG008T/HG008T.purple.somatic.vcf.gz", "HG008T"),
    },
    "SBX": {
        "SAGE": (VCF_DIR / "SBX/HG008T/HG008T.sage.somatic.vcf.gz", "HG008T"),
        "PURPLE": (VCF_DIR / "SBX/HG008T/HG008T.purple.somatic.vcf.gz", "HG008T"),
    },
    "Ultima": {
        "SAGE": (VCF_DIR / "Ultima/HG008T/HG008T.sage.somatic.vcf.gz", "HG008T"),
        "PURPLE": (VCF_DIR / "Ultima/HG008T/HG008T.purple.somatic.vcf.gz", "HG008T"),
    },
}
# No downsampled DRAGEN VCF exists -- omitted entirely, not just left blank.
PIPELINES_DOWNSAMPLED = {
    "Ilmn": {
        "SAGE": (VCF_DIR / "Ilmn_ds/HG008T/HG008T.sage.somatic.vcf.gz", "HG008T"),
        "PURPLE": (VCF_DIR / "Ilmn_ds/HG008T/HG008T.purple.somatic.vcf.gz", "HG008T"),
    },
    "SBX": {
        "SAGE": (VCF_DIR / "SBX_ds/HG008T/HG008T.sage.somatic.vcf.gz", "HG008T"),
        "PURPLE": (VCF_DIR / "SBX_ds/HG008T/HG008T.purple.somatic.vcf.gz", "HG008T"),
    },
    "Ultima": {
        "SAGE": (VCF_DIR / "Ultima_ds/HG008T/HG008T.sage.somatic.vcf.gz", "HG008T"),
        "PURPLE": (VCF_DIR / "Ultima_ds/HG008T/HG008T.purple.somatic.vcf.gz", "HG008T"),
    },
}

# "FP" = FILTER exactly one of these, at the pipeline's FINAL stage
# (PURPLE for the SAGE+PURPLE pipelines, the only stage for DRAGEN) --
# matches the finalized full-rescue methodology's definition of "counts
# as PASS." Whole-string match, not substring (e.g. NOT "PON;minMQF").
RESCUE_FILTERS = {"PASS", "PON", "PON;PONGnomad", "PONGnomad"}

# widest-tier run per pipeline -- source of the aardvark BD/EA columns
# (see module docstring for why widest-tier tags generalize to every tier)
WIDEST_RUN_ID_ORIGINAL = {
    "Dragen": "manuscript_dragen_widest",
    "Ilmn": "manuscript_ilmnpurple_widest",
    "SBX": "manuscript_sbxrf_widest",
    "Ultima": "manuscript_ultimapurple_widest",
}
WIDEST_RUN_ID_DOWNSAMPLED = {
    "Ilmn": "manuscript_ds_ilmnpurple_widest",
    "SBX": "manuscript_ds_sbxr_widest",
    "Ultima": "manuscript_ds_ultimapurple_widest",
}

REF_CONTIGS = None


def ref_contigs():
    global REF_CONTIGS
    if REF_CONTIGS is None:
        with open(f"{REF}.fai") as fh:
            REF_CONTIGS = ",".join(line.split("\t")[0] for line in fh)
    return REF_CONTIGS


def run(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def normalized_records(vcf_path, sample):
    """bcftools view (subset sample) | restrict to ref contigs | split multiallelics + atomize -> list of raw VCF lines."""
    cmd = (
        f"bcftools view -s {sample} '{vcf_path}' "
        f"| bcftools view -t '{ref_contigs()}' "
        f"| bcftools norm -m -any --atomize -f '{REF}' 2>/dev/null "
        f"| bcftools view -H"
    )
    out = run(cmd)
    return out.splitlines() if out else []


def build_lookup(vcf_path, sample, want_repc):
    """key (chrom,pos,ref,alt) -> dict(filter=..., af=..., repc=...)

    A caller can emit two candidate records that atomize to the same
    (chrom,pos,ref,alt) key (e.g. two slightly different local
    realignments of the same event, one PASS and one filtered). On
    collision, keep the PASS record -- matches PURPLE's own observed
    behavior of carrying forward the PASS SAGE candidate when
    duplicates exist. Ties (same PASS-ness) keep whichever was seen
    first.
    """
    lookup = {}
    for line in normalized_records(vcf_path, sample):
        f = line.split("\t")
        chrom, pos, _id, ref, alt, _qual, filt, info = f[:8]
        fmt_keys = f[8].split(":")
        fmt_vals = f[9].split(":")
        rec = dict(zip(fmt_keys, fmt_vals))
        af = rec.get("AF", "")
        repc = ""
        if want_repc:
            for kv in info.split(";"):
                if kv.startswith("REP_C="):
                    repc = kv.split("=", 1)[1]
                    break
        key = (chrom, pos, ref, alt)
        existing = lookup.get(key)
        if existing is not None and existing["filter"] == "PASS":
            continue  # already have the preferred PASS record, don't overwrite
        lookup[key] = {"filter": filt, "af": af, "repc": repc}
    return lookup


def build_bed_keyset(bed_path, vcf_path):
    """Set of (chrom,pos,ref,alt) keys from vcf_path whose FULL reference
    span (POS to POS+len(REF)-1) is entirely contained within this BED --
    not just any overlap.

    aardvark's own region-restriction requires full span containment: a
    multi-base deletion straddling the boundary between two adjacent BED
    intervals (e.g. one ends at position X, the next starts at X+2,
    leaving a 1bp gap the deletion's span crosses) is excluded from that
    run's truth.vcf.gz entirely. A naive any-overlap check (e.g.
    `bcftools view -R`, which tests the record's full REF span against
    the BED) instead counts it as "in region" whenever any part of the
    span touches either interval, which silently miscounted 391/10648
    truth records at the Ultima-HCR tier alone before this fix --
    verified by diffing this function's output against the actual `BD`
    tags in a real aardvark run's truth.vcf.gz at that tier.
    """
    out = run(f"bcftools view '{vcf_path}' -H")
    keys = []
    span_lines = []
    for line in out.splitlines():
        f = line.split("\t")
        chrom, pos, ref, alt = f[0], f[1], f[3], f[4]
        keys.append((chrom, pos, ref, alt))
        start = int(pos) - 1
        end = start + len(ref)
        span_lines.append(f"{chrom}\t{start}\t{end}\t{len(keys) - 1}\n")
    if not keys:
        return set()
    spans_path = "/tmp/_archive_spans_for_coverage.bed"
    with open(spans_path, "w") as fh:
        fh.writelines(span_lines)
    run(f"sort -k1,1 -k2,2n -o {spans_path} {spans_path}")
    cov_out = run(f"bedtools coverage -a {spans_path} -b '{bed_path}'")
    result = set()
    for line in cov_out.splitlines():
        f = line.rstrip("\n").split("\t")
        idx = int(f[3])
        frac = float(f[-1])
        if frac >= 0.999:  # full containment, allow for float rounding
            result.add(keys[idx])
    return result


def write_minimal_vcf(keys, out_path, chrom_rank):
    """Write a minimal sorted, bgzipped+tabix-indexed VCF for these keys (for -R region lookups)."""
    rows = sorted(keys, key=lambda k: (chrom_rank.get(k[0], 999), int(k[1])))
    plain = "/tmp/_archive_fp_keys.vcf"
    with open(plain, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in rows:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    run(f"bgzip -f {plain}")
    run(f"mv {plain}.gz {out_path}")
    run(f"tabix -f -p vcf {out_path}")


def load_aardvark_tags(run_id):
    """(chrom,pos,ref,alt) -> (BD,EA) from a run's aardvark truth.vcf.gz
    and query.vcf.gz. See module docstring: this is the source of truth
    for TP/FN/FP status in the archive, not the FILTER/AF columns."""
    run_dir = RESULTS_DIR / "runs" / run_id / "aardvark"

    def load(vcf_name):
        d = {}
        out = run(f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%BD]\\t[%EA]\\n' '{run_dir / vcf_name}'")
        for line in out.splitlines():
            c, p, r, a, bd, ea = line.split("\t")
            d[(c, p, r, a)] = (bd, ea)
        return d

    return load("truth.vcf.gz"), load("query.vcf.gz")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--downsampled", action="store_true",
                         help="use the ~30x/90x downsampled VCFs (3 pipelines, no DRAGEN) instead of the original full-depth VCFs")
    args = parser.parse_args()

    if args.downsampled:
        pipelines = PIPELINES_DOWNSAMPLED
        widest_run_id = WIDEST_RUN_ID_DOWNSAMPLED
        out_tsv = RESULTS_DIR / "Variant_archive_truthset_downsampled.tsv"
    else:
        pipelines = PIPELINES_ORIGINAL
        widest_run_id = WIDEST_RUN_ID_ORIGINAL
        out_tsv = RESULTS_DIR / "Variant_archive_truthset.tsv"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    sage_purple_pipelines = [p for p in pipelines if p != "Dragen"]
    final_stage = {p: ("" if p == "Dragen" else "PURPLE") for p in pipelines}

    with open(f"{REF}.fai") as fh:
        chrom_rank = {line.split("\t")[0]: i for i, line in enumerate(fh)}

    print("[1/6] Splitting truth VCF multiallelics...", file=sys.stderr)
    truth_split = "/tmp/truth_split_for_archive.vcf.gz"
    run(
        f"bcftools norm -m -any -f '{REF}' '{TRUTH_VCF}' 2>/dev/null "
        f"| bcftools sort -Oz -o {truth_split} 2>/dev/null"
    )
    run(f"tabix -f -p vcf {truth_split}")
    truth_records = run(f"bcftools view -H {truth_split}").splitlines()
    truth_keys = set()
    for line in truth_records:
        f = line.split("\t")
        truth_keys.add((f[0], f[1], f[3], f[4]))
    print(f"  {len(truth_records)} truth alleles after splitting", file=sys.stderr)

    print("[2/6] Building per-pipeline lookups from raw VCFs...", file=sys.stderr)
    lookups = {}  # (pipeline, stage) -> dict
    for pipeline, stages in pipelines.items():
        for stage, (vcf_path, sample) in stages.items():
            want_repc = pipeline != "Dragen"
            print(f"  {pipeline} {stage or '(single-stage)'} ({vcf_path.name})...", file=sys.stderr)
            lookups[(pipeline, stage)] = build_lookup(vcf_path, sample, want_repc)

    print("[3/6] Finding cross-tech FP candidates absent from truth...", file=sys.stderr)
    fp_keys = set()
    for pipeline, stage in final_stage.items():
        d = lookups[(pipeline, stage)]
        n_before = len(fp_keys)
        for key, rec in d.items():
            if rec["filter"] in RESCUE_FILTERS:
                fp_keys.add(key)
        print(f"  {pipeline}: cumulative union now {len(fp_keys)} keys (+{len(fp_keys)-n_before})", file=sys.stderr)
    fp_only_keys = fp_keys - truth_keys
    print(f"  {len(fp_keys)} total PASS/PON/PONGnomad candidates across all included techs, "
          f"{len(fp_only_keys)} not in truth set", file=sys.stderr)

    print("[4/6] Building a minimal VCF for the FP-only keys (for region lookups)...", file=sys.stderr)
    fp_vcf = "/tmp/fp_only_for_archive.vcf.gz"
    write_minimal_vcf(fp_only_keys, fp_vcf, chrom_rank)

    print("[5/6] Computing region membership (truth rows + FP rows)...", file=sys.stderr)
    bed_keysets = {}
    for label, bed_path in BEDS.items():
        if not bed_path.exists():
            print(f"  WARNING: missing {bed_path}, skipping {label}", file=sys.stderr)
            bed_keysets[label] = set()
            continue
        truth_hits = build_bed_keyset(bed_path, truth_split)
        fp_hits = build_bed_keyset(bed_path, fp_vcf)
        bed_keysets[label] = truth_hits | fp_hits
        print(f"  {label}: {len(truth_hits)} truth + {len(fp_hits)} FP alleles overlap", file=sys.stderr)

    print("[5.5/6] Loading aardvark's own BD/EA verdicts from widest-tier runs...", file=sys.stderr)
    aardvark_tags = {}  # pipeline -> (truth_bdea, query_bdea)
    for pipeline, run_id in widest_run_id.items():
        print(f"  {pipeline} ({run_id})...", file=sys.stderr)
        aardvark_tags[pipeline] = load_aardvark_tags(run_id)

    print("[6/6] Writing archive TSV...", file=sys.stderr)
    header = ["CHROM", "POS", "REF", "ALT", "Variant_Type", "in_truth_set"]
    header += list(BEDS.keys())
    header += ["in_homopolymer"]
    if "Dragen" in pipelines:
        header += ["Dragen_Filter", "AF_Dragen"]
    for pipeline in sage_purple_pipelines:
        header += [
            f"{pipeline}_SAGE_Filter",
            f"{pipeline}_PURPLE_Filter",
            f"AF_{pipeline}_SAGE",
            f"REP_C_{pipeline}_SAGE",
        ]
    for pipeline in pipelines:  # respects dict order; excludes Dragen when downsampled
        header += [f"{pipeline}_BD", f"{pipeline}_EA"]

    def variant_type(ref, alt):
        if len(ref) == 1 and len(alt) == 1:
            return "SNV"
        if len(alt) > len(ref):
            return "INS"
        return "DEL"

    def row_for(chrom, pos, ref, alt, is_truth):
        key = (chrom, pos, ref, alt)
        sage_by_pipeline = {p: lookups[(p, "SAGE")].get(key) for p in sage_purple_pipelines}
        purple_by_pipeline = {p: lookups[(p, "PURPLE")].get(key) for p in sage_purple_pipelines}
        # "in a homopolymer" if ANY pipeline's SAGE call carries a REP_C annotation there
        in_homopolymer = any(
            sage_by_pipeline[p] and sage_by_pipeline[p]["repc"] != "" for p in sage_purple_pipelines
        )

        row = [chrom, pos, ref, alt, variant_type(ref, alt), "1" if is_truth else "0"]
        for label in BEDS:
            row.append("1" if key in bed_keysets[label] else "0")
        row.append("1" if in_homopolymer else "0")
        if "Dragen" in pipelines:
            d = lookups[("Dragen", "")].get(key)
            row.append(d["filter"] if d else "")
            row.append(d["af"] if d else "")
        for pipeline in sage_purple_pipelines:
            sage = sage_by_pipeline[pipeline]
            purple = purple_by_pipeline[pipeline]
            row.append(sage["filter"] if sage else "")
            row.append(purple["filter"] if purple else "")
            row.append(sage["af"] if sage else "")
            row.append(sage["repc"] if sage else "")
        for pipeline in pipelines:
            truth_bdea, query_bdea = aardvark_tags[pipeline]
            bd, ea = (truth_bdea if is_truth else query_bdea).get(key, ("", ""))
            row.append(bd)
            row.append(ea)
        return row

    with open(out_tsv, "w") as out:
        out.write("\t".join(header) + "\n")
        for line in truth_records:
            f = line.split("\t")
            out.write("\t".join(row_for(f[0], f[1], f[3], f[4], True)) + "\n")
        for chrom, pos, ref, alt in sorted(
            fp_only_keys, key=lambda k: (chrom_rank.get(k[0], 999), int(k[1]))
        ):
            out.write("\t".join(row_for(chrom, pos, ref, alt, False)) + "\n")

    total = len(truth_records) + len(fp_only_keys)
    print(f"Wrote {total} rows ({len(truth_records)} truth + {len(fp_only_keys)} FP-only) to {out_tsv}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
