#!/usr/bin/env python3
"""Appends one run's results to the results.tsv ledger: the raw
GT/BASEPAIR rows straight from aardvark's summary.tsv, plus EA-aware
dosage-corrected GT rows.

Many "FP" calls under GT-mode scoring are not independent false calls:
they're sites aardvark's truth side already counts as TP, where the
query's reported genotype claims more alt alleles than truth expects
(EA=1,OA=2 -- typically a real heterozygous somatic variant plus tumor
LOH, which collapses to a homozygous-looking call from bulk reads). A
record with EA=0 has no truth expectation at that position at all --
that is a genuine, independent false positive. dosage_fp (EA>=1)
records are truth-matched correct calls, not errors, so they're added
to TP (not just dropped from FP -- dropping from both sides would
desync TP from truth_tp). See RESULTS.md for the full mechanism.

Before appending, the on-disk ledger header (if the file already
exists) is checked against this run's column set -- if a previous run
appended under a different schema (a column added/removed/renamed
since), this refuses to append rather than silently writing
misaligned rows.

Usage: append_ledger_rows.py <run_dir> <results_tsv> <run_id> <tech_label> <query_vcf> <query_sample> <bed_tier>
"""
import csv
import subprocess
import sys
from pathlib import Path

FIXED_COLS = ["run_id", "tech", "date", "bed_tier", "query_vcf", "query_sample"]
TRAILING_COLS = ["tier", "fn_root_cause", "dosage_scoring"]


def classify(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "Snv"
    if len(alt) > len(ref) and alt.startswith(ref):
        return "Insertion"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "Deletion"
    return None


def run_date(query_vcf):
    proc = subprocess.run(["date", "-r", query_vcf, "+%Y-%m-%d"], capture_output=True, text=True)
    return proc.stdout.strip() if proc.returncode == 0 else "unknown"


def main():
    run_dir, results_tsv, run_id, tech_label, query_vcf, query_sample, bed_tier = sys.argv[1:8]
    run_dir = Path(run_dir)
    results_tsv = Path(results_tsv)
    summary_tsv = run_dir / "aardvark" / "summary.tsv"
    query_vcf_scored = run_dir / "aardvark" / "query.vcf.gz"

    with open(summary_tsv, newline="") as fh:
        summary_rows = list(csv.DictReader(fh, delimiter="\t"))
    summary_fields = [c for c in summary_rows[0].keys() if c != "compare_label"] if summary_rows else []

    fieldnames = FIXED_COLS + summary_fields + TRAILING_COLS

    if results_tsv.exists():
        with open(results_tsv, newline="") as fh:
            existing_header = fh.readline().rstrip("\n").split("\t")
        if existing_header != fieldnames:
            sys.exit(
                "results.tsv header does not match this run's expected columns -- "
                "ledger schema drift detected, refusing to append.\n"
                f"  existing: {existing_header}\n  expected: {fieldnames}"
            )
        write_header = False
    else:
        write_header = True

    date_str = run_date(query_vcf)

    with open(results_tsv, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        if write_header:
            writer.writeheader()

        # -- raw rows: one per summary.tsv row, straight passthrough --
        for row in summary_rows:
            out = {
                "run_id": run_id, "tech": tech_label, "date": date_str,
                "bed_tier": bed_tier, "query_vcf": query_vcf, "query_sample": query_sample,
            }
            out.update({k: row[k] for k in summary_fields})
            out["tier"] = "ALL"
            out["fn_root_cause"] = "ALL"
            out["dosage_scoring"] = "raw"
            writer.writerow(out)

        # -- EA-aware dosage-corrected GT rows --
        truth_stats = {}
        for row in summary_rows:
            if row["comparison"] != "GT":
                continue
            truth_stats[row["variant_type"]] = {
                "truth_total": row["truth_total"],
                "truth_tp": row["truth_tp"],
                "truth_fn": row["truth_fn"],
                "metric_recall": row["metric_recall"],
            }

        counts = {vt: {"tp": 0, "novel_fp": 0, "dosage_fp": 0} for vt in ["ALL", "Snv", "Insertion", "Deletion"]}
        proc = subprocess.run(["bcftools", "view", "-H", str(query_vcf_scored)], capture_output=True, text=True, check=True)
        for line in proc.stdout.splitlines():
            f = line.rstrip("\n").split("\t")
            ref, alt, fmt, sample = f[3], f[4], f[8], f[9]
            rec = dict(zip(fmt.split(":"), sample.split(":")))
            bd, ea = rec.get("BD"), rec.get("EA")
            if bd not in ("TP", "FP"):
                continue
            vt = classify(ref, alt)
            for t in ["ALL"] + ([vt] if vt else []):
                if bd == "TP":
                    counts[t]["tp"] += 1
                elif ea == "0":
                    counts[t]["novel_fp"] += 1
                else:
                    counts[t]["dosage_fp"] += 1

        for vt in ["ALL", "Snv", "Insertion", "Deletion"]:
            c = counts[vt]
            tp = c["tp"] + c["dosage_fp"]
            novel_fp = c["novel_fp"]
            query_total = tp + novel_fp
            precision = tp / query_total if query_total else ""
            ts = truth_stats.get(vt, {"truth_total": "", "truth_tp": "", "truth_fn": "", "metric_recall": ""})
            recall = ts["metric_recall"]
            f1 = ""
            if precision != "" and recall not in ("", None):
                try:
                    r = float(recall)
                    if (r + precision) > 0:
                        f1 = 2 * r * precision / (r + precision)
                except ValueError:
                    pass
            out = {
                "run_id": run_id, "tech": tech_label, "date": date_str,
                "bed_tier": bed_tier, "query_vcf": query_vcf, "query_sample": query_sample,
                "comparison": "GT", "region_label": "ALL", "filter": "ALL", "variant_type": vt,
                "truth_total": ts["truth_total"], "truth_tp": ts["truth_tp"], "truth_fn": ts["truth_fn"],
                "query_total": str(query_total), "query_tp": str(tp), "query_fp": str(novel_fp),
                "metric_recall": str(recall), "metric_precision": str(precision), "metric_f1": str(f1),
                "truth_fn_gt": "", "query_fp_gt": "",
                "tier": "ALL", "fn_root_cause": "ALL", "dosage_scoring": "ea_corrected",
            }
            writer.writerow(out)


if __name__ == "__main__":
    main()
