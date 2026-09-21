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

# The variant_type labels this script expects aardvark's summary.tsv to use. The
# truth-side numbers (truth_total/truth_tp/truth_fn/metric_recall) of every
# ea_corrected row are looked up by these exact strings, so if aardvark ever
# renames them the lookups would all miss and the rows would come out blank
# rather than wrong-looking -- hence the explicit check in validate_variant_types.
VARIANT_TYPES = ["ALL", "Snv", "Insertion", "Deletion"]


def classify(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "Snv"
    if len(alt) > len(ref) and alt.startswith(ref):
        return "Insertion"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "Deletion"
    return None


def validate_variant_types(summary_rows, summary_tsv):
    """Check aardvark's GT variant_type labels against VARIANT_TYPES.

    A missing per-type row is legitimate (a narrow region tier can genuinely
    contain no insertions), so that is only reported. A missing "ALL" row, or a
    label set that shares nothing with what we expect, means the labels have
    moved and every truth-side lookup below would silently blank out -- fail.
    """
    gt_types = {r["variant_type"] for r in summary_rows if r["comparison"] == "GT"}
    if not gt_types:
        sys.exit(f"{summary_tsv} has no comparison=GT rows -- cannot build ea_corrected rows.")
    if "ALL" not in gt_types:
        sys.exit(
            f"{summary_tsv} has no comparison=GT row with variant_type='ALL' -- the main "
            f"manuscript table needs it.\n  found: {sorted(gt_types)}"
        )
    missing = [vt for vt in VARIANT_TYPES if vt not in gt_types]
    unknown = sorted(gt_types - set(VARIANT_TYPES))
    # A label we expect is absent *while* aardvark reports labels we don't know:
    # that is a rename, not an empty variant class, and every truth-side lookup
    # for those types would blank out. Fail, because these notes go to stderr,
    # which generate_manuscript_table.sh sends to a per-run log the user is not
    # watching -- a warning here would be effectively silent.
    if missing and unknown:
        sys.exit(
            f"{summary_tsv} variant_type labels have changed -- refusing to write rows whose "
            "truth-side columns would all be blank.\n"
            f"  expected but absent: {missing}\n  present but unknown: {unknown}\n"
            "  Update VARIANT_TYPES in append_ledger_rows.py to match aardvark."
        )
    if missing:
        # No unknown labels, so this region tier genuinely contains none of these.
        print(f"NOTE: {summary_tsv} has no GT row for {missing}; those ledger rows get blank "
              "truth-side columns and metrics.", file=sys.stderr)
    if unknown:
        print(f"NOTE: aardvark also reported variant_type {unknown}, which this ledger does not "
              "table.", file=sys.stderr)


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
    if not summary_rows:
        sys.exit(f"{summary_tsv} has no rows -- aardvark produced no summary for this run.")
    summary_fields = [c for c in summary_rows[0].keys() if c != "compare_label"]
    validate_variant_types(summary_rows, summary_tsv)

    fieldnames = FIXED_COLS + summary_fields + TRAILING_COLS

    if results_tsv.exists():
        run_id_col = fieldnames.index("run_id")
        with open(results_tsv, newline="") as fh:
            existing_header = fh.readline().rstrip("\n").split("\t")
            if existing_header != fieldnames:
                sys.exit(
                    "results.tsv header does not match this run's expected columns -- "
                    "ledger schema drift detected, refusing to append.\n"
                    f"  existing: {existing_header}\n  expected: {fieldnames}"
                )
            # The ledger is append-only, so a re-run of a run_id that already has rows
            # would double it up. build_manuscript_tables.py keys on run_id and would
            # silently keep whichever copy came last, so refuse instead.
            already = sum(1 for line in fh
                          if line.split("\t")[run_id_col:run_id_col + 1] == [run_id])
        if already:
            sys.exit(
                f"{results_tsv} already holds {already} row(s) for run_id '{run_id}' -- "
                "refusing to append duplicates.\n"
                f"  Remove those rows (and results/runs/{run_id}/) to redo this run."
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

        counts = {vt: {"tp": 0, "novel_fp": 0, "dosage_fp": 0} for vt in VARIANT_TYPES}
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

        for vt in VARIANT_TYPES:
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
