#!/usr/bin/env python3
"""Builds the manuscript TSVs (original + downsampled, ALL + per-variant-
type) from the results.tsv ledger, for a given list of run_ids. Called
by generate_manuscript_table.sh -- not meant to be run standalone unless
you already know which run_ids you want.
"""
import argparse
import sys

TIER_LABELS = {
    "all": "Widest truth set",
    "all_giabhcr": "GIAB HCR",
    "all_hcr": "Ultima HCR",
    "nogermlinewithin50bp_giabhcr": "GIAB HCR + nogermlinewithin50bp",
    "nogermlinewithin50bp_hcr": "Ultima HCR + nogermlinewithin50bp",
}
PIPELINE_LABELS = {
    "illumina_dragen": "ILMN DRAGEN",
    "illumina_purple": "ILMN SAGE/PURPLE",
    "sbx_purple": "SBX SAGE/PURPLE",
    "ultima_purple": "Ultima SAGE/PURPLE",
    "illumina_purple_ds": "ILMN SAGE/PURPLE",
    "sbx_purple_ds": "SBX SAGE/PURPLE",
    "ultima_purple_ds": "Ultima SAGE/PURPLE",
}
# variant_type in results.tsv -> (suffix for output filename, or "" for the main table)
VARIANT_TYPES = [("ALL", ""), ("Snv", "_SNV"), ("Insertion", "_Insertion"), ("Deletion", "_Deletion")]


def load_rows(results_tsv, wanted_run_ids):
    """variant_type -> {run_id: fields}, plus the column-name -> index map."""
    rows = {}
    with open(results_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[idx["run_id"]] not in wanted_run_ids:
                continue
            if f[idx["dosage_scoring"]] != "ea_corrected":
                continue
            vt = f[idx["variant_type"]]
            rows.setdefault(vt, {})[f[idx["run_id"]]] = f
    return rows, idx


def write_table(rows, idx, out_tsv, run_ids, label):
    for vt, suffix in VARIANT_TYPES:
        vt_rows = rows.get(vt, {})
        missing = [rid for rid in run_ids if rid not in vt_rows]
        if missing:
            print(f"WARNING [{label}]: missing {vt} ledger rows for: {missing}", file=sys.stderr)
        this_out = out_tsv if suffix == "" else out_tsv.replace(".tsv", f"{suffix}.tsv")
        with open(this_out, "w") as out:
            out.write("pipeline\tbed_tier\tTP\tFP\tFN\trecall\tprecision\tf1\n")
            n = 0
            for rid in run_ids:
                if rid not in vt_rows:
                    continue
                f = vt_rows[rid]
                tech = f[idx["tech"]]
                bed_tier = f[idx["bed_tier"]]
                tp = f[idx["truth_tp"]]
                fn = f[idx["truth_fn"]]
                fp = f[idx["query_fp"]]
                recall = float(f[idx["metric_recall"]]) * 100
                precision = float(f[idx["metric_precision"]]) * 100
                f1 = float(f[idx["metric_f1"]]) * 100
                out.write(
                    f"{PIPELINE_LABELS.get(tech, tech)}\t{TIER_LABELS.get(bed_tier, bed_tier)}\t"
                    f"{tp}\t{fp}\t{fn}\t{recall:.2f}\t{precision:.2f}\t{f1:.2f}\n"
                )
                n += 1
        print(f"Wrote {n} rows to {this_out}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-tsv", required=True, help="path to results.tsv (the ledger)")
    parser.add_argument("--out-tsv", required=True,
                         help="path for the original-dataset manuscript table; "
                              "the downsampled table is written alongside it with a _downsampled suffix")
    parser.add_argument("--original-run-ids", nargs="+", required=True)
    parser.add_argument("--downsampled-run-ids", nargs="+", required=True)
    args = parser.parse_args()

    wanted = set(args.original_run_ids) | set(args.downsampled_run_ids)
    rows, idx = load_rows(args.results_tsv, wanted)

    write_table(rows, idx, args.out_tsv, args.original_run_ids, "original")
    write_table(rows, idx, args.out_tsv.replace(".tsv", "_downsampled.tsv"), args.downsampled_run_ids, "downsampled")


if __name__ == "__main__":
    main()
