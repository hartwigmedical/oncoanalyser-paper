#!/usr/bin/env python3
"""Builds the manuscript TSVs (original + downsampled, ALL + per-variant-
type) from the results.tsv ledger, for a given list of run_ids. Called
by generate_manuscript_table.sh -- not meant to be run standalone unless
you already know which run_ids you want.

Columns: pipeline, bed_tier, truth_TP, query_TP, FP, FN, recall, precision, f1.

Both TP columns are reported because recall and precision are counted on
opposite sides of the comparison and are NOT two views of one TP number:

  recall    = truth_TP / (truth_TP + FN)   -- truth alleles aardvark matched
  precision = query_TP / (query_TP + FP)   -- query records aardvark matched

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


REQUIRED_COLS = ["run_id", "tech", "bed_tier", "variant_type", "dosage_scoring",
                 "truth_tp", "truth_fn", "query_tp", "query_fp",
                 "metric_recall", "metric_precision", "metric_f1"]


def load_rows(results_tsv, wanted_run_ids):
    """variant_type -> {run_id: fields}, plus the column-name -> index map."""
    rows = {}
    with open(results_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        absent = [c for c in REQUIRED_COLS if c not in idx]
        if absent:
            sys.exit(f"{results_tsv} is missing required column(s) {absent}.\n"
                     f"  header: {header}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[idx["run_id"]] not in wanted_run_ids:
                continue
            if f[idx["dosage_scoring"]] != "ea_corrected":
                continue
            vt = f[idx["variant_type"]]
            rows.setdefault(vt, {})[f[idx["run_id"]]] = f
    return rows, idx


def pct(value, where):
    """A ledger metric as a percentage string, or 'NA'.

    append_ledger_rows.py writes an empty metric when aardvark reported no GT
    row for that variant type (a narrow region tier can contain no insertions),
    so these fields are legitimately blank sometimes. float() on a blank used to
    abort the whole table build with a bare ValueError.
    """
    try:
        return f"{float(value) * 100:.2f}"
    except (TypeError, ValueError):
        print(f"WARNING: {where}: metric is not a number ({value!r}), writing NA", file=sys.stderr)
        return "NA"


def count(value):
    """A ledger count verbatim, or 'NA' when blank -- matching pct()'s output for
    the metrics, so a row with no truth-side numbers reads the same in every
    column instead of leaving empty TSV fields."""
    value = (value or "").strip()
    return value if value else "NA"


def write_table(rows, idx, out_tsv, run_ids, label):
    for vt, suffix in VARIANT_TYPES:
        vt_rows = rows.get(vt, {})
        missing = [rid for rid in run_ids if rid not in vt_rows]
        if missing:
            print(f"WARNING [{label}]: missing {vt} ledger rows for: {missing}", file=sys.stderr)
        this_out = out_tsv if suffix == "" else out_tsv.replace(".tsv", f"{suffix}.tsv")
        with open(this_out, "w") as out:
            out.write("pipeline\tbed_tier\ttruth_TP\tquery_TP\tFP\tFN\trecall\tprecision\tf1\n")
            n = 0
            for rid in run_ids:
                if rid not in vt_rows:
                    continue
                f = vt_rows[rid]
                tech = f[idx["tech"]]
                bed_tier = f[idx["bed_tier"]]
                # recall pairs with truth_tp/truth_fn, precision with query_tp/query_fp
                # -- see the module docstring on why both TP columns are printed.
                truth_tp = count(f[idx["truth_tp"]])
                query_tp = count(f[idx["query_tp"]])
                fn = count(f[idx["truth_fn"]])
                fp = count(f[idx["query_fp"]])
                where = f"[{label}] {vt} {rid}"
                recall = pct(f[idx["metric_recall"]], f"{where} recall")
                precision = pct(f[idx["metric_precision"]], f"{where} precision")
                f1 = pct(f[idx["metric_f1"]], f"{where} f1")
                out.write(
                    f"{PIPELINE_LABELS.get(tech, tech)}\t{TIER_LABELS.get(bed_tier, bed_tier)}\t"
                    f"{truth_tp}\t{query_tp}\t{fp}\t{fn}\t{recall}\t{precision}\t{f1}\n"
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
