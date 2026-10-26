#!/usr/bin/env python3
"""Builds a BED of reference tandem repeats long enough that SAGE would
annotate a variant there with REP_C > 8: repeat unit 1-5bp, >= 9 copies,
padded 1bp each side. Called by generate_manuscript_table.sh, which sorts
and merges the output into data/reference_data/regions/repc_gt8_mask.bed
and subtracts it from each confidence tier to make the *_repcle8 tiers
(see README.md, "Region tiers").

Usage: build_repc_mask.py <reference.fa.gz> <out_raw.bed>
(output is unsorted/unmerged, one row per run: chrom, start, end, unit, copies)

A region mask rather than a filter on SAGE's REP_C INFO field, because
REP_C only exists on records SAGE called: truth FNs and DRAGEN calls carry
none, so filtering on the field would drop variants unevenly across
pipelines. Validated against SAGE's own REP_C on every Illumina/SBX/Ultima
PURPLE call: 100% of REP_C > 8 calls fall inside the mask. Units that are
themselves repeats of a shorter unit (e.g. ATAT) are skipped, since the
shorter unit already covers them. Primary chromosomes only. Takes ~18 min
on the full GRCh38 reference.
"""
import gzip
import re
import sys

MINCOUNT = 9  # REP_C > 8
PAD = 1       # variant anchor base sits 1bp before the repeat it's in
MAX_UNIT = 5  # SAGE's REP_S never exceeds 5bp
CHROMS = {f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]}


def scan(name, seq, out):
    n = 0
    for k in range(1, MAX_UNIT + 1):
        pat = re.compile(r"([ACGT]{%d})\1{%d,}" % (k, MINCOUNT - 1))
        pos = 0
        while True:
            m = pat.search(seq, pos)
            if not m:
                break
            u = m.group(1)
            reducible = any(k % d == 0 and u == u[:d] * (k // d) for d in range(1, k))
            if reducible:
                pos = m.start() + 1
                continue
            out.write(f"{name}\t{max(0, m.start() - PAD)}\t{m.end() + PAD}\t{u}\t{(m.end() - m.start()) // k}\n")
            n += 1
            pos = m.end() - k + 1  # allow an overlapping run in another phase
    return n


def main():
    ref, out_bed = sys.argv[1], sys.argv[2]
    with gzip.open(ref, "rt") as fh, open(out_bed, "w") as out:
        name, buf = None, []

        def flush():
            if name in CHROMS:
                print(name, scan(name, "".join(buf).upper(), out), file=sys.stderr, flush=True)

        for line in fh:
            if line.startswith(">"):
                flush()
                name, buf = line[1:].split()[0], []
            elif name in CHROMS:
                buf.append(line.rstrip())
        flush()


if __name__ == "__main__":
    main()
