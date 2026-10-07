#!/usr/bin/env python3

"""Compare query loops with reference loops after PAF-based coordinate liftover.

Requirements: Python 3, intervaltree, and rustybam available on PATH (Unix/Linux).
Use the software versions from the original analysis when reproducing results.

Example (replace all /path/to/... entries):
    python compare_loops_liftover.py /path/to/reference.loops.bedpe \
        /path/to/query.loops.bedpe /path/to/alignment.paf \
        /path/to/loop_comparison.tsv --tolerance 10000

Loop files: tab-delimited text with at least six columns:
    chrom1  start1  end1  chrom2  start2  end2  [optional extra columns]
Coordinates must be 0-based, half-open, with start < end. Blank lines and lines
starting with # are ignored. Other header lines must be removed before use.
HiCExplorer loop files are accepted when their first six columns follow this
layout, regardless of the filename extension. Extra query columns are retained.

Liftover direction: rustybam liftover --bed projects PAF target to query.
The script's query_file must use the PAF target coordinates (PAF column 6), and
ref_file must use the PAF query coordinates (PAF column 1). These names refer to
two different roles: script query/reference files versus PAF query/target sides.
Chromosome names must match their respective PAF sides.
Use a PAF with the CIGAR tags required by rustybam (generated with -c --eqx).
The rustybam output must retain id:Z:<BED-column-4> tags to identify each anchor.
The script reads the lifted PAF query coordinates from columns 1, 3 and 4.

Matching procedure:
* Expand each reference anchor by --tolerance bp on each side (default 10000).
  Find overlaps between lifted query intervals and expanded reference anchors.
* Collect the reference loop IDs matched by each query anchor.
* BothMatch: the two query anchors have at least one matched reference loop ID
  in common.
* SingleMatch: no common reference loop ID, but at least one query anchor has
  a match.
* NoMatch: neither query anchor has a match.
All returned liftover segments are considered. Anchors without a liftover result
are treated as unmatched. Each query loop receives one classification relative
to the reference loop set.

Output: each original query data row with a tab and label appended; no header.
The output parent directory must exist; an existing output file is overwritten.
"""

import subprocess
import argparse
from intervaltree import IntervalTree
from collections import defaultdict

# =========================
# Command-line arguments
# =========================
parser = argparse.ArgumentParser(description="Compare query loops to reference loops via liftover")
parser.add_argument("ref_file", help="Reference loops: plain-text BEDPE-like file, at least 6 columns (PAF query coordinates)")
parser.add_argument("query_file", help="Query loops: plain-text BEDPE-like file, at least 6 columns (PAF target coordinates)")
parser.add_argument("paf_file", help="PAF with required CIGAR tags; target=query_file genome, query=ref_file genome")
parser.add_argument("output_file", help="Output TSV: original query columns plus BothMatch/SingleMatch/NoMatch, no header")
parser.add_argument("--tolerance", type=int, default=10000, help="Expand reference anchors by this many bp on each side (default: 10000)")
args = parser.parse_args()
if args.tolerance < 0:
    parser.error("--tolerance must be nonnegative")

ref_file = args.ref_file
query_file = args.query_file
paf_file = args.paf_file
output_file = args.output_file
TOLERANCE = args.tolerance


# =========================
# 1. Load reference loops into interval trees
# =========================
def load_reference_loops(ref_file):
    trees = {}
    loop_id = 0

    with open(ref_file) as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cols = line.strip().split()
            chr1, s1, e1, chr2, s2, e2 = cols[:6]

            s1, e1 = int(s1), int(e1)
            s2, e2 = int(s2), int(e2)

            if chr1 not in trees:
                trees[chr1] = IntervalTree()
            if chr2 not in trees:
                trees[chr2] = IntervalTree()

            trees[chr1].addi(s1 - TOLERANCE, e1 + TOLERANCE, loop_id)
            trees[chr2].addi(s2 - TOLERANCE, e2 + TOLERANCE, loop_id)

            loop_id += 1

    return trees


# =========================
# 2. Lift all query anchors in one rustybam call
# =========================
def liftover_batch(regions, paf_file):

    if not regions:
        return defaultdict(list)

    bed_lines = []
    for i, (c, s, e) in enumerate(regions):
        bed_lines.append(f"{c}\t{s}\t{e}\t{i}")

    bed_input = "\n".join(bed_lines) + "\n"

    cmd = [
        "rustybam",
        "liftover",
        "--bed", "/dev/stdin",
        paf_file
    ]

    result = subprocess.run(
        cmd,
        input=bed_input,
        text=True,
        capture_output=True
    )

    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError("rustybam failed")

    mapping_dict = defaultdict(list)
    total_lines = 0
    skipped_lines = 0

    for line in result.stdout.strip().split("\n"):
        if not line.strip():
            continue

        total_lines += 1
        cols = line.split()

        # PAF query coordinates are the lifted coordinates in ref_file space.
        lifted_chr = cols[0]
        lifted_start = int(cols[2])
        lifted_end = int(cols[3])

        idx = None
        for c in cols:
            if c.startswith("id:Z:"):
                parts = c.split(":")
                try:
                    idx = int(parts[-1])
                except ValueError:
                    continue
                break

        if idx is None:
            skipped_lines += 1
            continue

        mapping_dict[idx].append((lifted_chr, lifted_start, lifted_end))

    print(f"Liftover: {total_lines} lines parsed, {skipped_lines} lines skipped (no id:Z: field)")

    # Report how many original anchors have no liftover result.
    total_regions = len(regions)
    mapped_regions = len(mapping_dict)
    unmapped_regions = total_regions - mapped_regions
    print(f"Liftover: {mapped_regions}/{total_regions} anchors mapped, {unmapped_regions} anchors unmapped")

    return mapping_dict


# =========================
# 3. Find reference loop IDs hit by each query anchor
# =========================
def endpoint_match_loops(mappings, trees):
    matched_loops = set()

    for (chr_, start, end) in mappings:
        if chr_ not in trees:
            continue

        overlaps = trees[chr_].overlap(start, end)

        for iv in overlaps:
            matched_loops.add(iv.data)

    return matched_loops


# =========================
# 4. Compare query loops with the reference loop set
# =========================
def main():
    print("Loading reference loops...")
    trees = load_reference_loops(ref_file)

    print("Reading query loops...")
    query_loops = []
    regions = []

    with open(query_file) as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cols = line.strip().split()

            chr1, s1, e1 = cols[0], int(cols[1]), int(cols[2])
            chr2, s2, e2 = cols[3], int(cols[4]), int(cols[5])

            idx1 = len(regions)
            regions.append((chr1, s1, e1))

            idx2 = len(regions)
            regions.append((chr2, s2, e2))

            query_loops.append((line.strip(), idx1, idx2))

    print("Running liftover...")
    mapping_dict = liftover_batch(regions, paf_file)

    print("Comparing loops...")
    with open(output_file, "w") as fout:
        for line, idx1, idx2 in query_loops:

            map1 = mapping_dict.get(idx1, [])
            map2 = mapping_dict.get(idx2, [])

            loops1 = endpoint_match_loops(map1, trees)
            loops2 = endpoint_match_loops(map2, trees)

            # Find reference loop IDs shared by the two query-anchor match sets.
            common = loops1 & loops2

            if common:
                label = "BothMatch"
            elif loops1 or loops2:
                label = "SingleMatch"
            else:
                label = "NoMatch"

            fout.write(line + "\t" + label + "\n")

    print("Done!")


if __name__ == "__main__":
    main()
