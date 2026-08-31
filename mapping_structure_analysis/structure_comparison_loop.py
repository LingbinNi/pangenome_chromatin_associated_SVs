#!/usr/bin/env python3

import subprocess
import sys
import argparse
from intervaltree import IntervalTree
from collections import defaultdict

# =========================
# 参数解析
# =========================
parser = argparse.ArgumentParser(description="Compare query loops to reference loops via liftover")
parser.add_argument("ref_file", help="Reference loop file (BED-like, 6 columns)")
parser.add_argument("query_file", help="Query loop file (BED-like, 6 columns)")
parser.add_argument("paf_file", help="PAF alignment file for liftover")
parser.add_argument("output_file", help="Output file")
parser.add_argument("--tolerance", type=int, default=10000, help="Tolerance in bp for anchor matching (default: 10000)")
args = parser.parse_args()

ref_file = args.ref_file
query_file = args.query_file
paf_file = args.paf_file
output_file = args.output_file
TOLERANCE = args.tolerance


# =========================
# 1. Reference loops → interval tree
# =========================
def load_reference_loops(ref_file):
    trees = {}
    loop_id = 0

    with open(ref_file) as f:
        for line in f:
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
# 2. 批量 liftover
# =========================
def liftover_batch(regions, paf_file):

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

        target_chr = cols[0]
        target_start = int(cols[2])
        target_end = int(cols[3])

        idx = None
        for c in cols:
            if c.startswith("id:Z:"):
                parts = c.split(":")
                try:
                    idx = int(parts[-1])
                except:
                    continue
                break

        if idx is None:
            skipped_lines += 1
            continue

        mapping_dict[idx].append((target_chr, target_start, target_end))

    print(f"Liftover: {total_lines} lines parsed, {skipped_lines} lines skipped (no id:Z: field)")

    # 检查有多少原始 anchor 没有任何映射结果
    total_regions = len(regions)
    mapped_regions = len(mapping_dict)
    unmapped_regions = total_regions - mapped_regions
    print(f"Liftover: {mapped_regions}/{total_regions} anchors mapped, {unmapped_regions} anchors unmapped")

    return mapping_dict


# =========================
# 3. endpoint 匹配 loop
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
# 4. 主程序
# =========================
def main():
    print("Loading reference loops...")
    trees = load_reference_loops(ref_file)

    print("Reading query loops...")
    query_loops = []
    regions = []

    with open(query_file) as f:
        for line in f:
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
