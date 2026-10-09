#!/usr/bin/env python3
# Aggregate local contact matrices by SV genotype and plot group comparisons.
# Requires Python 3, cooler (analysis version: 0.9.3), numpy, pandas and matplotlib.
# Edit SV_ID_FILE, GENOTYPE_FILE and COOL_DIRS before running.
# Run: python /path/to/aggregate_contact_matrices_by_sv_genotype.py
# OUT_BASE defaults to the current working directory; one subdirectory per SV.
# SV list: no header, one chromosome-start-type-length ID per line, with a
# 0-based start and a positive length in bp; interval end = start + length.
# Genotype TSV (.gz supported): sample, sv_id, hap1_carrier, hap2_carrier.
# Carrier columns contain 0, 1 or NA. Their sum defines dosage 0/1/2,
# corresponding to no_carrier (0/0), single_carrier (0/1), both_carrier (1/1).
# Exclude samples with missing dosage or no genotype record for the SV.
# Keep the first record for each sample/SV and the first cool file per sample.
# Cool input: <root>/<sample>/Merge/*.25kb.corrected.cool by default.
# Sample directory names must match the genotype table. Supply contact maps
# with matching CHM13 coordinates, chromosome names and 25-kb bin definitions.
# USE_CORRECTED selects the filename suffix; matrices are read from stored
# count values with balance=False, without applying Cooler balancing weights.
# Extract [max(0, SV_start - 500000), SV_start + length + 500000).
# Divide each individual matrix by its total signal, then average within each
# genotype group. Retain groups with at least MIN_GROUP_N valid matrices.
# Top row: group mean contact density with a shared LogNorm color scale.
# Use the 25th and 99th percentiles of positive finite group-mean values as
# color limits and the white-yellow-red-black palette from Integration_TAD.heatmap.2.py.
# Bottom row: log2((group_mean + eps)/(reference_mean + eps)), eps = vmax / 1000.
# Plot 0/1 vs 0/0, 1/1 vs 0/0 and 1/1 vs 0/1, in that order, when both groups exist.
# With all three groups retained, the figure has two rows and three columns.
# Both rows use positions relative to the SV start; dashed lines mark start/end.
# Outputs per SV: <sv_id>_group_mean_matrices.npz (no_carrier, single_carrier,
# both_carrier when present, plus region), <sv_id>_aggregate_hic_bygroup.pdf
# and <sv_id>_aggregate_hic_bygroup.png. NPZ keys follow the density source script.

import glob
import os
import sys

import numpy as np
import pandas as pd
import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
CMAP_HIC = LinearSegmentedColormap.from_list(
    "fall", ["white", "#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c", "#800026", "black"])

## ------------------ Settings ------------------
SV_ID_FILE = "/path/to/sv_ids.txt"
GENOTYPE_FILE = "/path/to/SV_CHM13_all_genotypes.tsv.gz"
OUT_BASE = "."
FLANK_BP = 500_000
RESOLUTION = "25kb"          # Filename resolution; supply matching bin sizes
USE_CORRECTED = True         # Select *.corrected.cool when True

COOL_DIRS = [
    "/path/to/cohort1/contact_matrices",
    "/path/to/cohort2/contact_matrices",
]

GROUP_LEVELS = ["no_carrier", "single_carrier", "both_carrier"]
MIN_GROUP_N = 2              # Minimum valid matrices per plotted group
VMAX_PCT = 99                # Shared upper percentile of positive group-mean values
VMIN_PCT = 25                # Shared lower percentile of positive group-mean values
DIFF_LIM = 1.0               # Bottom-row log2 ratio color limits

## ------------------ 1. Parse SV IDs ------------------
def parse_sv_id(sv_id):
    chrom, start, svtype, svlen = sv_id.split("-")
    start, svlen = int(start), int(svlen)
    return chrom, start, start + svlen, svtype, svlen

## ------------------ 2. Find sample cool files ------------------
def find_sample_cools():
    suffix = f".{RESOLUTION}.corrected.cool" if USE_CORRECTED else f".{RESOLUTION}.cool"
    sample_cools = {}
    dups = []
    for d in COOL_DIRS:
        for f in sorted(glob.glob(os.path.join(d, "*", "Merge", f"*{suffix}"))):
            sample = os.path.basename(os.path.dirname(os.path.dirname(f)))
            if sample in sample_cools:
                dups.append(sample)
                continue  # Keep the first file for each sample
            sample_cools[sample] = f
    if dups:
        print(f"Multiple files found for these samples; keeping the first: {', '.join(sorted(set(dups)))}")
    return sample_cools

sample_cools = find_sample_cools()
print(f"Found {RESOLUTION} cool files for {len(sample_cools)} samples")
all_samples = set(sample_cools)

## ------------------ 3. Read genotypes and assign groups ------------------
print("Reading genotypes...")
geno = pd.read_csv(GENOTYPE_FILE, sep="\t",
                   usecols=["sample", "sv_id", "hap1_carrier", "hap2_carrier"])
geno["hap1_carrier"] = pd.to_numeric(geno["hap1_carrier"], errors="coerce")
geno["hap2_carrier"] = pd.to_numeric(geno["hap2_carrier"], errors="coerce")
geno["dosage"] = geno["hap1_carrier"] + geno["hap2_carrier"]

def group_samples_for_sv(sv_id):
    """Return grouped samples with nonmissing genotypes for this SV."""
    sub = geno.loc[geno["sv_id"] == sv_id, ["sample", "dosage"]].drop_duplicates("sample")
    sub = sub[sub["sample"].isin(all_samples)].dropna(subset=["dosage"])
    n_excluded = len(all_samples) - len(sub)
    if n_excluded > 0:
        print(f"  Excluded {n_excluded} samples without genotype information for this SV")
    groups = {g: [] for g in GROUP_LEVELS}
    for _, row in sub.iterrows():
        d = int(row["dosage"])
        g = GROUP_LEVELS[min(d, 2)]
        groups[g].append(row["sample"])
    return groups

## ------------------ 4. Extract and normalize individual matrices ------------------
def fetch_norm_matrix(cool_path, region):
    """Fetch the local contact matrix and divide by its total signal."""
    clr = cooler.Cooler(cool_path)
    m = clr.matrix(balance=False).fetch(region).astype(float)
    # Read stored count values for either selected filename suffix.
    total = np.nansum(m)
    if not np.isfinite(total) or total <= 0:
        return None
    return m / total

## ------------------ 5. Process each SV ------------------
with open(SV_ID_FILE) as fh:
    sv_ids = [l.strip() for l in fh if l.strip()]

for sv_id in sv_ids:
    print(f"===> Processing {sv_id}")
    chrom, sv_start, sv_end, svtype, svlen = parse_sv_id(sv_id)
    region_start = max(0, sv_start - FLANK_BP)
    region_end = sv_end + FLANK_BP
    region = f"{chrom}:{region_start}-{region_end}"

    sv_dir = os.path.join(OUT_BASE, sv_id)
    os.makedirs(sv_dir, exist_ok=True)

    groups = group_samples_for_sv(sv_id)
    print("Genotyped samples per group: " + ", ".join(f"{g}={len(s)}" for g, s in groups.items() if s))

    ## --- Average valid matrices within each genotype group ---
    group_mats = {}
    for g, samples in groups.items():
        if len(samples) < MIN_GROUP_N:
            if samples:
                print(f"  Skipping {g}: {len(samples)} samples (<{MIN_GROUP_N})")
            continue
        mats, shape0 = [], None
        for s in samples:
            try:
                m = fetch_norm_matrix(sample_cools[s], region)
            except Exception as e:
                print(f"  [WARNING] Skipping {s}: matrix extraction failed: {e}")
                continue
            if m is None:
                print(f"  [WARNING] Skipping {s}: total matrix signal is nonpositive or nonfinite")
                continue
            if shape0 is None:
                shape0 = m.shape
            if m.shape != shape0:
                print(f"  [WARNING] Skipping {s}: matrix shape {m.shape} differs from {shape0}")
                continue
            mats.append(m)
        if len(mats) >= MIN_GROUP_N:
            group_mats[g] = np.nanmean(np.stack(mats), axis=0)
            print(f"  {g}: {len(mats)} valid matrices used in the group mean")

    if not group_mats:
        print("  No group meets the minimum valid matrix count; skipping this SV")
        continue

    ## --- Save group mean matrices ---
    np.savez_compressed(os.path.join(sv_dir, f"{sv_id}_group_mean_matrices.npz"),
                        **group_mats, region=np.array([region]))

    ## --- Plot ---
    present = [g for g in GROUP_LEVELS if g in group_mats]
    n_bins = group_mats[present[0]].shape[0]
    # Shared LogNorm limits from positive finite values across all group means.
    allvals = np.concatenate([m[np.isfinite(m) & (m > 0)].ravel()
                              for m in group_mats.values()])
    vmax = np.percentile(allvals, VMAX_PCT)
    vmin = np.percentile(allvals, VMIN_PCT)

    diff_pairs = [
        ("single_carrier", "no_carrier"),  # 0/1 vs 0/0
        ("both_carrier", "no_carrier"),    # 1/1 vs 0/0
        ("both_carrier", "single_carrier"),  # 1/1 vs 0/1
    ]
    diff_pairs = [(g, ref) for g, ref in diff_pairs
                  if g in group_mats and ref in group_mats]
    ncol = max(len(present), len(diff_pairs)) if diff_pairs else len(present)
    nrow = 2 if diff_pairs else 1
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4.4 * nrow),
                             squeeze=False)

    extent_kb = [(region_start - sv_start) / 1000, (region_end - sv_start) / 1000]
    sv_kb = [0, (sv_end - sv_start) / 1000]

    def decorate(ax, title):
        for v in sv_kb:
            ax.axvline(v, color="black", ls="--", lw=0.6)
            ax.axhline(v, color="black", ls="--", lw=0.6)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Distance from SV start (kb)", fontsize=8)
        ax.tick_params(labelsize=7)

    # Top row: group mean contact density on the shared logarithmic color scale.
    for j, g in enumerate(present):
        ax = axes[0][j]
        im = ax.imshow(group_mats[g], cmap=CMAP_HIC,
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       extent=[extent_kb[0], extent_kb[1], extent_kb[1], extent_kb[0]],
                       interpolation="none")
        # Title n counts genotyped samples; the log reports valid matrices used.
        decorate(ax, f"{g} (n={len(groups[g])})")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    for j in range(len(present), ncol):
        axes[0][j].axis("off")

    # Bottom row: pairwise log2 ratios of normalized group mean matrices.
    if diff_pairs:
        for j, (g, ref) in enumerate(diff_pairs):
            ax = axes[1][j]
            eps = vmax / 1000  # Preserve the pseudocount used in the density version.
            ratio = np.log2((group_mats[g] + eps) / (group_mats[ref] + eps))
            im = ax.imshow(ratio, cmap="bwr", vmin=-DIFF_LIM, vmax=DIFF_LIM,
                           extent=[extent_kb[0], extent_kb[1], extent_kb[1], extent_kb[0]],
                           interpolation="none")
            decorate(ax, f"log2( {g} / {ref} )")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        for j in range(len(diff_pairs), ncol):
            axes[1][j].axis("off")

    fig.suptitle(f"{sv_id}  ({RESOLUTION}, {'corrected' if USE_CORRECTED else 'raw'})",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_pdf = os.path.join(sv_dir, f"{sv_id}_aggregate_hic_bygroup.pdf")
    fig.savefig(out_pdf, dpi=300)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=300)
    plt.close(fig)
    print(f"  Written: {out_pdf}")

print("Done. Results are in the subdirectory for each SV.")
