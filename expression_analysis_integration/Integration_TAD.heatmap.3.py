#!/usr/bin/env python3
# ============================================================
# 功能：对每个 SV，按基因型（no/single/both carrier）把样本分三组，
#       从各样本 25kb corrected .cool 中取 [SV ± 500kb] 的互作子矩阵，
#       每个样本按区域总信号归一化后求组内均值，画：
#         第一行：三组 log2(Observed/Expected) 热图
#                 （按对角线距离归一化，消除距离衰减背景，
#                   红=高于期望互作，蓝=低于期望；TAD呈红色方块、
#                   边界呈蓝色隔断，结构差异一目了然）
#         第二行：log2(single/no) 和 log2(both/no) 差异热图
#       分组逻辑与 Integration_TAD_v2.R 完全一致：
#         dosage 0/1/2 分三组，基因型缺失(NA)的样本剔除
#
# 依赖：pip install cooler numpy pandas matplotlib
# 用法：python Integration_HiC_heatmap.py
# ============================================================

import glob
import os
import sys

import numpy as np
import pandas as pd
import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

## ------------------ 参数（与 R 脚本对齐） ------------------
SV_ID_FILE = "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/expression_analysis_integration/01.Integration_TAD/example"
GENOTYPE_FILE = "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/structure_associated_SV_identification/00.SV_CHM13_ALL/SV_CHM13_all_genotypes.tsv.gz"
OUT_BASE = "."
FLANK_BP = 500_000
RESOLUTION = "25kb"          # 可改 "10kb" / "100kb"
USE_CORRECTED = True         # True 用 *.corrected.cool（ICE校正后）

COOL_DIRS = [
    "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/05.ContactMatrix-CHM13-HPRC-Y2-HPRC1",
    "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/05.ContactMatrix-CHM13-HPRC-Y2-HPRC2",
]

GROUP_LEVELS = ["no_carrier", "single_carrier", "both_carrier"]
MIN_GROUP_N = 2              # 组内样本数不足则该组不画（与 R 的 N>=2 一致）
OE_LIM = 0.8                 # 第一行 log2(O/E) 色标范围 [-OE_LIM, +OE_LIM]
DIFF_LIM = 1.0               # 第二行 log2 比值色标范围 [-DIFF_LIM, +DIFF_LIM]
EPS = 1e-6                   # 差异图伪计数，避免除零（矩阵已按总量归一化，量级约1e-4）

## ------------------ 1. 解析 SV id ------------------
def parse_sv_id(sv_id):
    chrom, start, svtype, svlen = sv_id.split("-")
    start, svlen = int(start), int(svlen)
    return chrom, start, start + svlen, svtype, svlen

## ------------------ 2. 找到所有样本的 cool 文件 ------------------
def find_sample_cools():
    suffix = f".{RESOLUTION}.corrected.cool" if USE_CORRECTED else f".{RESOLUTION}.cool"
    sample_cools = {}
    dups = []
    for d in COOL_DIRS:
        for f in sorted(glob.glob(os.path.join(d, "*", "Merge", f"*{suffix}"))):
            sample = os.path.basename(os.path.dirname(os.path.dirname(f)))
            if sample in sample_cools:
                dups.append(sample)
                continue  # 保留第一个（与 R 脚本一致）
            sample_cools[sample] = f
    if dups:
        print(f"以下样本在多个目录下重复出现，已保留第一个: {', '.join(sorted(set(dups)))}")
    return sample_cools

sample_cools = find_sample_cools()
print(f"共找到 {len(sample_cools)} 个样本的 {RESOLUTION} cool 文件")
all_samples = set(sample_cools)

## ------------------ 3. 读取完整基因型表，构建分组 ------------------
print("正在读取基因型文件...")
geno = pd.read_csv(GENOTYPE_FILE, sep="\t",
                   usecols=["sample", "sv_id", "hap1_carrier", "hap2_carrier"])
geno["hap1_carrier"] = pd.to_numeric(geno["hap1_carrier"], errors="coerce")
geno["hap2_carrier"] = pd.to_numeric(geno["hap2_carrier"], errors="coerce")
geno["dosage"] = geno["hap1_carrier"] + geno["hap2_carrier"]

def group_samples_for_sv(sv_id):
    """返回 {group_name: [samples]}，剔除基因型缺失/表中无记录的样本"""
    sub = geno.loc[geno["sv_id"] == sv_id, ["sample", "dosage"]].drop_duplicates("sample")
    sub = sub[sub["sample"].isin(all_samples)].dropna(subset=["dosage"])
    n_excluded = len(all_samples) - len(sub)
    if n_excluded > 0:
        print(f"  {n_excluded} 个样本无该位点基因型信息或不在基因型表中，已剔除")
    groups = {g: [] for g in GROUP_LEVELS}
    for _, row in sub.iterrows():
        d = int(row["dosage"])
        g = GROUP_LEVELS[min(d, 2)]
        groups[g].append(row["sample"])
    return groups

## ------------------ 4. 取单样本子矩阵并归一化 ------------------
def fetch_norm_matrix(cool_path, region):
    """取 region 的互作子矩阵；除以区域内总信号，消除测序深度差异"""
    clr = cooler.Cooler(cool_path)
    m = clr.matrix(balance=False).fetch(region).astype(float)
    # corrected.cool 已做 ICE 校正，直接用 count 值；未校正数据同样按总量归一化
    total = np.nansum(m)
    if not np.isfinite(total) or total <= 0:
        return None
    return m / total

## ------------------ 5. Observed/Expected 归一化 ------------------
def obs_over_exp(m):
    """按对角线距离归一化：每条对角线除以其均值，消除距离衰减背景"""
    n = m.shape[0]
    oe = np.full_like(m, np.nan)
    for d in range(n):
        diag_vals = np.diagonal(m, offset=d)
        exp = np.nanmean(diag_vals)
        if np.isfinite(exp) and exp > 0:
            idx = np.arange(n - d)
            oe[idx, idx + d] = m[idx, idx + d] / exp
            oe[idx + d, idx] = oe[idx, idx + d]   # 对称补全下三角
    return oe

## ------------------ 6. 主循环 ------------------
with open(SV_ID_FILE) as fh:
    sv_ids = [l.strip() for l in fh if l.strip()]

for sv_id in sv_ids:
    print(f"===> 处理 {sv_id}")
    chrom, sv_start, sv_end, svtype, svlen = parse_sv_id(sv_id)
    region_start = max(0, sv_start - FLANK_BP)
    region_end = sv_end + FLANK_BP
    region = f"{chrom}:{region_start}-{region_end}"

    sv_dir = os.path.join(OUT_BASE, sv_id)
    os.makedirs(sv_dir, exist_ok=True)

    groups = group_samples_for_sv(sv_id)
    print("分组样本数: " + ", ".join(f"{g}={len(s)}" for g, s in groups.items() if s))

    ## --- 逐组求均值矩阵 ---
    group_mats = {}
    for g, samples in groups.items():
        if len(samples) < MIN_GROUP_N:
            if samples:
                print(f"  组 {g} 仅 {len(samples)} 个样本(<{MIN_GROUP_N})，跳过该组")
            continue
        mats, shape0 = [], None
        for s in samples:
            try:
                m = fetch_norm_matrix(sample_cools[s], region)
            except Exception as e:
                print(f"  [警告] 样本 {s} 取矩阵失败，已跳过: {e}")
                continue
            if m is None:
                print(f"  [警告] 样本 {s} 区域总信号为 0，已跳过")
                continue
            if shape0 is None:
                shape0 = m.shape
            if m.shape != shape0:
                print(f"  [警告] 样本 {s} 矩阵尺寸 {m.shape} 与 {shape0} 不一致，已跳过")
                continue
            mats.append(m)
        if len(mats) >= MIN_GROUP_N:
            group_mats[g] = np.nanmean(np.stack(mats), axis=0)
            print(f"  组 {g}: 实际用于集成的样本数 = {len(mats)}")

    if not group_mats:
        print("  没有任何组满足条件，跳过该 SV")
        continue

    ## --- 计算各组 O/E 矩阵，连同均值矩阵一起保存，便于复用 ---
    group_oe = {g: obs_over_exp(m) for g, m in group_mats.items()}
    np.savez_compressed(
        os.path.join(sv_dir, f"{sv_id}_group_mean_matrices.npz"),
        **{f"mean_{g}": m for g, m in group_mats.items()},
        **{f"oe_{g}": m for g, m in group_oe.items()},
        region=np.array([region]))

    ## --- 画图 ---
    present = [g for g in GROUP_LEVELS if g in group_mats]
    has_diff = "no_carrier" in group_mats and len(present) > 1
    diff_pairs = [(g, "no_carrier") for g in ("single_carrier", "both_carrier")
                  if has_diff and g in group_mats]
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

    # 第一行：三组 log2(O/E) 热图（红=高于期望互作，蓝=低于期望）
    for j, g in enumerate(present):
        ax = axes[0][j]
        with np.errstate(divide="ignore", invalid="ignore"):
            oe_log = np.log2(group_oe[g])
        im = ax.imshow(oe_log, cmap="RdBu_r", vmin=-OE_LIM, vmax=OE_LIM,
                       extent=[extent_kb[0], extent_kb[1], extent_kb[1], extent_kb[0]],
                       interpolation="none")
        decorate(ax, f"{g} (n={len(groups[g])})  log2(O/E)")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    for j in range(len(present), ncol):
        axes[0][j].axis("off")

    # 第二行：log2 比值差异热图（携带者 vs 非携带者）
    if diff_pairs:
        for j, (g, ref) in enumerate(diff_pairs):
            ax = axes[1][j]
            ratio = np.log2((group_mats[g] + EPS) / (group_mats[ref] + EPS))
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
    print(f"  已输出: {out_pdf}")

print("全部完成，结果在各 SV 对应的子文件夹中。")