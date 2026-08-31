#!/usr/bin/env python3
import argparse
import gzip
import csv

###############################
# 解析 SV ID
###############################
def parse_sv_id(sv_id):
    chrom, pos, svtype, svlen = sv_id.split("-")
    return chrom, int(pos), svtype, int(svlen)


def get_flanking_regions(sv_id, flank=500_000):
    chrom, pos, svtype, svlen = parse_sv_id(sv_id)
    if svtype == "DEL":
        sv_start = pos
        sv_end = pos + svlen
    elif svtype == "INS":
        sv_start = pos
        sv_end = pos + 1
    else:
        raise ValueError("Unsupported SV type: " + svtype)

    left = (chrom, max(1, sv_start - flank), sv_start - 1)
    right = (chrom, sv_end + 1, sv_end + flank)
    return left, right


###############################
# 读取 GTF
###############################
def read_gtf(gtf_path):
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            yield line.strip().split("\t")


def extract_genes_in_region(gtf_path, region):
    chrom, start, end = region
    genes = []

    for fields in read_gtf(gtf_path):
        if len(fields) < 9:
            continue
        c, source, feature, gstart, gend, score, strand, frame, attrs = fields
        if feature != "gene" or c != chrom:
            continue
        gstart = int(gstart)
        gend = int(gend)
        if gend < start or gstart > end:
            continue

        gene_id = ""
        gene_name = ""
        for item in attrs.split(";"):
            item = item.strip()
            if item.startswith("gene_id"):
                gene_id = item.split('"')[1]
            elif item.startswith("gene_name"):
                gene_name = item.split('"')[1]

        genes.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "start": gstart,
            "end": gend,
            "strand": strand
        })
    return genes


###############################
# 读取 TPM
###############################
def load_tpm(tpm_path):
    tpm = {}
    with open(tpm_path) as f:
        next(f)  # skip header
        for line in f:
            if not line.strip():
                continue
            feature, val = line.split()
            tpm[feature] = float(val)
    return tpm


def attach_tpm_to_genes(genes, tpm_dict):
    results = []
    for g in genes:
        key = g["gene_name"] if g["gene_name"] else g["gene_id"]
        g2 = g.copy()
        g2["TPM"] = tpm_dict.get(key, "NA")
        results.append(g2)
    return results


###############################
# 解析 VCF 获取原始 GT（按 ID 匹配）
###############################
def parse_vcf_for_sv(vcf_path, sv_id, sample_name):
    opener = gzip.open if vcf_path.endswith(".gz") else open

    with opener(vcf_path, "rt") as f:
        for line in f:
            # 读取 header 找到 sample index
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                if sample_name not in header:
                    raise ValueError(f"Sample {sample_name} not found in VCF")
                sample_index = header.index(sample_name)
                continue

            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            v_id = fields[2]  # VCF ID 列

            # 直接根据 ID 完全匹配
            if v_id != sv_id:
                continue

            # 提取 GT
            gt_field = fields[sample_index].split(":")[0]
            return gt_field, sample_name

    # 未找到
    return "not_found", sample_name


###############################
# 主函数
###############################
def main():
    parser = argparse.ArgumentParser(description="SV flanking genes with TPM and GT info")
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--tpm", required=True)
    parser.add_argument("--sv", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True, help="Output TSV file")
    args = parser.parse_args()

    # 解析 SV
    left, right = get_flanking_regions(args.sv)

    # 读取 TPM
    tpm = load_tpm(args.tpm)

    # 读取 VCF 获取原始 GT
    gt_field, sample_name = parse_vcf_for_sv(args.vcf, args.sv, args.sample)

    # 提取左右基因
    left_genes = attach_tpm_to_genes(extract_genes_in_region(args.gtf, left), tpm)
    right_genes = attach_tpm_to_genes(extract_genes_in_region(args.gtf, right), tpm)

    # 写入 TSV
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "gene_id", "gene_name", "start", "end", "strand",
            "TPM", "haplotype", "sample"
        ], delimiter="\t")
        writer.writeheader()

        for g in left_genes + right_genes:
            g_out = g.copy()
            g_out["haplotype"] = gt_field
            g_out["sample"] = sample_name
            writer.writerow(g_out)

    print(f"Output written to {args.out}")


if __name__ == "__main__":
    main()
