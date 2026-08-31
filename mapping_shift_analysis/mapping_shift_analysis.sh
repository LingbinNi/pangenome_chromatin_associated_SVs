#!/bin/bash
# written by Lingbin
set -euo pipefail

HAP_DIR="/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/mapping_shift_analysis_chm13_haps/mapping_shift_analysis_chm13_hap1"
WORKDIR="$HAP_DIR/Analysis"
JOB="statHap1"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

for dir in $(ls -d "$HAP_DIR"/{HG,NA}*)
do
    sample=$(basename "$dir")
    input="$dir/common_sequences.txt"
    if [ ! -s "$input" ]; then
        echo "WARNING: missing or empty: $input" >&2
        continue
    fi
    echo "$sample"
    mkdir -p "$sample"
    cd "$sample"
    rm -rf ./*

    # 头部用可展开 heredoc，把样本名和输入路径写死进去
    cat > script.sh << EOF
#!/bin/bash
#\$ -S /bin/bash
#\$ -cwd
#\$ -N $JOB
#\$ -l mfree=2G,heavy_io=1
#\$ -o output.log
#\$ -e error.log
set -euo pipefail
SAMPLE="$sample"
INPUT="$input"
EOF

    # 主体用带引号的 heredoc，awk 里的 \$4 才不会被 shell 吃掉
    cat >> script.sh << 'EOF'

awk -F'\t' -v sample="$SAMPLE" '
BEGIN { OFS = "\t" }
NR == 1 { next }
{
    total++
    # 第四列有有效数值，说明两种mapping中都为cis
    if ($4 != "NA" && $4 != "") {
        cis++
        d = $4 + 0
        if (d < 0) { d = -d }
        if (d > 5000)   { long5++   }
        if (d > 10000)  { long10++  }
        if (d > 25000)  { long25++  }
        if (d > 50000)  { long50++  }
        if (d > 100000) { long100++ }
    }
}
END { print sample, long5+0, long10+0, long25+0, long50+0, long100+0, cis+0, total+0 }
' "$INPUT" > stat.tsv.tmp

mv stat.tsv.tmp stat.tsv
wc -l "$INPUT"
echo "Done."
EOF

    qsub script.sh
    cd ..
done

# 所有样本跑完后自动合并
qsub -N "${JOB}_merge" -cwd -V -hold_jid "$JOB" \
     -l mfree=1G -o merge.log -e merge.err \
     "$WORKDIR/merge_stat.sh"

echo "Submitted. Final output: $WORKDIR/StatTotal.Hap1.tsv"
