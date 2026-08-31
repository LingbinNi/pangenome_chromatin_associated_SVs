#!/bin/bash
# written by Lingbin
set -u
work_dir="$PWD"
t2t_dir="/net/eichler/vol28/projects/hprc/nobackups/hic_hprc/01.ReferenceMapping-CHM13-HPRC-Y2-HPRC1"
hap_dir="/net/eichler/vol28/projects/hprc/nobackups/hic_hprc/04.SelfMapping-Hap2-HPRC-RagTag-Y2-HPRC1"
cd "$work_dir" || exit 1
for dir in "$t2t_dir"/{HG,NA}*/; do
    [ -d "$dir" ] || continue
    sample=$(basename "$dir")
    echo "$sample"
    mkdir -p "$work_dir/$sample"
    cd "$work_dir/$sample" || { echo "cd $sample failed, skip"; continue; }
    rm -rf ./*                       # 清理上次结果
    {
        # 这几行现在就展开(把样本名/路径写死进脚本)
        printf '%s\n' '#!/bin/bash'
        printf '%s\n' '#$ -S /bin/bash'
        printf '%s\n' '#$ -cwd'
        printf '%s\n' '#$ -l mfree=8G,heavy_io=1'
        printf '%s\n' '#$ -pe serial 4'
        printf '%s\n' '#$ -o output.log'
        printf '%s\n' '#$ -e error.log'
        printf 'SAMPLE=%q\n'    "$sample"
        printf 'CHM13_DIR=%q\n' "$t2t_dir"
        printf 'HAP1_DIR=%q\n'  "$hap_dir"
        # 下面 heredoc 是「原样」写入,awk/$SAMPLE/$TMP 都留到任务运行时才解析
        cat <<'EOF'
set -euo pipefail          # 注意加了 -e
export LC_ALL=C
TMP="$(pwd)/tmp_$SAMPLE"
mkdir -p "$TMP"

prep () {                  # $1=glob  $2=标签
  local files=( $1 )
  if [ ! -e "${files[0]}" ]; then
    echo "ERROR: no input files match: $1" >&2; exit 1
  fi
  echo "  $2 inputs: ${#files[@]} file(s)" >&2   # 记录 lane 数到 error.log
  zcat "${files[@]}" | awk -F'\t' '
    /^#/{next} NF==0{next}
    {if($2==$4)d=($5>$3)?$5-$3:$3-$5;else d="NA"; print $1"\t"d}' \
  | sort -k1,1 -S 4G --parallel=4 -T "$TMP"
}

prep "$CHM13_DIR/$SAMPLE/*/*.mapped.pairs.gz" CHM13 > "$TMP/f1.sorted"
prep "$HAP1_DIR/$SAMPLE/*/*.mapped.pairs.gz" HAP   > "$TMP/f2.sorted"

# 三重体检:非空 + 确实排好序(join 的前提)
for fx in f1 f2; do
  [ -s "$TMP/$fx.sorted" ] || { echo "ERROR: $fx.sorted empty" >&2; exit 1; }
  sort -c -k1,1 "$TMP/$fx.sorted" || { echo "ERROR: $fx.sorted not sorted" >&2; exit 1; }
done

join -t $'\t' "$TMP/f1.sorted" "$TMP/f2.sorted" \
  | awk -F'\t' 'BEGIN{OFS="\t"; print "Sequence","Distance_file1","Distance_file2","Distance_difference"}
      {d1=$2;d2=$3; if(d1=="NA"||d2=="NA")diff="NA";else diff=(d1>d2)?d1-d2:d2-d1; print $1,d1,d2,diff}' > common_sequences.txt
join -t $'\t' -v1 "$TMP/f1.sorted" "$TMP/f2.sorted" > unique_to_file1.txt
join -t $'\t' -v2 "$TMP/f1.sorted" "$TMP/f2.sorted" > unique_to_file2.txt
rm -rf "$TMP"

echo "Record counts:"
printf '  common_sequences.txt : %s\n' "$(( $(wc -l < common_sequences.txt) - 1 ))"
printf '  unique_to_file1.txt  : %s\n' "$(wc -l < unique_to_file1.txt)"
printf '  unique_to_file2.txt  : %s\n' "$(wc -l < unique_to_file2.txt)"
echo "Done."
EOF
    } > script.sh
    qsub script.sh
    cd "$work_dir" || exit 1
done
