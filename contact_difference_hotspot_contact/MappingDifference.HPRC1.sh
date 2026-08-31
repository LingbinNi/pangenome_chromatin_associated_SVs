#!/bin/bash
# written by Lingbin
set -u
work_dir="$PWD"
mapping_dir1="/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/03.SelfMapping-Hap1-HPRC-RagTag-Y2-HPRC1-RawData"
mapping_dir2="/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/04.SelfMapping-Hap2-HPRC-RagTag-Y2-HPRC1-RawData"
specific_dir="/net/eichler/vol28/projects/hic_cohorts/nobackups/00.Analysis/mapping_difference_analysis_hprc1"

for file in "$specific_dir"/*/*/common.readid; do
    [ -e "$file" ] || continue
    libdir=$(dirname "$file")               # .../sample/filename
    filename=$(basename "$libdir")          # 文库名
    sample=$(basename "$(dirname "$libdir")")   # 样本名
    echo "$sample-$filename"

    mkdir -p "$work_dir/$sample/$filename"
    cd "$work_dir/$sample/$filename" || { echo "cd failed, skip"; continue; }

    {
        printf '%s\n' '#!/bin/bash'
        printf '%s\n' '#$ -S /bin/bash'
        printf '%s\n' '#$ -cwd'
        printf '%s\n' '#$ -l mfree=4G,heavy_io=1'
        printf '%s\n' '#$ -pe serial 4'
        printf '%s\n' '#$ -o output.log'
        printf '%s\n' '#$ -e error.log'
        printf 'SAMPLE=%q\n'   "$sample"
        printf 'FILENAME=%q\n' "$filename"
        printf 'MAP1=%q\n'     "$mapping_dir1"
        printf 'MAP2=%q\n'     "$mapping_dir2"
        printf 'SPEC=%q\n'     "$specific_dir"
        cat <<'EOF'
module load pairtools/1.0.3
export LC_ALL=C

# ---- Hap1 ----
pbgzip -dc -n 4 "$MAP1/$SAMPLE/$FILENAME/$FILENAME.mapped.pairs.gz" \
  | awk -F'\t' 'NR==FNR{ids[$1]=1; next} !/^#/ && ($1 in ids)' \
      "$SPEC/$SAMPLE/$FILENAME/Hap1.specific.dedup.readid" - \
  > "$SAMPLE.Hap1.unique.pairs"
echo "Hap1 matched: $(wc -l < "$SAMPLE.Hap1.unique.pairs") lines"

# ---- Hap2 ----
pbgzip -dc -n 4 "$MAP2/$SAMPLE/$FILENAME/$FILENAME.mapped.pairs.gz" \
  | awk -F'\t' 'NR==FNR{ids[$1]=1; next} !/^#/ && ($1 in ids)' \
      "$SPEC/$SAMPLE/$FILENAME/Hap2.specific.dedup.readid" - \
  > "$SAMPLE.Hap2.unique.pairs"
echo "Hap2 matched: $(wc -l < "$SAMPLE.Hap2.unique.pairs") lines"
EOF
    } > script.sh

    qsub script.sh
    cd "$work_dir" || exit 1
done
