#!/bin/bash
# written by Lingbin

# Extract Hap1- and Hap2-specific pairs using prepared read-ID lists.
# Set mapping_dir1, mapping_dir2 and specific_dir to your absolute input paths.
# Mapping input: <mapping_dir>/<sample>/<unit>/<unit>.mapped.pairs.gz
# Pair records are tab-delimited, with readID in the first column.
# ID files: <specific_dir>/<sample>/<unit>/Hap1.specific.dedup.readid
#           <specific_dir>/<sample>/<unit>/Hap2.specific.dedup.readid
# Each ID list contains one readID per line, without a header.
# <specific_dir>/<sample>/<unit>/common.readid selects the units to process;
# its presence is used for traversal, while the two specific ID lists select pairs.
# Output per unit: <sample>.Hap1.unique.pairs and <sample>.Hap2.unique.pairs.
# Outputs retain all columns of the selected records and omit pairs header lines.
# Each sample/unit is processed separately in <working_directory>/<sample>/<unit>.
# Requires Bash, awk, pbgzip and SGE qsub; adjust module loads and SGE resources.
# Run in a new empty output directory: bash /path/to/extract_haplotype_specific_pairs.sh

set -u
work_dir="$PWD"
mapping_dir1="/path/to/hap1_mapping"
mapping_dir2="/path/to/hap2_mapping"
specific_dir="/path/to/haplotype_specific_readids"

for file in "$specific_dir"/*/*/common.readid; do
    [ -e "$file" ] || continue
    libdir=$(dirname "$file")               # .../<sample>/<unit>
    filename=$(basename "$libdir")          # Input unit name
    sample=$(basename "$(dirname "$libdir")")   # Sample name
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
