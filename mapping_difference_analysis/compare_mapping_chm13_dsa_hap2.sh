#!/bin/bash
# written by Lingbin

# Compare retained read-pair IDs between CHM13 and donor-specific Hap2 mappings.
# Replace the input roots below with absolute paths without whitespace.
# Input layout: <root>/<sample>/<unit>/*.mapped.pairs.gz
# Inputs must be tab-delimited pairs: readID, chrom1, pos1, chrom2, pos2, ...
# Use the same source reads and the same input units for both mapping targets.
# Read IDs must be unique across all input units within each sample/target.
# Compress upstream plain .pairs files before running this script.
# Run Hap1 and Hap2 comparisons in separate new output directories.
# Existing sample output directories are refused rather than cleared.
# Adjust the {HG,NA}* sample pattern and SGE resources for your system.
# Requires Bash, awk, gzip/zcat and GNU sort/join; qsub submits each sample.
# file1 = CHM13; file2 = the selected donor-specific haplotype.
# Sequence in the retained output header means readID, not nucleotide sequence.
# Distance_difference = abs(cis_distance_file1 - cis_distance_file2), in bp.
# If either contact is trans, Distance_difference is NA; no coordinate liftover is performed.

set -u
work_dir="$PWD"
t2t_dir="/path/to/chm13_mapping"
hap_dir="/path/to/hap2_mapping"
cd "$work_dir" || exit 1
for dir in "$t2t_dir"/{HG,NA}*/; do
    [ -d "$dir" ] || continue
    sample=$(basename "$dir")
    echo "$sample"
    [ ! -e "$work_dir/$sample" ] && [ ! -L "$work_dir/$sample" ] || {
        echo "ERROR: output directory exists: $work_dir/$sample; use a new working directory." >&2; exit 1;
    }
    mkdir -p "$work_dir/$sample"
    cd "$work_dir/$sample" || { echo "cd $sample failed, skip"; continue; }
    {
        # Write sample names and input paths into the generated job script.
        printf '%s\n' '#!/bin/bash'
        printf '%s\n' '#$ -S /bin/bash'
        printf '%s\n' '#$ -cwd'
        printf '%s\n' '#$ -l mfree=8G,heavy_io=1'
        printf '%s\n' '#$ -pe serial 4'
        printf '%s\n' '#$ -o output.log'
        printf '%s\n' '#$ -e error.log'
        printf 'SAMPLE=%q\n'    "$sample"
        printf 'CHM13_DIR=%q\n' "$t2t_dir"
        printf 'HAP_DIR=%q\n'  "$hap_dir"
        # Keep the job body literal; variables below expand when the job runs.
        cat <<'EOF'
set -euo pipefail          # Stop on command or pipeline errors.
export LC_ALL=C
TMP="$(pwd)/tmp_$SAMPLE"
mkdir -p "$TMP"

prep () {                  # $1=input glob; $2=label
  local files=( $1 )
  if [ ! -e "${files[0]}" ]; then
    echo "ERROR: no input files match: $1" >&2; exit 1
  fi
  echo "  $2 inputs: ${#files[@]} file(s)" >&2   # Report the number of matched input files.
  zcat "${files[@]}" | awk -F'\t' '
    /^#/{next} NF==0{next}
    {if($2==$4)d=($5>$3)?$5-$3:$3-$5;else d="NA"; print $1"\t"d}' \
  | sort -k1,1 -S 4G --parallel=4 -T "$TMP"
}

prep "$CHM13_DIR/$SAMPLE/*/*.mapped.pairs.gz" CHM13 > "$TMP/f1.sorted"
prep "$HAP_DIR/$SAMPLE/*/*.mapped.pairs.gz" HAP   > "$TMP/f2.sorted"

# Require nonempty inputs sorted by read ID before joining.
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
