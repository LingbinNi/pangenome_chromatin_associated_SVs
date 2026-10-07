#!/bin/bash
# written by Lingbin

# Summarize absolute contact-distance differences between CHM13 and Hap2.
# Set HAP_DIR to the output root of compare_mapping_chm13_dsa_hap2.sh.
# Input layout: <HAP_DIR>/<sample>/common_sequences.txt (tab-delimited, with header).
# Wait for ALL comparison jobs to finish successfully before running this script.
# Use absolute paths without whitespace or shell-special characters.
# Adjust the {HG,NA}* sample pattern, JOB and SGE resources for your system.
# Requires Bash, awk and SGE qsub. Run with: bash summarize_mapping_distance_differences_hap2.sh
# The fourth input column must contain a numeric distance difference in bp or NA.
# Thresholds are strictly > 5/10/25/50/100 kb and counts are cumulative.
# stat.tsv has no header: sample, gt5kb, gt10kb, gt25kb, gt50kb, gt100kb, both_cis, common_total.
# both_cis counts shared pairs that are cis in both references; common_total includes trans.
# These are counts, not percentages; use both_cis as denominator for fractions among shared cis pairs.
# A generated merge job combines the sample rows into StatTotal.Hap2.tsv with a header.
# WORKDIR must be new; existing results are never cleared automatically.

set -euo pipefail

HAP_DIR="/path/to/chm13_hap2_comparison"
WORKDIR="$HAP_DIR/Analysis"
JOB="statHap2"

[ -d "$HAP_DIR" ] || { echo "ERROR: input root does not exist: $HAP_DIR" >&2; exit 1; }
[ ! -e "$WORKDIR" ] && [ ! -L "$WORKDIR" ] || {
    echo "ERROR: output directory exists: $WORKDIR; set WORKDIR to a new directory." >&2; exit 1;
}
mkdir -p "$WORKDIR"
cd "$WORKDIR"

submitted=0
for dir in "$HAP_DIR"/{HG,NA}*/
do
    [ -d "$dir" ] || continue
    sample=$(basename "$dir")
    input="$dir/common_sequences.txt"
    if [ ! -s "$input" ]; then
        echo "WARNING: missing or empty: $input" >&2
        continue
    fi
    echo "$sample"
    mkdir -p "$sample"
    cd "$sample"

    # Expand sample names and input paths into the job header.
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

    # Keep the job body literal so awk fields are evaluated when the job runs.
    cat >> script.sh << 'EOF'

awk -F'\t' -v sample="$SAMPLE" '
BEGIN { OFS = "\t" }
NR == 1 { next }
{
    total++
    # A numeric fourth column identifies shared pairs that are cis in both mappings.
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
    submitted=$((submitted + 1))
    cd ..
done

[ "$submitted" -gt 0 ] || { echo "ERROR: no usable sample inputs found." >&2; exit 1; }

# Generate the merge helper so no separate user-provided script is required.
cat > "$WORKDIR/merge_stat.sh" << 'EOF'
#!/bin/bash
#$ -S /bin/bash
set -euo pipefail
export LC_ALL=C
printf 'Sample\tDifference_gt_5kb\tDifference_gt_10kb\tDifference_gt_25kb\tDifference_gt_50kb\tDifference_gt_100kb\tCommon_cis_pairs\tCommon_pairs_total\n' > StatTotal.Hap2.tsv.tmp
for dir in ./*/; do
    [ -d "$dir" ] || continue
    [ -s "$dir/stat.tsv" ] || {
        echo "ERROR: missing or empty: $dir/stat.tsv; check the sample job logs." >&2; exit 1;
    }
    cat "$dir/stat.tsv" >> StatTotal.Hap2.tsv.tmp
done
mv StatTotal.Hap2.tsv.tmp StatTotal.Hap2.tsv
EOF

# Submit the merge after all submitted sample jobs have finished.
qsub -N "${JOB}_merge" -cwd -V -hold_jid "$JOB" \
     -l mfree=1G -o merge.log -e merge.err \
     "$WORKDIR/merge_stat.sh"

echo "Submitted. Final output: $WORKDIR/StatTotal.Hap2.tsv"
