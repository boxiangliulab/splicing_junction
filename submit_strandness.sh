#!/bin/bash
#PBS -P 11003054
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=02:00:00
#PBS -N eQTLGen_splicing_strandness
#PBS -j oe
#PBS -o prepare/strandness.log

# Please include make sure squashfs can be used
#module load squashfs/4.4

cd "$PBS_O_WORKDIR"
mkdir -p prepare

# ----------------------------
# Paths
# ----------------------------
ROOT=$PWD/..
SIF=${ROOT}/singularity_img/eqtlgen_splicing_junction_strandness.sif
BED=${ROOT}/data/gencode.v48.annotation.bed12

# If you run the star alignment pipeline, use the default output folder 
BAM_DIR=${BAM_OUTPUT_DIR}/aligned_bam
# Otherwise, specify your BAM file directory
#BAM_DIR=[your BAM file directory]


# ----------------------------
# Select representative BAMs
# ----------------------------
# Sample 5 BAMs: first, last, middle, and 2 random
# This provides better coverage across the cohort and detects batch effects

all_bams=($(ls "${BAM_DIR}"/*.bam | sort))
n_bams=${#all_bams[@]}

if [ $n_bams -eq 0 ]; then
  echo "[ERROR] No BAM files found in ${BAM_DIR}"
  exit 1
fi

echo "[INFO] Found ${n_bams} BAM files in total"

# Select BAMs to test
bams_to_test=()
bam_labels=()

# Always test first and last
bams_to_test+=("${all_bams[0]}")
bam_labels+=("FIRST")

if [ $n_bams -gt 1 ]; then
  bams_to_test+=("${all_bams[$((n_bams-1))]}")
  bam_labels+=("LAST")
fi

# Add middle if we have at least 3 BAMs
if [ $n_bams -ge 3 ]; then
  mid_idx=$((n_bams / 2))
  bams_to_test+=("${all_bams[$mid_idx]}")
  bam_labels+=("MIDDLE")
fi

# Add 2 random BAMs if we have at least 5 BAMs
if [ $n_bams -ge 5 ]; then
  # Generate 2 random indices (avoiding first, last, middle)
  for i in {1..2}; do
    # Try up to 20 times to find a unique random BAM
    for attempt in {1..20}; do
      rand_idx=$((RANDOM % n_bams))
      rand_bam="${all_bams[$rand_idx]}"
      # Check if this BAM is already in our test list
      is_duplicate=0
      for test_bam in "${bams_to_test[@]}"; do
        if [ "$test_bam" = "$rand_bam" ]; then
          is_duplicate=1
          break
        fi
      done
      if [ $is_duplicate -eq 0 ]; then
        bams_to_test+=("$rand_bam")
        bam_labels+=("RANDOM${i}")
        break
      fi
    done
  done
fi

echo "[INFO] Testing ${#bams_to_test[@]} representative BAM files for strandness"
echo ""

# ----------------------------
# Run infer_experiment.py
# ----------------------------
out="prepare/strandness.txt"
raw_out="prepare/strandness_raw.txt"

echo "======================================================" > $out
echo "  Strandness Detection Report" >> $out
echo "  Generated: $(date)" >> $out
echo "======================================================" >> $out
echo "" >> $out

# Store results for consistency check
declare -a sense_fractions
declare -a antisense_fractions

for i in "${!bams_to_test[@]}"; do
  bam="${bams_to_test[$i]}"
  label="${bam_labels[$i]}"
  bam_dir=$(dirname "${bam}")

  echo "[INFO] Testing BAM ${label}: $(basename ${bam})"

  echo "------------------------------------------------------" >> $out
  echo "Sample ${label}: $(basename ${bam})" >> $out
  echo "------------------------------------------------------" >> $out

  # Run infer_experiment.py
  tmp_out=$(mktemp)
  singularity exec --userns --containall \
    -B "${ROOT}":/work \
    -B "${bam_dir}:${bam_dir}" \
    --pwd /work \
    "${SIF}" \
    infer_experiment.py \
    -s 1000000 \
    -r /work/data/$(basename "${BED}") \
    -i "${bam}" \
    > ${tmp_out} 2>&1

  cat ${tmp_out} >> ${raw_out}

  # Parse the output to extract fractions
  # For paired-end data, the format is:
  # Fraction of reads explained by "1++,1--,2+-,2-+": 0.4322  (sense/forward)
  # Fraction of reads explained by "1+-,1-+,2++,2--": 0.4293  (antisense/reverse)
  sense=$(grep -oP 'Fraction of reads explained by "1\+\+,1--,2\+-,2-\+":\s+\K[0-9.]+' ${tmp_out} || echo "N/A")
  antisense=$(grep -oP 'Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s+\K[0-9.]+' ${tmp_out} || echo "N/A")

  cat ${tmp_out} >> $out
  echo "" >> $out

  if [ "$sense" != "N/A" ] && [ "$antisense" != "N/A" ]; then
    sense_fractions+=("$sense")
    antisense_fractions+=("$antisense")
  fi

  rm -f ${tmp_out}
done

echo "" >> $out
echo "======================================================" >> $out
echo "  SUMMARY AND RECOMMENDATION" >> $out
echo "======================================================" >> $out
echo "" >> $out

# ----------------------------
# Analyze results and recommend strand
# ----------------------------
if [ ${#sense_fractions[@]} -eq 0 ]; then
  echo "[ERROR] Could not parse strandness results. Please check ${raw_out}" >> $out
  echo "[ERROR] Could not determine strandness. Check prepare/strandness.txt"
  exit 1
fi

# Calculate average sense and antisense fractions
sum_sense=0
sum_antisense=0
for val in "${sense_fractions[@]}"; do
  sum_sense=$(awk "BEGIN {print $sum_sense + $val}")
done
for val in "${antisense_fractions[@]}"; do
  sum_antisense=$(awk "BEGIN {print $sum_antisense + $val}")
done

n_samples=${#sense_fractions[@]}
avg_sense=$(awk "BEGIN {printf \"%.4f\", $sum_sense / $n_samples}")
avg_antisense=$(awk "BEGIN {printf \"%.4f\", $sum_antisense / $n_samples}")

echo "Tested ${n_samples} BAM files:" >> $out
for i in "${!bam_labels[@]}"; do
  echo "  - ${bam_labels[$i]}: $(basename ${bams_to_test[$i]})" >> $out
done
echo "" >> $out

echo "Average fraction of reads:" >> $out
echo "  Sense/Forward (1++,1--,2+-,2-+):     ${avg_sense}" >> $out
echo "  Antisense/Reverse (1+-,1-+,2++,2--): ${avg_antisense}" >> $out
echo "" >> $out

# Check consistency (coefficient of variation < 20%)
max_diff_sense=0
max_diff_antisense=0
for val in "${sense_fractions[@]}"; do
  diff=$(awk "BEGIN {d = $val - $avg_sense; if (d < 0) d = -d; print d}")
  max_diff_sense=$(awk "BEGIN {print ($diff > $max_diff_sense) ? $diff : $max_diff_sense}")
done
for val in "${antisense_fractions[@]}"; do
  diff=$(awk "BEGIN {d = $val - $avg_antisense; if (d < 0) d = -d; print d}")
  max_diff_antisense=$(awk "BEGIN {print ($diff > $max_diff_antisense) ? $diff : $max_diff_antisense}")
done

inconsistent=0
if (( $(awk "BEGIN {print ($max_diff_sense > 0.15)}") )) || (( $(awk "BEGIN {print ($max_diff_antisense > 0.15)}") )); then
  inconsistent=1
  echo "WARNING: Inconsistency detected across samples!" >> $out
  echo "  Max deviation from average: sense=${max_diff_sense}, antisense=${max_diff_antisense}" >> $out
  echo "  This may indicate different library prep batches or mixed strandness." >> $out
  echo "  Please review individual results above carefully." >> $out
  echo "" >> $out
fi

# Recommend strand parameter
# Rules for paired-end data:
# - Unstranded (XS): both sense and antisense ~0.5 (difference < 0.20)
# - Reverse stranded (RF): antisense > 0.70 (most common for dUTP/Illumina TruSeq)
# - Forward stranded (FR): sense > 0.70 (less common)

recommended_strand=""
abs_diff=$(awk "BEGIN {d = $avg_sense - $avg_antisense; if (d < 0) d = -d; print d}")

if (( $(awk "BEGIN {print ($abs_diff < 0.20)}") )); then
  recommended_strand="XS"
  explanation="Both sense and antisense fractions are similar (~50% each), indicating unstranded library."
elif (( $(awk "BEGIN {print ($avg_antisense > 0.70)}") )); then
  recommended_strand="RF"
  explanation="Antisense/Reverse fraction is dominant (>70%), indicating reverse stranded library (dUTP/Illumina TruSeq)."
elif (( $(awk "BEGIN {print ($avg_sense > 0.70)}") )); then
  recommended_strand="FR"
  explanation="Sense/Forward fraction is dominant (>70%), indicating forward stranded library (Ligation/SOLiD)."
else
  recommended_strand="UNCLEAR"
  explanation="Strandness pattern is ambiguous (neither fraction >70% and difference >20%). Manual review required."
fi

echo "RECOMMENDED --strand parameter: ${recommended_strand}" >> $out
echo "" >> $out
echo "Explanation: ${explanation}" >> $out
echo "" >> $out

if [ "$recommended_strand" = "UNCLEAR" ] || [ $inconsistent -eq 1 ]; then
  echo "ACTION REQUIRED:" >> $out
  echo "  Please manually review the results above and consult with your sequencing core" >> $out
  echo "  or bioinformatics team to confirm the library strandness." >> $out
else
  echo "Next step:" >> $out
  echo "  Use --strand ${recommended_strand} when running the junction analysis pipeline." >> $out
fi

echo "" >> $out
echo "Full raw output saved to: prepare/strandness_raw.txt" >> $out
echo "======================================================" >> $out

# Print summary to console
echo ""
echo "======================================================"
echo "STRANDNESS DETECTION COMPLETE"
echo "======================================================"
echo "Recommended --strand parameter: ${recommended_strand}"
if [ $inconsistent -eq 1 ]; then
  echo "WARNING: Inconsistency detected across samples!"
fi
echo ""
echo "Full report: prepare/strandness.txt"
echo "Raw output: prepare/strandness_raw.txt"
echo "======================================================"


