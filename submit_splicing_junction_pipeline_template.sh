#!/bin/bash
#PBS -P 11003054
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -m bea
#PBS -N 1_splicing_junction

# These are needed modules in UT HPC to get singularity and Nextflow running.
# Replace with appropriate ones for your HPC.
#module load java-1.8.0_40
#module load singularity/3.5.3
#module load squashfs/4.4

# If you follow the eQTLGen phase II cookbook and analysis folder structure,
# some of the following paths are pre-filled.
# https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook

# We set the following variables for nextflow to prevent writing to your home directory.
export SINGULARITY_CACHEDIR=../../singularitycache
export NXF_HOME=../../nextflowcache

# Disable pathname expansion. Nextflow handles pathname expansion by itself.
set -f

cd $PBS_O_WORKDIR
mkdir -p logs

# Load an environment with Singularity and Java>=11
source ~/.bashrc
conda activate nf

# =========================================================
# Common settings
# =========================================================
nextflow_path=../../tools                                        # folder where Nextflow executable is
container_path=../singularity_img/eqtlgen_splicing_junction_v2.sif  # junction analysis container (regtools + leafcutter + python)
STAR_CONTAINER=../singularity_img/eqtlgen_splicing_star_v1.sif 

cohort_name=[your cohort name]
strand=[strandness: XS (unstranded) / RF (reverse stranded) / FR (forward stranded)]
output_path=../output

# =========================================================
# Option 1: Start from BAM files
# =========================================================
# Uncomment and fill in if your cohort already has aligned BAM files.
# Generates junction files, clustered results, and summary statistics.
# Cohorts can decide whether to share full results or just summary statistics.
#
# Note on strandness:
# - If you are CERTAIN about your library strandness, specify it directly in ${strand}
# - If you are UNCERTAIN, run submit_strandness.sh helper first to determine ${strand}
#   from your BAM files before running this pipeline

# BAM_DIR=[path to directory containing *.bam files]

# NXF_VER=21.10.6 ${nextflow_path}/nextflow run splicing_junction.sum_stat.nf \
#   --bam_dir        ${BAM_DIR} \
#   --outdir         ${output_path} \
#   --container      ${container_path} \
#   --cohort_name    ${cohort_name} \
#   --strand         ${strand} \
#   -profile pbs,singularity \
#   -resume \
#   -with-report  logs/report.html \
#   -with-trace   logs/trace.txt \
#   -with-timeline logs/timeline.html

# =========================================================
# Option 2: Start from FASTQ files (two-step workflow)
# =========================================================
# This workflow uses two linked pipelines:
#   Step 1: star_alignment.nf         (FASTQ → BAM)
#   Step 2: splicing_junction.sum_stat.nf (BAM → Junctions + Summary Stats)
#
# A checkpoint between steps allows you to:
#   - Review alignment quality
#   - Determine strandness if unknown (using submit_strandness.sh)

# FASTQ_DIR=[path to directory containing paired-end FASTQ files]
# READ_LENGTH=[read length, e.g. 150 for PE150]  # sjdbOverhang = READ_LENGTH - 1
# BAM_OUTPUT_DIR=../output/bam_alignment  # intermediate BAM output directory

# -----------------------------------------------------------
# Step 2.1: FASTQ → BAM (STAR alignment)
# -----------------------------------------------------------
# Align FASTQ files to BAM using STAR.
# After completion, review alignment quality before proceeding.

# NXF_VER=21.10.6 ${nextflow_path}/nextflow run star_alignment.nf \
#   --fastq_dir      ${FASTQ_DIR} \
#   --outdir         ${BAM_OUTPUT_DIR} \
#   --container      ${STAR_CONTAINER} \
#   --read_length    ${READ_LENGTH} \
#   -profile pbs,singularity \
#   -resume \
#   -with-report  logs/star_alignment_report.html \
#   -with-trace   logs/star_alignment_trace.txt \
#   -with-timeline logs/star_alignment_timeline.html

# -----------------------------------------------------------
# Step 2.2: Determine strandness (OPTIONAL - only if unknown)
# -----------------------------------------------------------
# If you are UNCERTAIN about library strandness, run submit_strandness.sh
# to infer the ${strand} value from generated BAM files.
# Point BAM_DIR to: ${BAM_OUTPUT_DIR}/aligned_bam/
#
# If you are CERTAIN about strandness, SKIP this step and proceed to Step 2.3.

# Edit submit_strandness.sh to set:
#   BAM_DIR=${BAM_OUTPUT_DIR}/aligned_bam
# Then submit:
#   qsub submit_strandness.sh
# Check result in prepare/strandness.txt for recommended --strand parameter.

# -----------------------------------------------------------
# Step 2.3: BAM → Junctions (junction analysis + summary stats)
# -----------------------------------------------------------
# After determining strandness (or if already known), run junction analysis.
# Use the ${strand} value from Step 2.2 or your known strandness.
# The BAM files from Step 2.1 will be in ${BAM_OUTPUT_DIR}/aligned_bam/

# NXF_VER=21.10.6 ${nextflow_path}/nextflow run splicing_junction.sum_stat.nf \
#   --bam_dir        ${BAM_OUTPUT_DIR}/aligned_bam \
#   --outdir         ${output_path} \
#   --container      ${container_path} \
#   --cohort_name    ${cohort_name} \
#   --strand         ${strand} \
#   -profile pbs,singularity \
#   -resume \
#   -with-report  logs/junction_report.html \
#   -with-trace   logs/junction_trace.txt \
#   -with-timeline logs/junction_timeline.html

# =========================================================
# Optional parameters for both BAM and FASTQ modes
# =========================================================
#   --regtools_memory_gb 1                     # default
#   --cluster_memory_gb  8                     # default
#   --sumstats_memory_gb 4                     # default

# =========================================================
# Optional FASTQ mode parameters only
# =========================================================
#   --fastq_pattern    '*_R{1,2}*.fastq.gz'   # default; adjust if filenames differ
#   --star_genome_dir  /references/hg38/assembly/STAR_GRCh38_gencode_v49  # default
#   --star_threads     16                      # default
#   --star_memory_gb   64                      # default
