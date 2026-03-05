#!/bin/bash
#PBS -P 11003054
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -m bea
#PBS -N 0_star_alignment

# These are needed modules in UT HPC to get singularity and Nextflow running.
# Replace with appropriate ones for your HPC.
#module load java-1.8.0_40
#module load singularity/3.5.3
#module load squashfs/4.4

# Load an environment with Singularity and Java>=11
source ~/.bashrc
conda activate nf

# If you follow the cookbook and analysis folder structure,
# some of the following paths are pre-filled.
# https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook

# We set the following variables for nextflow to prevent writing to your home directory.
export SINGULARITY_CACHEDIR=../../singularitycache
export NXF_HOME=../../nextflowcache

# Disable pathname expansion. Nextflow handles pathname expansion by itself.
set -f

cd $PBS_O_WORKDIR
mkdir -p logs

# =========================================================
# settings
# =========================================================

nextflow_path=../../tools
STAR_CONTAINER=../singularity_img/eqtlgen_splicing_star_v1.sif 
BAM_OUTPUT_DIR=../output/bam_alignment  # intermediate BAM output directory

FASTQ_DIR=[path to directory containing paired-end FASTQ files]
READ_LENGTH=[read length, e.g. 150 for PE150]  # sjdbOverhang = READ_LENGTH - 1


NXF_VER=21.10.6 ${nextflow_path}/nextflow run star_alignment.nf \
  --fastq_dir      ${FASTQ_DIR} \
  --outdir         ${BAM_OUTPUT_DIR} \
  --container      ${STAR_CONTAINER} \
  --read_length    ${READ_LENGTH} \
  -profile pbs,singularity \
  -resume \
  -with-report  logs/star_alignment_report.html \
  -with-trace   logs/star_alignment_trace.txt \
  -with-timeline logs/star_alignment_timeline.html

# =========================================================
# Optional FASTQ mode parameters 
# =========================================================
#   --fastq_pattern    '*_R{1,2}*.fastq.gz'   # default; adjust if filenames differ
#   --star_genome_dir  /references/hg38/assembly/STAR_GRCh38_gencode_v49  # default
#   --star_threads     16                      # default
#   --star_memory_gb   64                      # default


# =========================================================
# QC Report
# =========================================================
# After completion, review alignment quality before proceeding.
# The pipeline automatically generates: ../output/bam_alignment/star_alignment_qc_summary.tsv
# This TSV contains comprehensive QC metrics for all samples including:
#   - Input reads, uniquely mapped %, multi-mapped %, unmapped %
#   - Mismatch/indel rates, splice junction statistics
#   - QC pass/fail flags (PASS, WARNING_LOW_UNIQUE, FAIL_LOW_MAPPED)

# Most importantly, uniquely_mapped_pct & total_splices

