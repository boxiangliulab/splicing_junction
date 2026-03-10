#!/bin/bash
#PBS -P 11003054
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -m bea
#PBS -N 3_splicing_junction

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
container_path=../singularity_img/eqtlgen_splicing_junction_v2.sif
output_path=../output

# Note on strandness:
# - If you are CERTAIN about your library strandness, specify it directly in ${strand}
# - If you are UNCERTAIN, run submit_strandness.sh helper first to determine ${strand}
#   from your BAM files before running this pipeline

BAM_DIR=[path to directory containing *.bam files]
cohort_name=[your cohort name]
strand=[strandness: XS (unstranded) / RF (reverse stranded) / FR (forward stranded)]


NXF_VER=21.10.6 ${nextflow_path}/nextflow run splicing_junction.sum_stat.nf \
  --bam_dir        ${BAM_DIR} \
  --outdir         ${output_path} \
  --container      ${container_path} \
  --cohort_name    ${cohort_name} \
  --strand         ${strand} \
  -profile pbs,singularity \
  -resume \
  -with-report  logs/report.html \
  -with-trace   logs/trace.txt \
  -with-timeline logs/timeline.html


# =========================================================
# Optional parameters for both BAM and FASTQ modes
# =========================================================
#   --regtools_memory_gb 1                     # default
#   --cluster_memory_gb  8                     # default
#   --sumstats_memory_gb 4                     # default

