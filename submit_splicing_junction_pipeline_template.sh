#!/bin/bash
#PBS -P 11003054
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -m bea
#PBS -N 1_splicing_junction

# These are needed modules in UT HPC to get singularity and Nextflow running. Replace with appropriate ones for your HPC.
#module load java-1.8.0_40
#module load singularity/3.5.3
#module load squashfs/4.4

# If you follow the eQTLGen phase II cookbook and analysis folder structure,
# some of the following paths are pre-filled.
# https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook

# We set the following variables for nextflow to prevent writing to your home directory (and potentially filling it completely)
# Feel free to change these as you wish.
export SINGULARITY_CACHEDIR=../../singularitycache
export NXF_HOME=../../nextflowcache

# Disable pathname expansion. Nextflow handles pathname expansion by itself.
set -f


cd $PBS_O_WORKDIR
mkdir -p logs

# load an environment with Singularity and Java>=11
source ~/.bashrc
conda activate nf

# Disable pathname expansion. Nextflow handles pathname expansion by itself.
set -f

# Define paths (keep the same naming style as eQTLGen phase2)
nextflow_path=../../tools  # folder where Nextflow executable is
container_path=../singularity_img/eqtlgen_splicing_junction_v1.sif  # pre-filled

# ---- Inputs ----
BAM_DIR=[your BAM file directory]

# ---- Metadata ----
cohort_name=[your cohort name]
strand=[strandness parameter for RegTools]   # XS or RF or FR (mandatory)

# ---- Output ----
output_path=../output  # pre-filled 

# Command:
NXF_VER=21.10.6 ${nextflow_path}/nextflow run splicing_junction.nf \
  --rnaseq_dir ${BAM_DIR} \
  --outdir ${output_path} \
  --container ${container_path} \
  --cohort_name ${cohort_name} \
  --strand ${strand} \
  -profile pbs,singularity \
  -resume \
  -with-report logs/report.html \
  -with-trace logs/trace.txt \
  -with-timeline logs/timeline.html

