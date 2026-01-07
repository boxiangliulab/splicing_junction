#!/bin/bash
#PBS -P 11003054
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=02:00:00
#PBS -N sQTLGen_strandness
#PBS -j oe
#PBS -o prepare/strandness.log

cd "$PBS_O_WORKDIR"
mkdir -p prepare

# ----------------------------
# Paths
# ----------------------------
ROOT=$PWD/..
SIF=${ROOT}/singularity_img/eqtlgen_splicing_junction_strandness.sif 
BAM_DIR=[your BAM file directory]

BED=${ROOT}/data/gencode.v48.annotation.bed12

# if you are using hg19, please use the bed below  
#BED=${ROOT}/data/gencode.v19.annotation.bed12

echo "[INFO] Host: $(hostname)"
echo "[INFO] PWD:  $(pwd)"
echo "[INFO] ROOT: ${ROOT}"
echo "[INFO] SIF:  ${SIF}"
echo "[INFO] BAM_DIR: ${BAM_DIR}"
echo "[INFO] BED: ${BED}"
echo "[INFO] Date: $(date)"

# pick the first and the last BAMs as representatives
bam1=$(ls "${BAM_DIR}"/*.bam | head -1)
bamN=$(ls "${BAM_DIR}"/*.bam | tail -1)

bam1_dir=$(dirname "${bam1}")
bamN_dir=$(dirname "${bamN}")



# output files
out="prepare/strandness.txt"

echo "[INFO] bam1 = ${bam1}" > $out
echo "========== RUN 1 ==========" >> $out
singularity exec --userns --containall \
  -B "${ROOT}":/work \
  -B "${bam1_dir}:${bam1_dir}" \
  --pwd /work \
  "${SIF}" \
  infer_experiment.py \
  -s 1000000 \
  -r /work/data/$(basename "${BED}") \
  -i "${bam1}" \
  >> ${out} 2>&1

echo "[INFO] bamN = ${bamN}" >> $out
echo "========== RUN 2 ==========" >> $out
singularity exec --userns --containall \
  -B "${ROOT}":/work \
  -B "${bamN_dir}:${bamN_dir}" \
  --pwd /work \
  "${SIF}" \
  infer_experiment.py \
  -s 1000000 \
  -r /work/data/$(basename "${BED}") \
  -i "${bamN}" \
  >> ${out} 2>&1


