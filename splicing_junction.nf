nextflow.enable.dsl=2

def helpMessage() {
  log.info """
  =======================================================
   JunctionStage1 v${workflow.manifest.version ?: 'dev'}
  =======================================================
  Usage:
    nextflow run JunctionStage1.nf \\
      --rnaseq_dir /path/to/BAM \\
      --outdir /path/to/output \\
      --container /path/to/sqtlgen_junction_v1.sif \\
      --cohort_name GTEx_LCL \\
      --strand XS \\
      -profile pbs,singularity -resume

  Mandatory:
    --rnaseq_dir
    --outdir
    --container
    --cohort_name
    --strand (XS, RF, FR)
  """.stripIndent()
}

def valid_strands = ['XS','RF','FR']

def help              = (params.get('help', false) as boolean)
def rnaseq_dir         = params.get('rnaseq_dir', '').toString().trim()
def outdir             = params.get('outdir', '').toString().trim()
def cohort_name        = params.get('cohort_name', '').toString().trim()
def strand             = params.get('strand', '').toString().trim()
def container_path_raw = params.get('container', null)

def regtools_memory_gb = (params.get('regtools_memory_gb', 1) as int)
def cluster_memory_gb  = (params.get('cluster_memory_gb', 8) as int)

if (help) { helpMessage(); System.exit(0) }

if (!rnaseq_dir)  exit 1, "[Pipeline error] --rnaseq_dir is required"
if (!outdir)      exit 1, "[Pipeline error] --outdir is required"
if (!cohort_name) exit 1, "[Pipeline error] --cohort_name is required"
if (!strand)      exit 1, "[Pipeline error] --strand is required (XS, RF, FR)"

if (!(strand in valid_strands))
  exit 1, "[Pipeline error] Invalid --strand: ${strand}. Allowed: XS, RF, FR"

if (!new File(rnaseq_dir).exists())
  exit 1, "[Pipeline error] rnaseq_dir not found: ${rnaseq_dir}"

if (container_path_raw == null)
  exit 1, "[Pipeline error] --container is required (or set params.container in nextflow.config)"

def cont = container_path_raw.toString().trim()
if (!cont)
  exit 1, "[Pipeline error] container is empty"

if (!new File(cont).exists())
  exit 1, "[Pipeline error] container not found: ${cont}"

log.info "========================================="
log.info "Pipeline Name     : JunctionStage1"
log.info "Pipeline Version  : ${workflow.manifest.version ?: 'dev'}"
log.info "cohort_name       : ${cohort_name}"
log.info "rnaseq_dir        : ${rnaseq_dir}"
log.info "outdir            : ${outdir}"
log.info "container         : ${cont}"
log.info "strand            : ${strand}"
log.info "regtools_mem_GB   : ${regtools_memory_gb}"
log.info "cluster_mem_GB    : ${cluster_memory_gb}"
log.info "workDir           : ${workflow.workDir}"
log.info "launchDir         : ${workflow.launchDir}"
log.info "profile           : ${workflow.profile}"
log.info "========================================="

process CHECK_UNIQ_MAPQ {
  tag { sample }

  publishDir "${outdir}/00a_qc_mapq", mode: 'copy', overwrite: true, pattern: "*.uniq.qc.txt"

  container cont
  cpus 1
  memory "4 GB"

  input:
    tuple val(sample), path(bam), val(bam_name)

  output:
    tuple val(sample), path(bam), path("${sample}.uniq.qc.txt")

  script:
  """
  set -euo pipefail

  # Sanity check: STAR should have NH tag
  has_nh=\$(samtools view -F 4 ${bam} | head -2000 | grep -c \$'\\tNH:i:' || true)
  if [ "\$has_nh" -eq 0 ]; then
    echo "[ERROR] NH tag not found in BAM: ${bam}. This pipeline requires NH:i: tag (e.g., STAR/HISAT2 output)." >&2
    exit 1
  fi

  total=\$(samtools view -c -F 4 ${bam})

  uniq=\$(samtools view -F 4 ${bam} \\
    | awk -F'\\t' '{
        nh="";
        for(i=12;i<=NF;i++){
          if(\$i ~ /^NH:i:/){ split(\$i,a,":"); nh=a[3]; break }
        }
        if(nh=="1") c++
      }
      END{ print c+0 }')

  if [ "\$total" -eq "\$uniq" ]; then
    echo "ALL_UNIQ" > ${sample}.uniq.qc.txt
  else
    echo "FILTER_REQUIRED" > ${sample}.uniq.qc.txt
  fi

  {
    echo "method\\tNH"
    echo "total_mapped\\t\$total"
    echo "uniq_mapped(NH=1)\\t\$uniq"
  } >> ${sample}.uniq.qc.txt
  """
}

process FILTER_UNIQ_BAM {
  tag { sample }

  container cont
  cpus 1
  memory "4 GB"

  input:
    tuple val(sample), path(bam), path(qcflag)

  output:
    tuple val(sample), path("${sample}.uniq.bam"), path("${sample}.uniq.bam.bai")

  script:
  """
  set -euo pipefail
  qc=\$(head -1 ${qcflag})

  # Require NH tag
  has_nh=\$(samtools view -F 4 ${bam} | head -2000 | grep -c \$'\\tNH:i:' || true)
  if [ "\$has_nh" -eq 0 ]; then
    echo "[ERROR] NH tag not found in BAM: ${bam}. This pipeline requires NH:i: tag (e.g., STAR/HISAT2 output)." >&2
    exit 1
  fi

  if [ "\$qc" = "ALL_UNIQ" ]; then
    ln -sfn ${bam} ${sample}.uniq.bam
  else
    # Keep mapped alignments with NH:i:1 (no MAPQ logic)
    samtools view -h -F 4 ${bam} \\
    | awk -F'\\t' 'BEGIN{OFS="\\t"}
        /^@/ {print; next}
        {
          nh="";
          for(i=12;i<=NF;i++){
            if(\$i ~ /^NH:i:/){ split(\$i,a,":"); nh=a[3]; break }
          }
          if(nh=="1") print
        }' \\
    | samtools view -b -o ${sample}.uniq.bam -
  fi

  samtools index ${sample}.uniq.bam
  """
}

process EXTRACT_JUNCTIONS {
  tag { sample }
  publishDir "${outdir}/01_junc", mode: 'copy', overwrite: true

  container cont
  cpus 1
  memory "${regtools_memory_gb} GB"

  input:
    tuple val(sample), path(ubam), path(bai)

  output:
    tuple val(sample), path("${sample}.junc")

  script:
  """
  set -euo pipefail

  regtools junctions extract \
    -a 8 -m 50 -M 500000 \
    -s ${strand} \
    -o ${sample}.junc \
    ${ubam}
  """
}

process CLUSTER_JUNCTIONS {
  publishDir "${outdir}/02_cluster", mode: 'copy', overwrite: true

  container cont
  cpus 1
  memory "${cluster_memory_gb} GB"

  input:
    path juncfile

  output:
    path "${cohort_name}_cluster_perind_numers.counts.gz"

  script:
  """
  set -euo pipefail

  leafcutter_cluster_regtools \
    -j ${juncfile} -m 0 -p 0 -l 500000 \
    -o ${cohort_name}_cluster
  """
}

workflow {

  def bam_ch = Channel
    .fromPath("${rnaseq_dir}/*.bam", checkIfExists: true)
    .filter { it != null }
    .ifEmpty { exit 1, "Input BAM files not found in: ${rnaseq_dir}" }
    .map { bam -> tuple(bam.baseName, bam, bam.getName()) }

  def qc_ch = CHECK_UNIQ_MAPQ(bam_ch)
  def uniq_bam_ch = FILTER_UNIQ_BAM(qc_ch)
  def junc_ch = EXTRACT_JUNCTIONS(uniq_bam_ch)

  def juncfile_ch = junc_ch
    .map { sample, junc -> junc.toString() }
    .collectFile(name: 'juncfile.txt', newLine: true)

  CLUSTER_JUNCTIONS(juncfile_ch)
}

workflow.onComplete {
  println(workflow.success ? "Pipeline finished!" : "Pipeline failed - check logs/workDir.")
}
