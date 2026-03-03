nextflow.enable.dsl=2

def helpMessage() {
  log.info """
  =======================================================
   Junction Analysis with Summary Statistics v${workflow.manifest.version ?: 'dev'}
  =======================================================
  Usage:
    nextflow run splicing_junction.sum_stat.nf \\
      --bam_dir /path/to/BAM \\
      --outdir /path/to/output \\
      --container /path/to/sqtlgen_junction_v1.sif \\
      --cohort_name GTEx_LCL \\
      --strand XS \\
      -profile slurm,singularity -resume

  Mandatory:
    --bam_dir              Directory containing BAM files (*.bam)
    --outdir               Output directory
    --container            Junction analysis container (regtools + leafcutter)
    --cohort_name          Cohort name used for output file prefix
    --strand               XS (unstranded), RF (reverse stranded), FR (forward stranded)

  Optional:
    --regtools_memory_gb   Memory for regtools per sample in GB (default: 1)
    --cluster_memory_gb    Memory for leafcutter clustering in GB (default: 8)
    --sumstats_memory_gb   Memory for summary statistics in GB (default: 4)

  Output:
    - QC reports in <outdir>/00a_qc_mapq/
    - Junction files in <outdir>/01_junc/
    - Clustered junctions in <outdir>/02_cluster/
    - Summary statistics in <outdir>/03_summary_stats/
  """.stripIndent()
}

def valid_strands = ['XS', 'RF', 'FR']

def help                 = (params.get('help', false) as boolean)
def bam_dir              = params.get('bam_dir', '').toString().trim()
def outdir               = params.get('outdir', '').toString().trim()
def cohort_name          = params.get('cohort_name', '').toString().trim()
def strand               = params.get('strand', '').toString().trim()
def container_path_raw   = params.get('container', null)
def regtools_memory_gb   = (params.get('regtools_memory_gb', 1) as int)
def cluster_memory_gb    = (params.get('cluster_memory_gb', 8) as int)
def sumstats_memory_gb   = (params.get('sumstats_memory_gb', 4) as int)

if (help) { helpMessage(); System.exit(0) }

// Validate parameters
if (!bam_dir)     exit 1, "[Pipeline error] --bam_dir is required"
if (!outdir)      exit 1, "[Pipeline error] --outdir is required"
if (!cohort_name) exit 1, "[Pipeline error] --cohort_name is required"
if (!strand)      exit 1, "[Pipeline error] --strand is required (XS, RF, FR)"

if (!(strand in valid_strands))
  exit 1, "[Pipeline error] Invalid --strand: ${strand}. Allowed: XS, RF, FR"

if (!new File(bam_dir).exists())
  exit 1, "[Pipeline error] bam_dir not found: ${bam_dir}"

if (container_path_raw == null)
  exit 1, "[Pipeline error] --container is required"

def cont = container_path_raw.toString().trim()
if (!cont || !new File(cont).exists())
  exit 1, "[Pipeline error] container not found: ${cont}"

log.info "========================================="
log.info "Pipeline       : Junction Analysis + Summary Stats"
log.info "Version        : ${workflow.manifest.version ?: 'dev'}"
log.info "cohort_name    : ${cohort_name}"
log.info "bam_dir        : ${bam_dir}"
log.info "outdir         : ${outdir}"
log.info "container      : ${cont}"
log.info "strand         : ${strand}"
log.info "regtools_mem_GB: ${regtools_memory_gb}"
log.info "cluster_mem_GB : ${cluster_memory_gb}"
log.info "sumstats_mem_GB: ${sumstats_memory_gb}"
log.info "workDir        : ${workflow.workDir}"
log.info "launchDir      : ${workflow.launchDir}"
log.info "profile        : ${workflow.profile}"
log.info "========================================="

// =====================
// Junction Processes
// =====================

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

  # Sanity check: require NH tag (STAR/HISAT2 output)
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

  has_nh=\$(samtools view -F 4 ${bam} | head -2000 | grep -c \$'\\tNH:i:' || true)
  if [ "\$has_nh" -eq 0 ]; then
    echo "[ERROR] NH tag not found in BAM: ${bam}. This pipeline requires NH:i: tag (e.g., STAR/HISAT2 output)." >&2
    exit 1
  fi

  if [ "\$qc" = "ALL_UNIQ" ]; then
    ln -sfn ${bam} ${sample}.uniq.bam
    if [ -f ${bam}.bai ]; then
      ln -sfn ${bam}.bai ${sample}.uniq.bam.bai
    else
      samtools index ${sample}.uniq.bam
    fi
  else
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
    samtools index ${sample}.uniq.bam
  fi
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
    path "${cohort_name}_cluster_perind_numers.counts.gz", emit: counts

  script:
  """
  set -euo pipefail

  leafcutter_cluster_regtools \
    -j ${juncfile} -m 0 -p 0 -l 500000 \
    -o ${cohort_name}_cluster
  """
}

process COMPUTE_SUMMARY_STATS {
  publishDir "${outdir}/03_summary_stats", mode: 'copy', overwrite: true

  container cont
  cpus 1
  memory "${sumstats_memory_gb} GB"

  input:
    path counts_gz

  output:
    path "${cohort_name}_summary_stats.tsv"

  script:
  """
  #!/usr/bin/env python3
  import gzip
  import sys
  import numpy as np

  # Read the cluster counts file
  with gzip.open('${counts_gz}', 'rt') as f:
      # Read header (sample names)
      header = f.readline().strip().split()
      n_samples = len(header)

      # Prepare output
      with open('${cohort_name}_summary_stats.tsv', 'w') as out:
          # Write header
          out.write('junction_id\\tjunction_coord\\tmissingness\\ttotal_counts\\tmean_count\\tstd_count\\tcohort\\n')

          # Process each junction
          for line in f:
              fields = line.strip().split()
              if len(fields) < n_samples + 1:
                  continue

              junction_id = fields[0]
              counts = np.array([float(x) for x in fields[1:n_samples+1]])

              # Extract junction coordinates
              # Format: chr1:14829:14970:clu_1_- -> chr1:14829-14970
              parts = junction_id.split(':')
              if len(parts) >= 3:
                  junction_coord = f"{parts[0]}:{parts[1]}-{parts[2]}"
              else:
                  junction_coord = junction_id

              # Calculate statistics
              total_counts = int(np.sum(counts))
              mean_count = np.mean(counts)
              std_count = np.std(counts, ddof=1)  # Sample standard deviation

	      # New Statistics
              median_count = np.median(counts)
              min_count    = np.min(counts)
              max_count    = np.max(counts)

	      # Calculate missingness (fraction of samples with zero counts)
              n_zero = np.sum(counts == 0)
              missingness = n_zero / n_samples

              # Write output
	      out.write(f"{junction_id}\\t{junction_coord}\\t{missingness}\\t{total_counts}\\t{mean_count:.6f}\\t{std_count:.6f}\\t{median_count:.6f}\\t{min_count:.6f}\\t{max_count:.6f}\\t${cohort_name}\\n")

  print(f"Summary statistics computed successfully for ${cohort_name}", file=sys.stderr)
  """
}

// =====================
// Workflow
// =====================

workflow {
  def bam_ch = Channel
    .fromPath("${bam_dir}/*.bam", checkIfExists: true)
    .ifEmpty { exit 1, "No BAM files found in: ${bam_dir}" }
    .map { bam -> tuple(bam.baseName, bam, bam.getName()) }

  def qc_ch      = CHECK_UNIQ_MAPQ(bam_ch)
  def uniq_bam_ch = FILTER_UNIQ_BAM(qc_ch)
  def junc_ch    = EXTRACT_JUNCTIONS(uniq_bam_ch)

  def juncfile_ch = junc_ch
    .map { sample, junc -> junc.toString() }
    .collectFile(name: 'juncfile.txt', newLine: true)

  def cluster_counts = CLUSTER_JUNCTIONS(juncfile_ch)
  COMPUTE_SUMMARY_STATS(cluster_counts.counts)
}

workflow.onComplete {
  println(workflow.success ? "Pipeline finished!" : "Pipeline failed — check logs/workDir.")
}
