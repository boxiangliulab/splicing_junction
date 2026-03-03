nextflow.enable.dsl=2

def helpMessage() {
  log.info """
  =======================================================
   STAR Alignment: FASTQ → BAM v${workflow.manifest.version ?: 'dev'}
  =======================================================
  Usage:
    nextflow run star_alignment.nf \\
      --fastq_dir /path/to/FASTQ \\
      --outdir /path/to/output \\
      --container /path/to/eqtlgen_star_v1.sif \\
      --read_length 150 \\
      -profile slurm,singularity -resume

  Mandatory:
    --fastq_dir            Directory containing paired-end FASTQ files
    --outdir               Output directory for aligned BAM files
    --container            STAR container (must contain hg38 FASTA + GENCODE V49 GTF)
    --read_length          Read length (e.g. 150 for PE150); sets sjdbOverhang = read_length - 1

  Optional:
    --fastq_pattern        Glob for paired FASTQ files (default: '*_R{1,2}*.fastq.gz')
    --star_genome_dir      Path inside container to genome index (default: /references/hg38/assembly/STAR_GRCh38_gencode_v49)
    --star_threads         Threads for STAR genome build and alignment (default: 16)
    --star_memory_gb       Memory in GB for STAR genome build and alignment (default: 64)
    --use_prebuilt_index   Use pre-built genome index instead of generating one (default: false)

  Output:
    - Aligned BAM files in <outdir>/aligned_bam/
    - STAR logs in <outdir>/aligned_bam/
    - Gene counts in <outdir>/aligned_bam/

  Next step:
    After BAM generation, use the main pipeline in BAM mode:
      nextflow run splicing_junction_with_star.nf \\
        --bam_dir <outdir>/aligned_bam \\
        --outdir <junction_outdir> \\
        --container <junction_container> \\
        --cohort_name <name> \\
        --strand <XS|RF|FR>
  """.stripIndent()
}

def help                 = (params.get('help', false) as boolean)
def fastq_dir            = params.get('fastq_dir', '').toString().trim()
def outdir               = params.get('outdir', '').toString().trim()
def star_container_raw   = params.get('container', null)
def read_length          = (params.get('read_length', 0) as int)
def fastq_pattern        = params.get('fastq_pattern', '*_R{1,2}*.fastq.gz').toString().trim()
def star_genome_dir      = params.get('star_genome_dir', '/references/hg38/assembly/STAR_GRCh38_gencode_v49').toString().trim()
def star_threads         = (params.get('star_threads', 16) as int)
def star_memory_gb       = (params.get('star_memory_gb', 64) as int)
def use_prebuilt_index   = (params.get('use_prebuilt_index', false) as boolean)

if (help) { helpMessage(); System.exit(0) }

// Validate parameters
if (!fastq_dir)           exit 1, "[Pipeline error] --fastq_dir is required"
if (!outdir)              exit 1, "[Pipeline error] --outdir is required"
if (read_length < 1)      exit 1, "[Pipeline error] --read_length must be a positive integer (e.g. 150)"
if (star_container_raw == null) exit 1, "[Pipeline error] --container is required"

if (!new File(fastq_dir).exists())
  exit 1, "[Pipeline error] fastq_dir not found: ${fastq_dir}"

def star_cont = star_container_raw.toString().trim()
if (!star_cont || !new File(star_cont).exists())
  exit 1, "[Pipeline error] container not found: ${star_cont}"

def sjdbOverhang = read_length - 1

log.info "========================================="
log.info "Pipeline       : STAR Alignment (FASTQ → BAM)"
log.info "Version        : ${workflow.manifest.version ?: 'dev'}"
log.info "fastq_dir      : ${fastq_dir}"
log.info "fastq_pattern  : ${fastq_pattern}"
log.info "outdir         : ${outdir}"
log.info "container      : ${star_cont}"
log.info "read_length    : ${read_length}"
log.info "sjdbOverhang   : ${sjdbOverhang}"
log.info "star_genome_dir: ${star_genome_dir}"
log.info "star_threads   : ${star_threads}"
log.info "star_memory_GB : ${star_memory_gb}"
log.info "use_prebuilt   : ${use_prebuilt_index}"
log.info "workDir        : ${workflow.workDir}"
log.info "launchDir      : ${workflow.launchDir}"
log.info "profile        : ${workflow.profile}"
log.info "========================================="

// =====================
// STAR Processes
// =====================

process STAR_BUILD_INDEX {
  tag "sjdbOverhang_${sjdbOverhang}"
  publishDir "${outdir}/00_star_genome_index", mode: 'copy', overwrite: true

  container star_cont
  cpus star_threads
  memory "${star_memory_gb} GB"

  output:
    path "genome_index", emit: index

  script:
  """
  set -euo pipefail

  mkdir genome_index

  STAR \\
    --runMode genomeGenerate \\
    --runThreadN ${task.cpus} \\
    --genomeDir genome_index \\
    --genomeFastaFiles /references/hg38/assembly/GRCh38.primary_assembly.genome.fa \\
    --sjdbGTFfile /references/hg38/annotations/gencode.v49.annotation.gtf \\
    --sjdbOverhang ${sjdbOverhang}
  """
}

process STAR_ALIGN {
  tag { sample }
  publishDir "${outdir}/aligned_bam", mode: 'copy', overwrite: true

  container star_cont
  cpus star_threads
  memory "${star_memory_gb} GB"

  input:
    tuple val(sample), path(reads)
    val genome_index

  output:
    path "${sample}.bam", emit: bam
    path "${sample}.Log.final.out", emit: log
    path "${sample}.ReadsPerGene.out.tab", emit: counts

  script:
  def genome_param = use_prebuilt_index ? star_genome_dir : genome_index
  """
  set -euo pipefail

  STAR \\
    --runMode alignReads \\
    --runThreadN ${task.cpus} \\
    --twopassMode Basic \\
    --genomeDir ${genome_param} \\
    --readFilesIn ${reads[0]} ${reads[1]} \\
    --readFilesCommand zcat \\
    --outFileNamePrefix ${sample}. \\
    --outSAMtype BAM SortedByCoordinate \\
    --outSAMunmapped Within \\
    --outSAMattributes NH HI AS nM NM ch \\
    --outSAMstrandField intronMotif \\
    --outSAMattrRGline ID:${sample} SM:${sample} \\
    --sjdbOverhang ${sjdbOverhang} \\
    --outFilterMultimapNmax 20 \\
    --alignSJoverhangMin 8 \\
    --alignSJDBoverhangMin 1 \\
    --outFilterMismatchNmax 999 \\
    --outFilterMismatchNoverLmax 0.1 \\
    --alignIntronMin 20 \\
    --alignIntronMax 1000000 \\
    --alignMatesGapMax 1000000 \\
    --outFilterType BySJout \\
    --outFilterScoreMinOverLread 0.33 \\
    --outFilterMatchNminOverLread 0.33 \\
    --outFilterMatchNmin 0 \\
    --limitSjdbInsertNsj 1200000 \\
    --outFilterIntronMotifs None \\
    --alignSoftClipAtReferenceEnds Yes \\
    --quantMode GeneCounts \\
    --chimSegmentMin 15 \\
    --chimJunctionOverhangMin 15 \\
    --chimOutType Junctions WithinBAM SoftClip \\
    --chimMainSegmentMultNmax 1 \\
    --genomeLoad NoSharedMemory

  mv ${sample}.Aligned.sortedByCoord.out.bam ${sample}.bam
  """
}

// =====================
// Workflow
// =====================

workflow {
  def fastq_ch = Channel
    .fromFilePairs("${fastq_dir}/${fastq_pattern}", checkIfExists: true)
    .ifEmpty { exit 1, "No paired FASTQ files found in: ${fastq_dir} (pattern: ${fastq_pattern})" }

  if (use_prebuilt_index) {
    // Use pre-built genome index
    def genome_index_val = star_genome_dir
    STAR_ALIGN(fastq_ch, genome_index_val)
  } else {
    // Build cohort-specific genome index (sjdbOverhang = read_length - 1)
    def genome_index = STAR_BUILD_INDEX().index
    STAR_ALIGN(fastq_ch, genome_index)
  }
}

workflow.onComplete {
  println(workflow.success ? "Pipeline finished! BAM files in ${outdir}/aligned_bam/" : "Pipeline failed - check logs/workDir.")
  if (workflow.success) {
    println("Next step: Run junction analysis using splicing_junction.sum_stat.nf with --bam_dir ${outdir}/aligned_bam")
  }
}
