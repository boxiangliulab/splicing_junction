# bulk sQTL workflow for eQTLGen: 1. splicing junction

We implement the sQTL workflow in two stages—splicing junction calling and sQTL calling—to ensure robust quality control of junctions. 

This is the cookbook for the first stage: splicing junction calling. In concordance with previous consortium efforts (e.g. GTEx, INTERVAL), we use the RegTools and LeafCutter pipeline. To harmonize junctions across all cohorts, filters for junction cluster (total reads) and junction (ratio in the group) are not applied in this workflow to generate per-cohort data. Instead, we will perform a stringent filtering based on the data from all cohorts and generate the harmonized splicing phenotypes for sQTL calling and all other downstream analysis. 

Below is the cookbook for per-cohort splicing junction file ({cohort_name}_cluster_perind_numers.counts.gz) preparation. 

For this cookbook, we follow the style of [eQTLGen Phase II cookbook](https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-cookbook.html) and aim to require minimal additional user actions. 



## Prerequisites

If you have run the workflow in eQTLGen Phase II cookbook, you should have all the following prerequisites. However, the current pipelines and configurations are tailored for use with the **PBS scheduler**. If your HPC uses a different scheduler, modify the job script **headers** to request the appropriate resources accordingly.


- High-Performance Computing (HPC) environment
  
- Singularity

- Nextflow 



## Getting Started 

### Setup

Make analysis root folder for this project, e.g. eQTLGen_splicing_phase1. Please ensure the folder is in data storage disk (not your home directory) to avoid issues. FYI, the output will take a similar srorage size of your input bam files. 

```
mkdir eQTLGen_splicing_phase1
```

Inside this folder make another folder for Nextflow executable, called tools. This path is specified in all the template scripts, so that scripts find Nextflow executable.

```
mkdir eQTLGen_splicing_phase1/tools
```

Inside eQTLGen_phase2/tools download and self-install Nextflow executable, as specified [here](https://www.nextflow.io/docs/latest/install.html#installation). You might need to load Java >=11 before running self-install (e.g. module load [Java >=11 module name in your HPC]).

```
cd eQTLGen_splicing_phase1/tools
wget -qO- https://get.nextflow.io | bash
```

You may refer to the [eQTLGen Phase II offline instructions](https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-offline-instructions.html) if your HPC has restrictions in connecting with internet. 

Inside eQTLGen_splicing_phase1, make additional separate folders for each step of this analysis plan:

```
cd ../
mkdir 1_splicing_junction
```

Inside each of those folders, you should clone/download the corresponding pipeline. If git is available in your HPC, easiest is to use command git clone [repo of the pipeline]. This yields a pipeline folder e.g. for the first pipeline it will look like that: eQTLGen_splicing_phase1/1_splicing_junction/splicing_junction.

```
cd 1_splicing_junction
git clone https://github.com/boxiangliulab/splicing_junction.git
```

We recommend to specify output (and, if needed, input etc.) folder(s) for each step. E.g. eQTLGen_phase2/1_splicing_junction/output.

```
cd ..
mkdir 1_splicing_junction/output
```

Inside each pipeline folder is the script template using name format submit_*_template.sh e.g. eQTLGen_splicing_phase1/1_splicing_junction/splicing_junction/submit_splicing_junction_pipeline_template.sh. You should adjust this according to your data (e.g. specify the path to your input folder, add required HPC modules), and save it. For better tracking, I usually rename it too: eQTLGen_splicing_phase1/1_splicing_junction/splicing_junction/submit_splicing_junction_pipeline_GTExWB_RNAseq.sh.

To ease the process, we have pre-filled some paths in the templates, assuming that you use the recommended folder structure.

You should submit each pipeline from inside each pipeline folder.

### Recources

Please download singularity images and gene annotation file from [Zenodo](https://zenodo.org/records/18067033) and check md5sum of the file. 

```
cd 1_splicing_junction
wget https://zenodo.org/records/18067033/files/eQTLGen_splicing_junction_v1.tar.gz
wget https://zenodo.org/records/18067033/files/eQTLGen_splicing_junction_v1.tar.gz.md5
# This should print eQTLGen_splicing_junction_v1.tar.gz: OK
md5sum -c eQTLGen_splicing_junction_v1.tar.gz.md5
```

Extract them in the current directory. 

```
tar -xzvf eQTLGen_splicing_junction_v1.tar.gz
```

You should be able to see two folders: data and singularity_img. Below is the full content by far: 

```
|-- 1_splicing_junction
|   |-- data
|   |   |-- gencode.v19.annotation.bed12
|   |   `-- gencode.v48.annotation.bed12
|   |-- eQTLGen_splicing_junction_v1.tar.gz
|   |-- eQTLGen_splicing_junction_v1.tar.gz.md5
|   |-- output
|   |-- singularity_img
|   |   |-- eqtlgen_splicing_junction_strandness.sif
|   |   `-- eqtlgen_splicing_junction_v1.sif
|   `-- splicing_junction
|       |-- README.md
|       |-- conf
|       |   |-- Danish_Computerome_profile.config
|       |   |-- JCTF_profile.config
|       |   |-- NSCC_profile.config
|       |   |-- base.config
|       |   |-- local_vm.config
|       |   |-- lsf_per_core.config
|       |   |-- lsf_per_job.config
|       |   |-- pbs.config
|       |   |-- sge.config
|       |   `-- slurm.config
|       |-- nextflow.config
|       |-- splicing_junction.nf
|       |-- submit_splicing_junction_pipeline_template.sh
|       `-- submit_strandness.sh
`-- tools
    `-- nextflow
```

Enter the main work directory, and prepare to run the pipelines.  

```
cd splicing_junction
```




## Data 

For this analysis plan you need the following information and data:

- Name of your cohort.
  
- Aligned BAM files

BAM files generated by STAR (recommended) or HISAT2 can be used. Importantly, if your library is unstranded, the BAM file must contain the XS tag on spliced reads. This is required for junction strand assignment in this workflow. To manually check if the XS tag is present in your BAM files, you can use ```samtools view [sample_name].bam | grep XS | head ```

If you are not sure whether your RNA-seq library is stranded or not, please run the short diagnostic script named submit_strandness.sh (please modify the headers according to your HPC scheduler). We use infer_experiment.py from [RSeQC](https://rseqc.sourceforge.net/) to check strandness of the first BAM and the last BAM files in your RNA-seq directory. 

```
# Please modify BAM_DIR before you run
# Adjust headers based on your scheduler
qsub submit_strandness.sh   # or equivalent command for your scheduler
```

You can check the results in prepare/strandness.log after the job completes. This usually takes less than 1 minute. For paired-end RNA-seq: 

       Usually, strand-specific RNA-seq libraries are prepared using a dUTP-based protocol, resulting in reverse-stranded reads (--strandness RF). In some cases (e.g. older datasets or non–strand-specific protocols), RNA-seq libraries may be unstranded (--strandness XS).
       If the fraction of reads explained by 1+-,1-+,2++,2-- is much larger than that explained by 1++,1--,2+-,2-+ and is close to 1             => RF (reverse-stranded, dUTP-based)
       If 1+-,1-+,2++,2-- is close to 1++,1--,2+-,2-+ and both are close to 0.5             => XS (unstranded)

       Below are two example outputs:

       Example 1: Paired-end, non–strand-specific: 
       This is PairEnd Data
       Fraction of reads failed to determine: 0.0172
       Fraction of reads explained by "1++,1--,2+-,2-+": 0.4903
       Fraction of reads explained by "1+-,1-+,2++,2--": 0.4925

       Example 2: Paired-end, strand-specific (reverse-stranded):
       This is PairEnd Data
       Fraction of reads failed to determine: 0.0072
       Fraction of reads explained by "1++,1--,2+-,2-+": 0.0487
       Fraction of reads explained by "1+-,1-+,2++,2--": 0.9441

If you encounter other scenarios or ambiguous results, please let us know. You can also refer to this [table](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/) for more details. 

### Genome Assembly

By default, we assume the genome assembly used for read alignment is hg38. Please inform us if you are using a different version (e.g. hg19). 



## Configure

modify the job file submit_splicing_junction_pipeline_[cohort name].sh in the main work directory. 

1. the absolute path of your input BAM directory. 

2. cohort name. 

3. strand parameter based on the results above. 

Now you are ready to go! 

```
qsub submit_splicing_junction_pipeline_[cohort name]   # or equivalent command for your scheduler
```



## Output

The pipeline will generate a output folder with junction files for analysis. 

```
output/
├── 00a_qc_mapq
|   |-- {sample_name}.uniq.qc.txt
|   `-- ... .uniq.qc.txt
├── 01_junc
|   |-- {sample_name}.junc
|   `-- ... .junc
`-- 02_cluster
    `-- {cohort_name}_cluster_perind_numers.counts.gz
```

Send output file to central analysis for joint junction quality control & filtering: 
Please upload {cohort_name}_cluster_perind_numers.counts.gz via your preferred file-sharing service (e.g., Dropbox/Globus) and share with us the link. 

Please also submit relevant metadata (if available) of your cohort. For example, ancestry, sex and age can be included. This can be a table in either csv/tsv or txt format. 



## Runtime, memory, and storage

Resource usage estimates below were obtained from test runs on the **GTEx LCL cohort (N = 327 samples)**, with a **total input BAM size of ~2 TB**. All benchmarks were performed on a shared HPC filesystem (NFS-based /scratch storage), which is optimized for high-throughput I/O and concurrent access across compute nodes.

Actual resource usage may vary depending on BAM file size, sequencing depth, read length, and filesystem performance.

| Task | Threads | Memory | Runtime | Disk Usage |
|------|--------:|-------:|--------:|-----------:|
| `check_unique_mapq` | 1 | ~0.1 GB | ~1–2 min | – |
| `filter_unique_bam` | 1 | ~0.1 GB | ~5–10 min | ~6 GB per sample (~2 TB total) |
| `extract_junctions` (RegTools) | 1 | ~0.2 GB | ~1–2 min | ~20 MB per sample (~6.5 GB total) |
| `cluster_junctions` (LeafCutter) | 1 | ~0.6 GB | ~24 min (cohort-level) | ~1 GB total |



## Troubleshooting

You can check the scheduler-generated stdout/stderr log files (e.g. 1_splicing_junction.o* in PBS) or inspect .nextflow.log.

Please also refer to Section "Running, monitoring and debugging Nextflow pipelines" in [eQTLGen Phase II cookbook](https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-cookbook.html). 




## Citation

Nextflow:
[Di Tommaso, P., Chatzou, M., Floden, E. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316–319 (2017). https://doi.org/10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820)

RSeQC: 
[Wang L, Wang S, Li W. RSeQC: quality control of RNA-seq experiments. Bioinformatics. 2012 Aug 15;28(16):2184-5. doi: 10.1093/bioinformatics/bts356. Epub 2012 Jun 27. PMID: 22743226.](https://academic.oup.com/bioinformatics/article/28/16/2184/325191?login=false)

Regtools: 
[Cotto KC, Feng YY, Ramu A, Richters M, Freshour SL, Skidmore ZL, Xia H, McMichael JF, Kunisaki J, Campbell KM, Chen TH, Rozycki EB, Adkins D, Devarakonda S, Sankararaman S, Lin Y, Chapman WC, Maher CA, Arora V, Dunn GP, Uppaluri R, Govindan R, Griffith OL, Griffith M. Integrated analysis of genomic and transcriptomic data for the discovery of splice-associated variants in cancer. Nature Communications. 2023 Mar. pmid: 36949070; doi: 10.1038/s41467-023-37266-6](https://www.nature.com/articles/s41467-023-37266-6)

LeafCutter:
[Li YI, Knowles DA, Humphrey J, Barbeira AN, Dickinson SP, Im HK, Pritchard JK. Annotation-free quantification of RNA splicing using LeafCutter. Nat Genet. 2018 Jan;50(1):151-158. doi: 10.1038/s41588-017-0004-9. Epub 2017 Dec 11. PMID: 29229983; PMCID: PMC5742080.](https://www.nature.com/articles/s41588-017-0004-9)



## Contact

Please contact e1101920@u.nus.edu if you have any other questions. 

