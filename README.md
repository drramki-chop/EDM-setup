**** Under development *****


## Background

[ExomeDepth](https://cran.r-project.org/web/packages/ExomeDepth/index.html) is one of the most sensitive tools for detecting copy-number variants (CNVs) from exome sequencing data. ExomeDepth is an R package and it is a challenge to deal with large datasets in R/ ExomeDepth. To address this issue, EDM implements the ExomeDepth pipeline with two major changes built in. First, EDM uses mosdepth to compute coverage over exons, and implements an approximation to the algorithm to select the reference panel for performance reasons. EDM uses [clustermq](https://cran.r-project.org/web/packages/clustermq/index.html) to run the pipeline in parallel (counting reads and variant calling). 

## Step 1: Setup the environment
EDM depends on several software tools including R packages. We use a conda environment to manage the dependencies. You can create a conda environment with all the dependencies using the edm_environment.yml file provided here. If you would like a different name for your environment, edit the yml file.

```
git clone https://github.com/drramki-chop/EDM.git
conda env create -f edm_conda_env.yml
conda activate edm_environment
```

If you'd like to install conda from scratch, follow the protocol here:
https://conda.io/projects/conda/en/latest/user-guide/install/linux.html

## Step 2: Install EDM package

Within R, install the EDM package from source.

```r
devtools::install_github("drramki/edm")
```

## Step 3: Setup a template for [clustermq](https://cran.r-project.org/web/packages/clustermq/vignettes/userguide.html)

clustermq is installed as part of the conda environment. EDM comes with default templates for clustermq (and specific HPC environments) that can be specified in the configuration file (.yml in step 6). However, you can setup the template for resources for the cluster jobs independently (~/.clustermq_template).

For SGE:

```
#$ -N {{ job_name }}               # job name
#$ -q default                      # submit to queue named "default"
#$ -j y                            # combine stdout/error in one file
#$ -o {{ log_file | /dev/null }}   # output file
#$ -cwd                            # use pwd as work dir
#$ -V                              # use environment variable
#$ -t 1-{{ n_jobs }}               # submit jobs as array
#$ -pe {{ cores | 1 }}             # number of cores to use per job

ulimit -v $(( 1024 * {{ memory | 4096 }} ))
conda activate edm_environment
CMQ_AUTH={{ auth }} R --no-save --no-restore -e 'clustermq:::worker("{{ master }}")'
```

Please refer to [clustermq](https://cran.r-project.org/web/packages/clustermq/vignettes/userguide.html) documentation for other HPC environments.

## Step 4: Remove exons with low mean mappability (optional; recommended)

We recently published a workflow ([Rajagopalan R et. al., 2020](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-0712-0) demonstrating that excluding exons with low mean mappability reduces the number of false-positives originating from the repetitive regions of the exome while maintaining the same sensitivity.

```diff
- This excludes some ~4.5% of the exons incuding 0.6% of the exons that may be clinically-relevant.
```
We provide the workflow to filter the exons with low mean mappability if you have your own bed file or use the exon definitions in ExomeDepth. However, you can simply use the `exons.hg19.mappability.filtered` object provided in the EDM package (`data(exons.hg19.mappability.filtered`).

### Workflow to filter the exons with low mean mappability (< 0.7)

In R (with default exon definitions from ExomeDepth):

```r
> data("exons.hg19", package="ExomeDepth")
> data("exons.hg19.X", package = "ExomeDepth")
> exons.hg19 <- rbind(exons.hg19, exons.hg19.X)
> exons.hg19$name <- paste0(exons.hg19$chromosome,":",exons.hg19$start,"-",exons.hg19$end,"_",exons.hg19$name)
> write.table(exons.hg19,"exons.hg19.bed",row.names =F,sep="\t",quote=F,col.names=F)

```
In shell:

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig
bigWigAverageOverBed wgEncodeCrgMapabilityAlign36mer.bigWig exons.hg19.bed exons.hg19.mappability.tab
cat exons.hg19.mappability.tab | awk '$NF >= 0.7' | cut -f1-3 > exons.hg19.mappability.bed  #for use in EDM workflow
```

If you have your own exon definitions/ BED file, you can directly use the bigWigAverageOverBed on the file.

Note: Using mappability or alignability should produce similar results while the difference between the choice of k-mer length may be noticeable (36-mer vs 100-mer).

## Step 5: Creat a manifest file for the samples

EDM forces a certain format for the workflow with minimal mandated metadata (bam, sampleID, sex). 
```diff
- Column names should be the same in the manifest file (as the functions use them).
```

| bam   |      sampleID      |  sex |
|:--------:|:-------------:|:-----:|
| bams/sample1.bam | ALGS-1P| F |
| bams/sample2.bam | ALGS-1M| M |
| bams/sample3.bam | ALGS-1F| F |

 - F - female; M - male
 - Absolute paths are required for bam files
 
## Step 6: Create a configuration file (.yml)

EDM requires a configuration file with four mandated fields (cohort.name, manifest, outputDir, scheduler). Providing `transition.probability` is optional. We recommend a transition probability of 1e-8 for spcificity and 1e-4 for more sensitivity.

```yaml
# required
cohort.name: Epi4k

# required
manifest: pediseq.edm.manifest.txt

# required
output.directory: ./

# required
transition.probability: 1e-4

# required
bed.file: # bedfile path or leave blank

## required
hpc.scheduler: sge

##optional 
clustermq.template: ## user defined template - if provided certain options may change.
```

## Step 7: Run the pipeline

EDM is specifically designed to run in high-performance (HPC) computing environments and suited for large cohorts. One can run smaller cohorts as well but EDM may not make a difference in terms of computational performance. But it provides a seamless automated workflow.

```
## edm.submit.script.sh (SGE)

#$ -V
#$ -cwd
#$ -l mem_free=12g,h_vmem=12g

conda activate edm_env
Rscript -e 'EDM::wrapper.script('input.yaml')'

```

```
qsub ./edm.submit.script.sh
```


## Considerations

### General
1. If EDM fails, think about following reasons:
    - See if the input files conforms to the requirements (input.yml/ manifest/ bed file) and appropriate permissions exist for the output folders.
    - Absolte paths are provided for bam files and the index files should be in the same folder.
    - See if there is enough memory/ resources provided in the template.

2. ExomeDepth performs best when the cohort is homogeneous in terms of sample preparation, library prep, exome capture kit,  sequencing platform and sequencing center. If you gather samples from different sources, this may not be the best tool. Always test the pipeline first with known variants.

3. CNV calling in chromosome X is challenging as it depends only on ~6000 exons. A decent number of samples is required to get reliable results.

4. We recommend a transition probability of 1e-8 for specificity. If you prefer more sensitivity, 1e-4 is the default recommended by the ExomeDepth.

###  For larger cohorts

1. We have tested this workflow with ~2000 samples. We anticipate this pipeline will work for several thousands of samples but with more memory. The submit script keeps a master thread alive (that runs worker threads/ array jobs). The memory requirement for the master thread increases linearly with the number of samples. 

2. If your cohort has a large number of samples, perhaps, you can consider splitting the cohort by some factor. For example, by sequencing time, batch, or center. We are exploring ways to scale this pipeline for larger cohorts as well.

## To do

 - Gene-based annotations for downstream analysis
 - Case-control analysis
 - VCF export
 - Pedigree-based analysis for recessive variants

### Reference

1. [A Robust Model for Read Count Data in Exome Sequencing Experiments and Implications for Copy Number Variant Calling ](https://pubmed.ncbi.nlm.nih.gov/22942019), Bioinformatics, 28(21), 2747-54, 2012. 
2. [A highly sensitive and specific workflow for detecting rare copy-number variants from exome sequencing data](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-0712-0), Genome Medicine, 12(1), 14, 2020.

## Authors
 - Ramakrishnan Rajagopalan (rajagopalanr@email.chop.edu) [Google Scholar](https://scholar.google.com/citations?hl=en&user=OCUSDOkAAAAJ&view_op=list_works&sortby=pubdate)
 - Shiva Ganesan (ganesans@email.chop.edu)
