## Background

[ExomeDepth](https://cran.r-project.org/web/packages/ExomeDepth/index.html) is one of the most sensitive tools for detecting copy-number variants (CNVs) from exome sequencing data. ExomeDepth is an R package and it is a challenge to deal with large datasets in R/ ExomeDepth. To address this issue, EDM implements the ExomeDepth pipeline with two major changes built in. First, EDM uses mosdepth to compute coverage over exons, and implements an approximation to the algorithm to select the reference panel for performance reasons. EDM uses [clustermq](https://cran.r-project.org/web/packages/clustermq/index.html) to run the pipeline in parallel (counting reads and variant calling). 

## Step 1: Setup the environment
EDM depends on several software tools including R packages. We use a conda environment to manage the dependencies. You can create a conda environment with all the dependencies using the edm_environment.yml file provided here. If you would like a different name for your environment, edit the yml file.

`git clone https://github.com/drramki-chop/EDM.git` 

`conda env create -f edm_environment.yml

conda activate edm_environment`

## Step 2: Install EDM package

Within R, install the EDM package from source.

`> install.packages("EDM_0.0.1.tar.gz",type="source",repos=NULL)`

## Step 3: Setup a template for [clustermq](https://cran.r-project.org/web/packages/clustermq/vignettes/userguide.html)

clustermq is installed as part of the conda environment. Now, we have to setup a template for resources for the cluster jobs.

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
CMQ_AUTH={{ auth }} R --no-save --no-restore -e 'clustermq:::worker("{{ master }}")'
```

Please refer to clustermq documentation for other HPC environments.

## Step 4: Remove exons with low mean mappability (optional; recommended)

We recently published a workflow ([Rajagopalan R et. al., 2020](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-0712-0) demonstrating that excluding exons with low mean mappability reduces the number of false-positives originating from the repetitive regions of the exome.

<span style="color:red">Caution:</span> This excludes some ~4.5% of the exons <span style="text-decoration:underline">incuding 0.6% of the exons that may be clinically-relevant</span>.

We provide the workflow to filter the exons with low mean mappability if you have your own bed file or use the exon definitions in ExomeDepth. However, you can simply use the `exons.hg19.mappability.filtered` object provided in the EDM package (`data(exons.hg19.mappability.filtered`).

### Workflow to filter the exons with low mean mappability (< 0.7)

In R:

```
> data("exons.hg19", package="ExomeDepth")
> data("exons.hg19.X", package = "ExomeDepth")
> exons.hg19 <- rbind(exons.hg19, exons.hg19.X)
> write.table(exons.hg19,"exons.hg19.bed",row.names =F,sep="\t",quote=F,col.names=F)

```
In shell:

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig
bigWigAverageOverBed wgEncodeCrgMapabilityAlign36mer.bigWig exons.hg19.bed exons.hg19.mappability.tab
cat exons.hg19.mappability.tab | awk '$NF >= 0.7' | cut -f1-4 > exons.hg19.mappability.bed  #for use in EDM workflow
```

## Step 5: Creat a manifest file for the samples

EDM forces a certain format for the workflow with minimal mandated metadata (bam, sampleID, sex). Column names should be the same in the manifest file (as the functions use them).

| bam   |      Are      |  Cool |
|:--------:|:-------------:|:-----:|
| bams/sample1.bam | ALGS-1P| F |
| bams/sample2.bam | ALGS-1M| M |
| bams/sample3.bam | ALGS-1F| F |










