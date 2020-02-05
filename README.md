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
