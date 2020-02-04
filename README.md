## Background

[ExomeDepth]("https://cran.r-project.org/web/packages/ExomeDepth/index.html") is one of the most sensitive tools for detecting copy-number variants (CNVs) from exome sequencing data. EDM implements the ExomeDepth pipeline with two major changes built in. First, EDM uses mosdepth to compute coverage over exons, and implements an approximation to the algorithm to select the reference panel for performance reasons. EDM uses [clustermq]("https://cran.r-project.org/web/packages/clustermq/index.html") to run the pipeline in parallel (counting reads and variant calling). 

## Step 1: Setup the environment
EDM depends on several software tools including R packages. We use a conda environment to manage the dependencies. You can create a conda environment with all the dependencies using the edm_environment.yml file provided here. If you would like a different name for your environment, edit the yml file.

conda env create -f edm_environment.yml


## Step 2: Install EDM package
