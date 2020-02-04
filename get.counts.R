library(ExomeDepth)
library(clustermq)
library(EDM)
nSamples = NROW(manifest);

foo = input.yaml

 Q(get.bam.counts.mosdepth,pkgs="EDM",task_id=1:nSamples,n_jobs=nSamples, max_calls_worker = 1, template=list(job_name="compute_coverage"), export = list(input = foo))
