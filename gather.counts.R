library(ExomeDepth)
library(clustermq)
library(EDM)

foo = input.yaml
Q(gather.mosdepth.coverage,pkgs="EDM",task_id=1,n_jobs=1, max_calls_worker = 1,template=list(job_name="gather_coverage"), export = list(input = foo))
