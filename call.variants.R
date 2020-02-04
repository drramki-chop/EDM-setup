library(EDM)
library(ExomeDepth)
library(clustermq)

nSamples = NROW(manifest);

foo = input.yaml

Q(call.variants,pkgs=list("EDM","ExomeDepth"),columnIndex=1:nSamples,n_jobs=nSamples, max_calls_worker = 1,template=list(job_name="call_variants"), export = list(input =foo))


