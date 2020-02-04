#$ -V
#$ -cwd
#$ -l mem_free=12g,h_vmem=12g

source activate edm_env

Rscript wrapper.R --input example.input.yaml
