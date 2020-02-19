## edm.submit.script.sh (SGE)

#$ -V
#$ -cwd
#$ -l mem_free=12g,h_vmem=12g

conda activate edm_env
Rscript -e 'EDM::wrapper.script(input.yaml)'
