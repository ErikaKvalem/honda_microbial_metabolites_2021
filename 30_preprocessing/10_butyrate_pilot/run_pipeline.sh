nextflow run nf-core/rnaseq \
  -r 3.1 \
  -profile icbi,singularity \
  -c rnaseq.config \
  -w /data/scratch/sturm/projects/2021/microbial-metabolites/work
