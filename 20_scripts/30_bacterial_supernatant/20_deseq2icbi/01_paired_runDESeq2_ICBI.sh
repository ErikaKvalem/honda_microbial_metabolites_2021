Rscript runDESeq2_ICBI.R /data/projects/2021/MicrobialMetabolites/bacterial-supernatant/10_rnaseq_pipeline/pipeline_info/samplesheet_group_organoid.valid.csv /data/projects/2021/MicrobialMetabolites/bacterial-supernatant/10_rnaseq_pipeline/star_salmon/salmon.merged.gene_counts.tsv  --c1=11mix --c2=10mix --nfcore  --condition_col=group --sample_col=sample --gtf_file=/data/genomes/mm39/annotation/gencode/gencode.vM27.primary_assembly.annotation.gtf --organism=mouse  --result_dir=/data/projects/2021/MicrobialMetabolites/bacterial-supernatant/20_deseq2icbi/paired_grp/deseq2_11mix_vs_10mix --remove_batch_effect --batch_col=organoid
