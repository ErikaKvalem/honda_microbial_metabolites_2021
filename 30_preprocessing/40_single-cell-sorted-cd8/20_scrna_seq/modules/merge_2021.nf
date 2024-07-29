nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process MERGE_2021{
    publishDir "${out_dir}", mode: "$mode"

    input:
    tuple val(experiment)
    

    output:
    path("*adata.h5ad"), emit: adata

    
	script:
	"""    
     001_merge.py \\
    --count_dir= \\
    --experiment=${experiment} 
    """

}