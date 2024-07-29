nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process FILTER_SAMPLE_ADATA{
    publishDir "${out_dir}", mode: "$mode"

    input:
    tuple val(sample_id), path(input_path)


    output:
    path("*.h5ad"), emit: adata_filtered

    
	script:
	"""    
    002_filter_sample_adata.py \\
    --adata=${input_path}  \\
    --sample_id=${sample_id}
    """

}