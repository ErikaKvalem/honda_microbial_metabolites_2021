nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process SCAR_BY_SAMPLE{
    publishDir "${out_dir}", mode: "$mode"

    input:
    tuple val(id), file(ch_adata_filtered)

    output:
    path("*.h5ad"), emit: adata_denoised
    
	script:    
    """
    003_scar.py --filtered_adata=${ch_adata_filtered} --id=${id}
	"""

}