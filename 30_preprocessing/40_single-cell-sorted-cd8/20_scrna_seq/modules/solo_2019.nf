nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process SOLO_2019{
    publishDir "${out_dir}", mode: "$mode"


    input:
    path(denoised_filtered_adata)

    

    output:
     path("*.h5ad"), emit: adata_nodoublet


    
	script:
	"""    
    004_solo.py --denoised_filtered_adata=${denoised_filtered_adata}
    """

}