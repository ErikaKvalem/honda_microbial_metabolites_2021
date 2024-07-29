nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process SCAR_2019{
    publishDir "${out_dir}", mode: "$mode"


    input:
    path(raw_adata)
    path(filtered_adata)
    

    output:
     path("*.h5ad"), emit: adata_denoised


    
	script:
	"""    
    002_scAR.py --raw_adata=${raw_adata} --filtered_adata=${filtered_adata}
    """

}