nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process QC{
    publishDir "${out_dir}", mode: "$mode"


    input:
    path(raw_adata)
    

    output:
     path("*adata.h5ad"), emit: adata_denoised


	script:
	"""    
    003_quality_control.py --adata_denoised=${adata_denoised}
    """
}