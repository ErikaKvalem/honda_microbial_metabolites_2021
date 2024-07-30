nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode


process QUALITY_CONTROL{
    publishDir "${out_dir}", mode: "$mode"
  

    input:
    tuple val(id), file(adata_denoised)
    

    output:
    path("*.h5ad")
    val true, emit: start_proc

	script:
	"""    
    004_quality_control.py --adata_denoised=${adata_denoised} --id=${id} 
    """
}




