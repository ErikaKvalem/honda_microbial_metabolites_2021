nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process MERGE_AND_SOLO_SAMPLES{
    publishDir "${out_dir}", mode: "$mode"
    label "gpu"

    input:
    val ready
    val path_adata_denoised_after_qc
    

    output:
     path("*.h5ad"), emit: adata_merged

	script:
	"""    
    005_006_merge_solo.py --path_adata_denoised_after_qc=${path_adata_denoised_after_qc}
    """
}