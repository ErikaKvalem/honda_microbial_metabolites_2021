#!/usr/bin/env nextflow

nextflow.enable.dsl=2



include { FILTER_SAMPLE_ADATA } from "./modules/filter_sample_adata"
include { SCAR_BY_SAMPLE } from "./modules/scar_by_sample"
include { QUALITY_CONTROL } from "./modules/quality_control"
include { MERGE_AND_SOLO_SAMPLES } from "./modules/merge_samples"


workflow {
    samplesheet_csv = Channel.fromPath(params.samplesheet_csv)
    //samplesheet_csv.splitCsv( skip:1).view()
   

    //input_path = Channel.fromPath(params.input_path)
    path_adata_denoised_after_qc = Channel.fromPath(params.path_adata_denoised_after_qc)
    
    
    // start workflow

    FILTER_SAMPLE_ADATA(samplesheet_csv.splitCsv( skip:1))

    ch_adata_filtered = FILTER_SAMPLE_ADATA.out.adata_filtered.flatten().map { 
       it -> [it.baseName.replace(".h5ad", ""), it]
    }

    SCAR_BY_SAMPLE(ch_adata_filtered)

    
   ch_adata_denoised = SCAR_BY_SAMPLE.out.adata_denoised.flatten().map { 
    it -> [it.baseName.replace(".h5ad", ""), it]
    }
    
   QUALITY_CONTROL(ch_adata_denoised)

   //PROPER WAY TO PASS ALL ADATAS TO MERGE
   //ch_adata_to_merge = QUALITY_CONTROL.out.adata_denoised_after_qc.collect() 
   //ch_adata_to_merge.view()

   MERGE_AND_SOLO_SAMPLES(QUALITY_CONTROL.out.start_proc,path_adata_denoised_after_qc)

 


}


