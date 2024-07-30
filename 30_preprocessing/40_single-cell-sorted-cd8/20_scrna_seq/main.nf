#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//include { SCAR_2019 } from "./modules/scar_2019"
//include { SCAR_2021 } from "./modules/scar_2021"
//include { SOLO_2019 } from "./modules/solo_2019"
//include { SOLO_2021 } from "./modules/solo_2021"
include { SOLO } from "./modules/solo"


workflow {
    //raw_adata_2019 = Channel.fromPath(params.raw_input_path_2019)
    //filtered_adata_2019 = Channel.fromPath(params.filtered_input_path_2019)
    //raw_adata_2021 = Channel.fromPath(params.raw_input_path_2021)
    //filtered_adata_2021 = Channel.fromPath(params.filtered_input_path_2021)

    //denoised_filtered_adata_2019 = Channel.fromPath(params.denoised_filtered_input_path_2019)
    //denoised_filtered_adata_2021 = Channel.fromPath(params.denoised_filtered_input_path_2021)
    denoised_colon_and_tumor = Channel.fromPath(params.denoised_filtered_input_path_colon_and_tumor)
    

    // start workflow
    //SCAR_2019(raw_adata_2019,filtered_adata_2019)
    //SCAR_2021(raw_adata_2021,filtered_adata_2021)
    //SOLO_2019(denoised_filtered_adata_2019)
    //SOLO_2021(denoised_filtered_adata_2021)
    SOLO(denoised_colon_and_tumor)
    
}


