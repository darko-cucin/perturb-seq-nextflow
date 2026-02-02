#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/perturbseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/perturbseq
    Website: https://nf-co.re/perturbseq
    Slack  : https://nfcore.slack.com/channels/perturbseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PERTURBSEQ  } from './workflows/perturbseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_perturbseq_pipeline'
//include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_perturbseq_pipeline'
//include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_perturbseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES CHECKING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//

workflow PERTURBSEQ_MAIN{

    take:
    samplesheet

    main:

    //
    // WORKFLOW: Run pipeline
    //

    PERTURBSEQ (
        samplesheet
    )
  
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    PIPELINE_INITIALISATION (
        params.version,           // boolean: Display version and exit
        params.validate_params,  // boolean: Boolean whether to validate parameters against the schema at runtime
        params.monochrome_logs,   // boolean: Do not use coloured log outputs
        args, //   array: List of positional nextflow CLI args
        params.outdir,
    )

    if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

    if (!params.index_reference_crispr && !(params.fasta_reference_crispr && params.gtf_reference_crispr)) {
      exit 1, "Either --index_reference_crispr or both --fasta_reference_crispr and --gtf_reference_crispr must be provided"
    }

    if (!params.index_reference_gex && !(params.fasta_reference_gex && params.gtf_reference_gex)) {
      exit 1, "Either --index_reference_gex or both --fasta_reference_gex and --gtf_reference_gex must be provided"
    }

   

    PERTURBSEQ_MAIN(
        ch_input
    )


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/




