#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/variantannotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/variantannotation
    Website: https://nf-co.re/variantannotation
    Slack  : https://nfcore.slack.com/channels/variantannotation
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta             = getGenomeAttribute('fasta')
params.snpeff_db         = getGenomeAttribute('snpeff_db')
params.snpeff_genome     = getGenomeAttribute('snpeff_genome')
params.vep_cache_version = getGenomeAttribute('vep_cache_version')
params.vep_genome        = getGenomeAttribute('vep_genome')
params.vep_species       = getGenomeAttribute('vep_species')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANTANNOTATION         } from './workflows/variantannotation'
include { PREPARE_VARIANTANNOTATION } from './subworkflows/local/prepare_variantannotation'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_variantannotation_pipeline'
include { GATHER_REPORTS_VERSIONS   } from './subworkflows/local/utils_nfcore_variantannotation_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_variantannotation_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow NFCORE_VARIANTANNOTATION {

    take:
    samplesheet // channel: samplesheet read in from --input
    bcfann_enabled
    bcftools_annotations
    bcftools_annotations_tbi
    bcftools_header_lines
    snpeff_enabled
    snpeff_genome
    snpeff_cache
    merge_enabled
    vep_enabled
    vep_cache
    vep_cache_version
    vep_extra_files
    vep_fasta
    vep_genome
    vep_species

    main:
    //
    // WORKFLOW: Run pipeline
    //
    VARIANTANNOTATION(
        samplesheet,
        bcfann_enabled,
        bcftools_annotations,
        bcftools_annotations_tbi,
        bcftools_header_lines,
        snpeff_enabled,
        snpeff_genome,
        snpeff_cache,
        merge_enabled,
        vep_enabled,
        vep_cache,
        vep_cache_version,
        vep_extra_files,
        vep_fasta,
        vep_genome,
        vep_species)

    emit:
    reports  = VARIANTANNOTATION.out.reports  // channel: [ path(reports) ]
    versions = VARIANTANNOTATION.out.versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    versions = Channel.empty()

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // SUBWORKFLOW: Handle local or cloud annotation cache
    //     Or alternatevly download cache
    //     And index file if needed
    //
    PREPARE_VARIANTANNOTATION(
        params.fasta,
        (params.tools && (params.tools.split(',').contains("bcfann")) && params.bcftools_annotations),
        params.bcftools_annotations,
        params.bcftools_annotations_tbi,
        params.bcftools_header_lines,
        params.download_cache,
        "Please refer to https://nf-co.re/variantannotation/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.",
        (params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
        params.snpeff_cache,
        params.snpeff_genome,
        params.snpeff_db,
        (params.tools && (params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
        params.vep_cache,
        params.vep_species,
        params.vep_cache_version,
        params.vep_include_fasta,
        params.vep_genome,
        params.vep_custom_args,
        params.dbnsfp,
        params.dbnsfp_tbi,
        params.spliceai_snv,
        params.spliceai_snv_tbi,
        params.spliceai_indel,
        params.spliceai_indel_tbi)

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_VARIANTANNOTATION(
        PIPELINE_INITIALISATION.out.samplesheet.map{ meta, vcf -> [ meta + [ id:meta.sample ], vcf ] },
        (params.tools && (params.tools.split(',').contains("bcfann"))),
        PREPARE_VARIANTANNOTATION.out.bcftools_annotations,
        PREPARE_VARIANTANNOTATION.out.bcftools_annotations_tbi,
        PREPARE_VARIANTANNOTATION.out.bcftools_header_lines,
        (params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
        params.snpeff_genome ? "${params.snpeff_genome}.${params.snpeff_db}" : "${params.genome}.${params.snpeff_db}",
        PREPARE_VARIANTANNOTATION.out.snpeff_cache,
        (params.tools && (params.tools.split(',').contains("merge"))),
        (params.tools && (params.tools.split(',').contains("vep"))),
        PREPARE_VARIANTANNOTATION.out.vep_cache,
        params.vep_cache_version,
        PREPARE_VARIANTANNOTATION.out.vep_extra_files,
        PREPARE_VARIANTANNOTATION.out.vep_fasta,
        params.vep_genome,
        params.vep_species
    )

    //
    // GATHER REPORTS/VERSIONS AND RUN MULTIC
    //
    GATHER_REPORTS_VERSIONS(
        params.outdir,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        NFCORE_VARIANTANNOTATION.out.versions,
        NFCORE_VARIANTANNOTATION.out.reports
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        GATHER_REPORTS_VERSIONS.out.report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
