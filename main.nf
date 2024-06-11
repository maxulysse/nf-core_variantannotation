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
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANTANNOTATION               } from './workflows/variantannotation'

include { PIPELINE_INITIALISATION         } from './subworkflows/local/utils_nfcore_variantannotation_pipeline'
include { ANNOTATION_CACHE_INITIALISATION } from './subworkflows/local/annotation_cache_initialisation'
include { DOWNLOAD_CACHE_SNPEFF_VEP       } from './subworkflows/local/download_cache_snpeff_vep'
include { GATHER_REPORTS_VERSIONS         } from './subworkflows/local/gather_reports_versions'
include { PIPELINE_COMPLETION             } from './subworkflows/local/utils_nfcore_variantannotation_pipeline'

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

fasta = params.fasta ? Channel.fromPath(params.fasta).map{ fasta -> [ [id:fasta.baseName], fasta ] }.collect() : Channel.empty()
vep_fasta = params.vep_include_fasta ? fasta : [ [id: 'null'], [] ]

bcftools_annotations  = params.bcftools_annotations  ? Channel.fromPath(params.bcftools_annotations).collect()  : Channel.empty()
bcftools_header_lines = params.bcftools_header_lines ? Channel.fromPath(params.bcftools_header_lines).collect() : Channel.empty()

vep_extra_files = []

if (params.dbnsfp && params.dbnsfp_tbi) {
    vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
    vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
}

if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi) {
    vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
    vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
    vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
    vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_VARIANTANNOTATION {

    take:
    samplesheet // channel: samplesheet read in from --input
    tools
    bcftools_annotations
    bcftools_annotations_tbi
    bcftools_header_lines
    snpeff_genome
    snpeff_cache
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
        tools,
        bcftools_annotations,
        bcftools_annotations_tbi,
        bcftools_header_lines,
        snpeff_genome,
        snpeff_cache,
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

    bcftools_annotations_tbi  = params.bcftools_annotations ? params.bcftools_annotations_tbi ? Channel.fromPath(params.bcftools_annotations_tbi).collect() : PREPARE_GENOME.out.bcftools_annotations_tbi : Channel.empty([])

    // Download cache
    if (params.download_cache) {
        // Assuming that even if the cache is provided, if the user specify download_cache, variantannotation will download the cache
        ensemblvep_info = Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
        snpeff_info     = Channel.of([ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], params.snpeff_genome, params.snpeff_db ])
        DOWNLOAD_CACHE_SNPEFF_VEP(ensemblvep_info, snpeff_info)
        snpeff_cache = DOWNLOAD_CACHE_SNPEFF_VEP.out.snpeff_cache
        vep_cache    = DOWNLOAD_CACHE_SNPEFF_VEP.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }

        versions = versions.mix(DOWNLOAD_CACHE_SNPEFF_VEP.out.versions)
    } else {
        // Looks for cache information either locally or on the cloud
        ANNOTATION_CACHE_INITIALISATION(
            (params.snpeff_cache && params.tools &&(params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
            params.snpeff_cache,
            params.snpeff_genome,
            params.snpeff_db,
            (params.vep_cache && params.tools &&(params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
            params.vep_cache,
            params.vep_species,
            params.vep_cache_version,
            params.vep_genome,
            "Please refer to https://nf-co.re/variantannotation/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.")

            snpeff_cache = ANNOTATION_CACHE_INITIALISATION.out.snpeff_cache
            vep_cache    = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache
    }
    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_VARIANTANNOTATION(
        PIPELINE_INITIALISATION.out.samplesheet.map{ meta, vcf -> [ meta + [ id:meta.sample ], vcf ] },
        params.tools,
        bcftools_annotations,
        bcftools_annotations_tbi,
        bcftools_header_lines,
        params.snpeff_genome ? "${params.snpeff_genome}.${params.snpeff_db}" : "${params.genome}.${params.snpeff_db}",
        snpeff_cache,
        vep_cache,
        params.vep_cache_version,
        vep_extra_files,
        vep_fasta,
        params.vep_genome,
        params.vep_species
    )

    // GATHER REPORTS/VERSIONS AND RUN MULTIC

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
