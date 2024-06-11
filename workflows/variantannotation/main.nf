/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VCF_ANNOTATE_ALL                            } from '../../subworkflows/local/vcf_annotate_all/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTANNOTATION {
    take:
    samplesheet              // channel (queue): samplesheet read in from --input
    tools                    // array
    bcftools_annotations     // channel (value):
    bcftools_annotations_tbi // channel (value):
    bcftools_header_lines    // string
    snpeff_genome            // string
    snpeff_cache             // channel (value):
    vep_cache                // channel (value):
    vep_cache_version        // string
    vep_extra_files          // array?
    vep_fasta                // channel (value)
    vep_genome               // string
    vep_species              // string

    main:
    reports  = Channel.empty()
    versions = Channel.empty()

    if (tools.split(',').contains('merge') || tools.split(',').contains('bcfann') || tools.split(',').contains('snpeff') || tools.split(',').contains('vep') ) {

        VCF_ANNOTATE_ALL(
            samplesheet.map{meta, vcf -> [ meta + [ file_name: vcf.baseName ], vcf ] },
            vep_fasta,
            tools,
            snpeff_genome,
            snpeff_cache,
            vep_genome,
            vep_species,
            vep_cache_version,
            vep_cache,
            vep_extra_files,
            bcftools_annotations,
            bcftools_annotations_tbi,
            bcftools_header_lines)

        // Gather used softwares versions
        reports  = reports.mix(VCF_ANNOTATE_ALL.out.reports)
        versions = versions.mix(VCF_ANNOTATE_ALL.out.versions)
    }


    emit:
    reports  // channel: [ path(reports) ]
    versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
