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
    ch_samplesheet // channel: samplesheet read in from --input
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
    ch_reports  = Channel.empty()
    ch_versions = Channel.empty()

    if (tools.split(',').contains('merge') || tools.split(',').contains('bcfann') || tools.split(',').contains('snpeff') || tools.split(',').contains('vep') ) {

        VCF_ANNOTATE_ALL(
            ch_samplesheet.map{meta, vcf -> [ meta + [ file_name: vcf.baseName ], vcf ] },
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
        ch_reports  = ch_reports.mix(VCF_ANNOTATE_ALL.out.reports)
        ch_versions = ch_versions.mix(VCF_ANNOTATE_ALL.out.versions)
    }


    emit:
    reports  = ch_reports  // channel: [ path(reports) ]
    versions = ch_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
