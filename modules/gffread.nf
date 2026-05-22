#!/usr/bin/env nextflow

/*
 * Extract FASTA sequences of transcripts from GTF file
 */
process gffread {
    label 'short'

    publishDir "/path/to/scallop2_output/${sample}", mode: 'copy'
    container "/path/to/gffread-0.12.7--h077b44d_6.img"

    input:
        val sample
        path scallop2_gtf
        path ref_gen

    output:
        path "${sample}_transcripts_novel.fa"

    script:
    """
    gffread ${scallop2_gtf} -g ${ref_gen} -w ${sample}_transcripts_novel.fa
    """
}
