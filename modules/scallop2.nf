#!/usr/bin/env nextflow

/*
 * Assembles transcripts from BAM file using Scallop2
 */
process scallop2 {
    label 'short'

    publishDir "/path/to/scallop2_output/${sample}", mode: 'copy'
    container '/path/to/scallop2-1.1.2--h5ca1c30_8.img'

    input:
        val sample
        path bamfile

    output:
        path "${sample}_transcripts_scallop2.gtf"

    script:
    """
    scallop2 -i ${bamfile} -o ${sample}_transcripts_scallop2.gtf
    """
}
