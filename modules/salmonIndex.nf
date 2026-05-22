#!/usr/bin/env nextflow

/*
 * Makes a Salmon index of a transcriptome FASTA
 */
process salmonIndex {
    label 'short'

    publishDir "/path/to/salmon_output/${sample}", mode: 'copy'
    container '/path/to/salmon-1.10.3--h45fbf2d_5.img'

    input:
        val sample
        path transcripts

    output:
        path "${sample}_salmonindex/"

    script:
    """
    salmon index -t ${transcripts} -i ${sample}_salmonindex
    """
}
