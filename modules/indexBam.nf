#!/usr/bin/env nextflow

/*
 * Makes an index of a BAM file
 */
process indexBam {
    label 'short'

    publishDir "/path/to/hisat2_output/${sample}", mode: 'copy'
    container '/path/to/samtools_1.17--hd87286a_1.img'

    input:
        val sample
        path bamfile

    output:
        path "${bamfile}.bai"

    script:
        """
        samtools index ${bamfile} ${bamfile}.bai
        """
}
