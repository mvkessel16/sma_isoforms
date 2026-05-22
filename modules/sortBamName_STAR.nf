#!/usr/bin/env nextflow

/*
 * Takes a STAR-aligned BAM file and sorts it by name
 */
process sortBamName_STAR {
    label 'short'

    publishDir "/path/to/bam_STAR/${sample}", mode: 'copy'
    container '/path/to/samtools_1.17--hd87286a_1.img'

    input:
        val sample
        path bamfile

    output:
        path "${sample}_STARaligned_name.bam"

    script:
        """
        samtools sort -n -o ${sample}_STARaligned_name.bam ${bamfile}
        """
}
