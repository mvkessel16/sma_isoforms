#!/usr/bin/env nextflow

/*
 * Takes a HISAT2-aligned BAM file and sorts it by name
 */
process sortBamName_HISAT2 {
    label 'short'

    publishDir "/path/to/bam_HISAT2/${sample}", mode: 'copy'
    container '/path/to/samtools_1.17--hd87286a_1.img'

    input:
        val sample
        path samfile

    output:
        path "${sample}_HISAT2aligned_name.bam"

    script:
        """
        samtools sort -n -o ${sample}_HISAT2aligned_name.bam ${samfile}
        """
}
