#!/usr/bin/env nextflow

/*
 * Takes a BAM file and sorts it by coordinates
 */
process sortBamCoord {
    label 'short'

    publishDir "/path/to/hisat2_output/${sample}", mode: 'copy'
    container '/path/to/samtools_1.17--hd87286a_1.img'

    input:
        val sample
        path samfile

    output:
        path "${sample}_HISAT2aligned_coord.bam"

    script:
        """
        samtools sort -o ${sample}_HISAT2aligned_coord.bam ${samfile} 
        """
}
