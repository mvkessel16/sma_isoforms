#!/usr/bin/env nextflow

/*
 * Takes a name-sorted STAR-aligned BAM file and converts it to FASTQ files
 */
process bamToFastq_STAR {
    label 'short'

    publishDir "path/to/fastq_STAR/${sample}", mode: 'copy'
    container 'path/to/samtools_1.17--hd87286a_1.img'

    input:
        val sample
        path bamfile

    output:
        path "${sample}_STAR_single_bam_to.fastq"
        path "${sample}_STAR_read1_bam_to.fastq", emit: read1
        path "${sample}_STAR_read2_bam_to.fastq", emit: read2

    script:
    """
    samtools fastq -n -s ${sample}_STAR_single_bam_to.fastq -1 ${sample}_STAR_read1_bam_to.fastq -2 ${sample}_STAR_read2_bam_to.fastq ${bamfile}
    """
}
