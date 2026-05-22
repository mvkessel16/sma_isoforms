#!/usr/bin/env nextflow

/*
 * Runs the script read_duplication.py from RSeQC on the HISAT2-aligned BAM files
 */
process RSeQC {
    label 'short'

    publishDir "/path/to/RSeQC/${sample}/", mode: 'copy'
    container "/path/to/rseqc_5.0.4--pyhdfd78af_1.img"

    input:
        val sample
        path bam

    output:
        path "${sample}_read_duplication*"

    script:
    """
    read_duplication.py -i ${bam} -o ${sample}_read_duplication
    """
}
