#!/usr/bin/env nextflow

/*
 * Marks duplicate reads
 */
process sambambaMarkDup {
    label 'short'

    publishDir "/path/to/sambamba/${sample}", mode: 'copy'
    container '/path/to/sambamba_1.0--h98b6b92_0.img'

    input:
        val sample
        path bamfile

    output:
        path "${sample}_markdup.bam"

    script:
    """
    sambamba markdup -t 4 --tmpdir . ${bamfile} ${sample}_markdup.bam
    """
}
