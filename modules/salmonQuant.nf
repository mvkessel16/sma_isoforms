#!/usr/bin/env nextflow

/*
 * Estimates the expression of each transcript using Salmon
 */
process salmonQuant {
    label 'long'

    publishDir "/path/to/salmon_output/${sample}", mode: 'copy'
    container '/path/to/salmon-1.10.3--h45fbf2d_5.img'

    input:
        val sample
        path salmon_index
        path read1
        path read2

    output:
        path "${sample}_salmonquant/*"

    script:
    """
    salmon quant -l A --numBootstraps 30 -i ${salmon_index} -1 ${read1} -2 ${read2} --seqBias --gcBias -o ${sample}_salmonquant
    """
}
