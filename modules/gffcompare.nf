#!/usr/bin/env nextflow

/*
 * Compares Scallop2 GTF with reference GTF
 */
process gffcompare {
    label 'short'

    publishDir "/path/to/scallop2_output/${sample}", mode: 'copy'
    container "/path/to/gffcompare-0.12.10--h9948957_0.img"

    input:
        val sample
        path scallop2_gtf
        path ref_gtf

    output:
        path "${sample}_gffcompare.${scallop2_gtf}.tmap", emit: tmap
        path "${sample}_gffcompare.*"

    script:
    """
    gffcompare -o ${sample}_gffcompare -r ${ref_gtf} ${scallop2_gtf}
    """
}
