#!/usr/bin/env nextflow

/*
 * Creates file with statistics about FLAG tags
 */
process sambambaFlagstat {
    label 'short'

    publishDir "/path/to/sambamba/${sample}", mode: 'copy'
    container '/path/to/sambamba_1.0--h98b6b92_0.img'

    input:
        val sample
        path bamfile

    output:
        path "${sample}_flagstat.txt"

    script:
    """
    sambamba flagstat -t 4 ${bamfile} > ${sample}_flagstat.txt
    """
}
