#!/usr/bin/env nextflow

/*
 * Combine reference and novel transcripts
 */
process catTranscriptome {
    label 'short'

    publishDir "/path/to/scallop2_output/${sample}", mode: 'copy'

    input:
        val sample
        path ref_transcripts
        path novel_transcripts

    output:
        path "${sample}_transcripts.fa"

    script:
    """
    cat ${ref_transcripts} ${novel_transcripts} > ${sample}_transcripts.fa
    """
}
