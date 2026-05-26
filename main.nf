#!/usr/bin/env nextflow

/*
 * Pipeline parameters: sample, bam_STAR, read1_unmapped, read2_unmapped, hisat2_index, ref_gtf, ref_gen, ref_transcripts
 */

// Include modules
include { sortBamName_STAR } from './modules/sortBamName_STAR.nf'
include { bamToFastq_STAR } from './modules/bamToFastq_STAR.nf'
include { hisat2 } from './modules/hisat2.nf'
include { sortBamCoord } from './modules/sortBamCoord.nf'
include { indexBam } from './modules/indexBam.nf'
include { scallop2 } from './modules/scallop2.nf'
include { gffcompare } from './modules/gffcompare.nf'
include { gtfcuff } from './modules/gtfcuff.nf'
include { gffread } from './modules/gffread.nf'
include { catTranscriptome } from './modules/catTranscriptome.nf'
include { salmonIndex } from './modules/salmonIndex.nf'
include { sortBamName_HISAT2 } from './modules/sortBamName_HISAT2.nf'
include { bamToFastq_HISAT2 } from './modules/bamToFastq_HISAT2.nf'
include { salmonQuant } from './modules/salmonQuant.nf'
include { sambambaMarkDup } from './modules/sambambaMarkDup.nf'
include { sambambaFlagstat } from './modules/sambambaFlagstat.nf'


workflow {

    // sort STAR-aligned file from facility by name
    sortBamName_STAR(params.sample,params.bam_STAR)

    // convert to FASTQ files
    bamToFastq_STAR(params.sample,sortBamName_STAR.out)

    // align with HISAT2
    hisat2(params.sample, bamToFastq_STAR.out.read1, params.read1_unmapped, bamToFastq_STAR.out.read2, params.read2_unmapped)

    // sort HISAT2-aligned file by coord
    sortBamCoord(params.sample,hisat2.out.sam)

    // index BAM file
    indexBam(params.sample,sortBamCoord.out)

    // assemble transcripts with Scallop2
    scallop2(params.sample,sortBamCoord.out)

    // compare Scallop2 GTF with reference GTF
    gffcompare(params.sample,scallop2.out,params.ref_gtf)

    // extract novel transcripts
    gtfcuff(params.sample,scallop2.out,gffcompare.out.tmap)

    // extract FASTA of novel transcripts
    gffread(params.sample,gtfcuff.out,params.ref_gen)

    // create transcriptome
    catTranscriptome(params.sample,params.ref_transcripts,gffread.out)

    // create Salmon index
    salmonIndex(params.sample,catTranscriptome.out)

    // sort HISAT2-aligned file by name
    sortBamName_HISAT2(params.sample,hisat2.out.sam)

    // convert to FASTQ files
    bamToFastq_HISAT2(params.sample,sortBamName_HISAT2.out)

    // quantify transcript expression
    salmonQuant(params.sample,salmonIndex.out,bamToFastq_HISAT2.out.read1,bamToFastq_HISAT2.out.read2)

    // mark duplicate reads
    sambambaMarkDup(params.sample,sortBamCoord.out)

    // get FLAG statistics
    sambambaFlagstat(params.sample,sambambaMarkDup.out)
    
}
