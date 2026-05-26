#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -o log/slurm_nextflow.%j.out
#SBATCH -e log/slurm_nextflow.%j.err
#SBATCH --mail-user "email@email.com"
#SBATCH --mail-type FAIL

sample=$1
bam_STAR=$2
read1_unmapped=$3
read2_unmapped=$4
hisat2_index=$5
ref_gtf=$6
ref_gen=$7
ref_transcripts=$8

nextflow run mainpipeline_poly.nf --sample ${sample} --bam_STAR ${bam_STAR} \
    --read1_unmapped ${read1_unmapped} --read2_unmapped ${read2_unmapped} \
    --hisat2_index ${hisat2_index} --ref_gtf ${ref_gtf} --ref_gen ${ref_gen} \
    --ref_transcripts ${ref_transcripts}
