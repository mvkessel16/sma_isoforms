# Analysing SMN2 mRNA isoforms in spinal muscular atrophy using short-read RNA sequencing data
This repository is associated with my MSc internship at UMCU. Here, I analyzed SMN2 transcripts isoforms in SMA. For this purpose, I created a Nextflow pipeline that contains the following main steps:
- Read alignment with HISAT2
- Transcript assembly with Scallop2
- Transcript quantification with Salmon

As input, prepare a CSV file which contains the following parameters per sample:
sample, bam_STAR, read1_unmapped, read2_unmapped

Then run this Bash code:

```
while IFS=',' read -r sample bam_STAR read1_unmapped read2_unmapped 
do
	echo ${sample}    
	sbatch run_main.sh ${sample} ${bam_STAR} ${read1_unmapped} ${read2_unmapped} /path/to/hisat2_index /path/to/ref_gtf /path/to/ref_gen /path/to/ref_transcripts
done < samples.csv > log/jobs.txt
```

Some samples will unfortunately not run, and instead give this error: "ERROR ~ No locks available". These should be run again.

Subsequent analysis is done using Sleuth (https://github.com/pachterlab/sleuth). See R_scripts for these additional scripts.
