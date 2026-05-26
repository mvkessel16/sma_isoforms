Before using MergeReferenceScallop2GFF.R, each 'transcripts_novel.gtf' file was filtered on SMN2 transcripts using:
```
awk '($1 == "chr5" && $4 >= 70808817 && $5 <= 70837732)' /path/to/scallop2_output/${sample}/${sample}_transcripts_novel.gtf  > /path/to/scallop2_output/${sample}/${sample}_transcripts_novel_SMN2.gtf
```
Then, transcript names were extracted using:
```
awk '{print $12}' /path/to/scallop2_output/${sample}/${sample}_transcripts_novel_SMN2.gtf | sed 's/^"\(.*\)";$/\1/' | sort -u > ${sample}/old_transcripts.txt
```

seqkit rename was used to rename the transcripts, so that there is no identical transcript name across the dataset.
