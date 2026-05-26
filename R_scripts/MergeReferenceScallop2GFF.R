#Load packages
library(dplyr)
library(rtracklayer)

#Rename isoforms in Scallop2 files
sample_names <- readLines("samples.txt")
scallop_files <- list.files(path="/path/to/scallop2_output/",pattern="novel_SMN2.gtf",full.names=T,recursive=T)
old_transcript_files <- list.files(path="/path/to/scallop2_output/",pattern="old_transcripts.txt",full.names=T,recursive=T)
fasta_files <- list.files(path="/path/to/scallop2_output/renamed_novelIsoforms_fa/",pattern="_transcripts_novel.fa",full.names=T)
names_gtf <- NULL
new_gene_name <- "CHM13_G0044653" #This is the gene SMN2

for (i in 1:32){
  sample <- sample_names[i]
  #If there are new SMN2 isoforms found by Scallop2
  if (file.size(scallop_files[i])>0){ 
    gtf_novel <- rtracklayer::import(scallop_files[i])
    gtf_novel <- dplyr::as_tibble(gtf_novel)
    old_transcripts <- readLines(old_transcript_files[i])
    fasta <- readLines(fasta_files[i])
    for (j in (1:length(old_transcripts))){
      new_transcriptname <- fasta[grepl(old_transcripts[j],fasta)]
      new_transcriptname <- sub("^>", "", new_transcript)
      gtf_novel <- gtf_novel %>% mutate(transcript_id = ifelse(transcript_id==old_transcripts[j], new_transcriptname, transcript_id))
      gtf_novel <- gtf_novel %>% mutate(gene_id = new_gene_name)
    }
    name_gtf <- paste0('/path/to/scallop2_output/renamed_novelSMN2Isoforms_gtf/',sample,'_changed_novel_isoforms.gtf')
    names_gtf <- c(names_gtf,name_gtf) #Collect names of renamed Scallop2 GFF files
    rtracklayer::export.gff2(gtf_novel,con=name_gtf) 
  }
}  

#Create one file with all novel SMN2 isoforms identified by Scallop2
all_novel_gtf <- NULL

for (i in 1:25){ #The names_gtf vector was 25 items long
  gtf <- readLines(names_gtf[i])
  gtf <- gtf[4:length(gtf)]
  all_novel_gtf <- c(all_novel_gtf,gtf)
}

writeLines(all_novel_gtf, "/path/to/AllNovelSMN2Isoforms.gtf")

#Combine UCSC reference with novel isoforms in GFF3
novel_gtf <- rtracklayer::import("/path/to/AllNovelSMN2Isoforms.gtf")
novel_gtf <- dplyr::as_tibble(novel_gtf)
ref_gtf <- rtracklayer::import("path/to/T2T_UCSC_GENCODEv35_CAT_Liftoffv2_SMN1masked.gff3")
ref_gtf <- dplyr::as_tibble(ref_gtf)
total_gtf <- bind_rows(ref_gtf,novel_gtf)
total_gtf$gene_name[total_gtf$source=="scallop2"]="SMN2"
rtracklayer::export.gff3(total_gtf,con="/path/to/UCSC_T2T_Isoforms_Plus_AllNovelSMN2Isoforms.gff3")
save(total_gtf,file="/path/to/UCSC_T2T_Isoforms_Plus_AllNovelSMN2Isoforms_gff3.RData")

#Combine UCSC reference with novel isoforms in GFF2/GTF
ref_gtf2 <- rtracklayer::import.gff2("/path/to/T2T_UCSC_GENCODEv35_CAT_Liftoffv2_SMN1masked.gtf")
ref_gtf2 <- dplyr::as_tibble(ref_gtf2)
total_gtf2 <- bind_rows(ref_gtf2,novel_gtf)
total_gtf2$gene_name[total_gtf2$source=="scallop2"]="SMN2"
rtracklayer::export.gff2(total_gtf2,con="/path/to/UCSC_T2T_Isoforms_Plus_AllNovelSMN2Isoforms.gtf")
