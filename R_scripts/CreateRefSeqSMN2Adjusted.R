library(dplyr)
library(tidyverse)
library(rtracklayer)
load("/path/to/chm13v2.0_RefSeq_Liftoff_v5.2_SMN1masked_gff3.Rdata")

##Adjust all 5' UTR to major TSS as in FL and d7
reference_gtf <- reference_gtf %>% mutate(end=ifelse(gene=="SMN2" & end==70837821,70837675,end))

##Delete one RefSeq 6b isoform
reference_gtf <- reference_gtf %>% filter_out(Parent=="XM_017009787.2")
reference_gtf <- reference_gtf %>% filter_out(ID=="XM_017009787.2")

#Create isoforms from Singh et al.
singh_isoforms <- reference_gtf %>% filter(Parent=="NM_017411.4" | ID=="NM_017411.4")
singh_isoforms <- singh_isoforms[rep(1:nrow(singh_isoforms), 2),]
singh_isoforms$source <- "Singh"
singh_isoforms <- singh_isoforms %>% mutate(Parent=ifelse(singh_isoforms$type=="transcript",Parent,as.list("Singh_6b_7_8")))
singh_isoforms <- singh_isoforms %>% mutate(ID=ifelse(singh_isoforms$type=="transcript","Singh_6b_7_8",ID))
singh_isoforms$product <- NA
singh_isoforms$tag <- NA
singh_isoforms$note <- NA
singh_isoforms$matches_ref_protein <- NA
singh_isoforms$protein_id <- NA
singh_isoforms$valid_ORF <- NA
singh_isoforms$Parent[20:nrow(singh_isoforms)] <- as.list("Singh_6b_8")
singh_isoforms$ID[19] <- "Singh_6b_8"
singh_isoforms$exon_number[c(2,20)] <- "10"
singh_isoforms$exon_number[c(3,11,21,29)] <- "9"

exonrow <- singh_isoforms[3,]
exonrow$end <- 70816056
exonrow$start <- 70815949
exonrow$width <- 108
exonrow$source <- "Singh"
exonrow$Parent[1] <- as.list("Singh_6b_7_8")
exonrow$product <- NA
exonrow$exon_number <- "8"
exonrow$tag <- NA
exonrow <- exonrow[rep(1, 2),]
exonrow$Parent[2] <- as.list("Singh_6b_8")

cdsrow <- singh_isoforms[12,]
cdsrow$end <- 70816056
cdsrow$start <- 70815949
cdsrow$width <- 108
cdsrow$source <- "Singh"
cdsrow$phase <- NA
cdsrow$Parent[1] <- as.list("Singh_6b_7_8")
cdsrow$product <- NA
cdsrow$exon_number <- "8"
cdsrow$protein_id <- NA
cdsrow$tag <- NA
cdsrow$note <- NA
cdsrow <- cdsrow[rep(1, 2),]
cdsrow$Parent[2] <- as.list("Singh_6b_8")
cdsrow$phase <- 0

singh_isoforms <- add_row(singh_isoforms,exonrow[1,],.before=4)
singh_isoforms <- add_row(singh_isoforms,exonrow[2,],.before=23)
singh_isoforms <- add_row(singh_isoforms,cdsrow[1,],.before=13)
singh_isoforms <- add_row(singh_isoforms,cdsrow[2,],.before=33)
singh_isoforms <- singh_isoforms[-c(23,32),]
singh_isoforms$exon_number[22] <- "9"
singh_isoforms$target_id[1:20] <- "Singh_6b_7_8"
singh_isoforms$target_id[21:38] <- "Singh_6b_8"

#Add Singh isoforms to reference
reference_gtf <- add_row(reference_gtf,singh_isoforms,.before=1240347)
rtracklayer::export.gff3(reference_gtf,con="/path/to/chm13v2.0_RefSeq_Liftoff_v5.2_SMN1masked_adjusted.gff3")
save(reference_gtf,file="/path/to/chm13v2.0_RefSeq_Liftoff_v5.2_SMN1masked_adjusted_gff3.Rdata")
