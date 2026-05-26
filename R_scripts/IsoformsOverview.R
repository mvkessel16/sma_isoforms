#Load packages
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggtranscript)
library(rtracklayer)
library(tidyr)

##UCSC reference + novel Scallop2 SMN2 isoforms
#Load GFF3 file
load("/path/to/UCSC_T2T_Isoforms_Plus_AllNovelSMN2Isoforms_gff3.RData")
reference_gtf <- total_gtf

#Extract SMN2 isoforms
SMN2_start <- 70808817 #based on SMN2 gene start in ref GFF3 file
SMN2_end <- 70837732 #highest coordinate in 'end' column of SMN2 novel transcripts 
SMN2_isoforms <- dplyr::filter(reference_gtf, start>=SMN2_start & end<=SMN2_end & strand == "-" & seqnames == "chr5")

#Plot
SMN2_exons <- dplyr::filter(SMN2_isoforms, type == "exon")
SMN2_introns <- to_intron(SMN2_exons, "transcript_id")
SMN2_rescaled <- shorten_gaps(exons=SMN2_exons,introns=SMN2_introns,group_var="transcript_id")

#To make the canonical, full-length transcript appear at the bottom of the plot, it is renamed
SMN2_rescaled$transcript_id[SMN2_rescaled$transcript_id=="CHM13_T0175899"]="canonical: CHM13_T0175899"

SMN2_rescaled_exons <- SMN2_rescaled[SMN2_rescaled$type=="exon",]
SMN2_rescaled_introns <- SMN2_rescaled[SMN2_rescaled$type=="intron",]

canonical <- dplyr::filter(SMN2_rescaled_exons, transcript_id == "canonical: CHM13_T0175899")
not_canonical <- dplyr::filter(SMN2_rescaled_exons, transcript_id != "canonical: CHM13_T0175899")

SMN2_rescaled_diffs <- to_diff(
  exons = not_canonical,
  ref_exons = canonical,
  group_var = "transcript_id"
)

plot <- SMN2_rescaled_exons %>%
  ggplot(aes(xstart = start,
             xend = end,
             y = transcript_id
  )) +
  geom_range() +
  geom_intron(
    data = SMN2_rescaled_introns,
    aes(strand=strand),
    arrow.min.intron.length = 300
  ) +
  geom_range(
    data = SMN2_rescaled_diffs,
    aes(fill = diff_type),
    alpha = 0.2
  )+
  labs(fill="Exon compared to SMN2-FL",y="isoform") +
  scale_fill_hue(labels=c("absent","cryptic"))

#Add exon numbers
with_exons <- add_exon_number(canonical, "transcript_id")
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="2","2a",exon_number))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="3","2b",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="4","3",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="5","4",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="6","5",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="7","6",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="8","7",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="9","8",exon))

plot2 <- plot +
  geom_text(
    data = with_exons,
    aes(x = (start + end) / 2, y=Inf, label = exon),
    size = 4) +
  coord_cartesian(clip = "off") +
  theme(axis.title.x=element_blank())

#Save the plot
ggsave("path/to/UCSCPlusAllNovelSMN2Isoforms.png", plot = plot2, width = 210, 
       height = 260, 
       units = "mm")

##RefSeq reference SMN2 isoforms
#Load GFF3 file
load("/path/to/chm13v2.0_RefSeq_Liftoff_v5.2_SMN1masked_adjusted_gff3.Rdata")

#Extract SMN2 isoforms
SMN2_isoforms <- reference_gtf %>% filter(ext_gene=="SMN2")

#Add corresponding ID to each exon and CDS
has_id <- SMN2_isoforms[!is.na(SMN2_isoforms$ID),]
type_has_id <- c(as.vector(unique(has_id$type)),"exon","CDS")
SMN2_isoforms_filled <- SMN2_isoforms %>% fill(ID, .direction = "down")

#Manually give names to isoforms
isoform_names_refseq <- read.table("refseq_isoform_names.txt",header=T)
colnames(isoform_names_refseq) <- c("ID","isoform")
SMN2_isoforms_filled <- full_join(SMN2_isoforms_filled, isoform_names_refseq, by = "ID")

#Plot
SMN2_exons <- dplyr::filter(SMN2_isoforms_filled, type == "exon")
SMN2_introns <- to_intron(SMN2_exons, "isoform")
SMN2_rescaled <- shorten_gaps(exons=SMN2_exons,introns=SMN2_introns,group_var="isoform")

#To make the canonical, full-length transcript appear at the bottom of the plot, it is renamed
SMN2_rescaled$isoform[SMN2_rescaled$ID=="NM_017411.4"]="canonical: FL"

SMN2_rescaled_exons <- SMN2_rescaled[SMN2_rescaled$type=="exon",]
SMN2_rescaled_introns <- SMN2_rescaled[SMN2_rescaled$type=="intron",]

canonical <- dplyr::filter(SMN2_rescaled_exons, ID == "NM_017411.4")
not_canonical <- dplyr::filter(SMN2_rescaled_exons, ID != "NM_017411.4")

SMN2_rescaled_diffs <- to_diff(
  exons = not_canonical,
  ref_exons = canonical,
  group_var = "isoform"
)

plot <- SMN2_rescaled_exons %>%
  ggplot(aes(xstart = start,
             xend = end,
             y = isoform
  )) +
  geom_range() +
  geom_intron(
    data = SMN2_rescaled_introns,
    aes(strand=strand),
    arrow.min.intron.length = 200
  ) +
  geom_range(
    data = SMN2_rescaled_diffs,
    aes(fill = diff_type),
    alpha = 0.2
  ) +
  labs(fill="Exon compared to SMN2-FL") +
  scale_fill_hue(labels=c("absent","cryptic"))

#Add exon numbers
with_exons <- add_exon_number(canonical, "transcript_id")
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="2","2a",exon_number))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="3","2b",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="4","3",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="5","4",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="6","5",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="7","6",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="8","7",exon))
with_exons <- mutate(with_exons, exon=ifelse(exon_number=="9","8",exon))

#Add exon numbers not present in canonical isoform
exon9and6b <- SMN_rescaled_exons[c(58,12),] 
exon9and6b$exon_number <- as.double(exon9and6b$exon_number)
with_exons2 <- add_row(with_exons,exon9and6b)
with_exons2$exon[10] <- as.double("9")
with_exons2$exon[11] <- "6b"

plot2 <- plot +
  geom_text(
    data = with_exons2,
    aes(x = (start + end) / 2, y=Inf, label = exon),
    size = 4) +
  coord_cartesian(clip = "off") +
  theme(axis.title.x=element_blank())

#Save the plot
ggsave("/path/to/isoformsRefSeq.png", plot = plot2, width = 210, 
       height = 130, 
       units = "mm")
