library(dplyr)
library(ggplot2)
library(paletteer)
library(patchwork)
library(ggpubr)
library(readr)

load("/path/to/data")
sample_anon <- read_delim("Anonymization_fibroblasts.txt")

load("/path/to/chm13v2.0_RefSeq_Liftoff_v5.2_SMN1masked_adjusted_gff3.Rdata")
smn2_transcripts <- unique(reference_gtf$target_id[reference_gtf$ens_gene=="SMN2"])
smn2_transcripts <- smn2_transcripts[2:length(smn2_transcripts)]

#Prepare data (whole-cell mRNA dataset)
data_total_SMN2 <- data_total[data_total$target_id %in% smn2_transcripts,]
data_total_SMN2$condition <- rep(c("C","R"))
data_total_SMN2$SMAtype <- c("")
data_total_SMN2 <- data_total_SMN2 %>% group_by(sample,condition) %>% mutate(total_tpm = sum(tpm))
data_total_SMN2$percoftotal <- data_total_SMN2$tpm/data_total_SMN2$total_tpm*100
isoform_names_refseq <- read.table("refseq_isoform_names.txt",header=T)
data_total_SMN2 <- na.omit(full_join(data_total_SMN2, isoform_names_refseq, by = "target_id"))
data_total_SMN2 <- left_join(data_total_SMN2, sample_anon, by = "sample")

#Analyze total SMN2 expression
total_expression <- data_total_SMN2 %>% select(sample, condition, total_tpm, SMAtype)
total_expression <- total_expression[1:32,]
t.test(total_expression$total_tpm[total_expression$condition=="C"],total_expression$total_tpm[total_expression$condition=="R"],paired=T)
total_ris <- ggplot(total_expression, aes(x=condition, y=total_tpm)) +
  geom_boxplot() +
  geom_point(aes(colour=SMAtype)) +
  labs(y="SMN2 expression (TPM)",tag="a",x="treatment",colour="SMA type") +
  scale_x_discrete(labels=c("U","R")) +
  scale_color_paletteer_d("waRhol::camo_87_2")+
  geom_bracket(xmin = "C", xmax = "R", y.position = 135, label = "") +
  theme(plot.title.position = "plot") +
  ylim(0,145)

#Analyze SMN2 isoform expression per sample
barplot_colours <- c(NA,NA,"#4478A0FF","#0F2F4FFF",NA,"#89C9D8FF","#E76254FF", "#EF8A47FF", "#F7AA58FF", "#FFD06FFF", "#FFE6B7FF", "#AADCE0FF", "#72BCD5FF", "#528FADFF", "#376795FF", "#1E466EFF") 
#Paletteer, MetBrewer::Hiroshige
names(barplot_colours) <- c("Δ3,5,7,8_exon9","Δ3,7,8_exon9","Δ5,7,8_exon9","Δ7,8_exon9","exon6b_Δ7,8_v1","Δ3,5,7","exon6b","exon6b_Δ7","exon6b_Δ7,8","FL","Δ3","Δ3,5","Δ3,7","Δ5","Δ5,7","Δ7")
sample_isoform <- data_total_SMN2 %>% filter(condition=="C" & percoftotal>0) %>% ggplot(aes(fill=isoform, y=percoftotal, x=Anonymized_ID)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~SMAtype, scales="free_x",ncol=4) +
  labs(x="ID",y="% of total SMN2 expression",title="a") +
  scale_fill_manual(values=barplot_colours) +
  theme(plot.title.position = "plot") +
  theme(axis.text.x = element_text(angle = 270))

#Analyze expression of isoform between conditions
t.test(SMN2_data$percoftotal[SMN2_data$target_id=="NM_022875.3" & SMN2_data$condition=="C"],SMN2_data$percoftotal[SMN2_data$target_id=="NM_022875.3" & SMN2_data$condition=="R"],paired=T)
d7_ris <- ggplot(data_total_SMN2[data_total_SMN2$target_id=="NM_022875.3",], aes(x=condition,y=percoftotal)) +
  geom_boxplot() +
  geom_point(aes(colour=SMAtype)) +
  ylim(0,119) +
  labs(y="% of total SMN2 expression",tag="d",title="Δ7",x="treatment",colour="SMA type") +
  scale_x_discrete(labels=c("U","R")) +
  scale_color_paletteer_d("waRhol::camo_87_2") +
  geom_bracket(xmin = "C", xmax = "R", y.position = 110, label = "")

d7_ris_tpm <- ggplot(data_total_SMN2[data_total_SMN2$target_id=="NM_022875.3",], aes(x=condition,y=tpm)) +
  geom_boxplot() +
  geom_point(aes(colour=SMAtype)) +
  ylim(0,105) +
  labs(y="Δ7 SMN2 expression (TPM)",tag="b",x="treatment") +
  scale_x_discrete(labels=c("U","R")) +
  scale_color_paletteer_d("waRhol::camo_87_2") +
  geom_bracket(xmin = "C", xmax = "R", y.position = 100, label = "")
