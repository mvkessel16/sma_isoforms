library(tidyverse)
library(sleuth)
library(wasabi)
library(annotables)
library(rtracklayer)
library(dplyr)

#Prepare the data for Sleuth
quant_dirs <- list.dirs(path="/path/to/salmon_output")
prepare_fish_for_sleuth(quant_dirs)

#Create metadata
samples <- readLines("samples.txt")
sample_groups <- data.frame(sampleID=samples,row.names=T)
sample_groups$sample <- row.names(sample_groups)
sample_groups$path <- quant_dirs
sample_groups$condition <- c("U","R")
sample_groups$SMAtype <- c("")
sample_groups$sex <- c("")
sample_groups$ageatbiopsy <- c("")
sample_groups$preferredtreatment <- c("")

names(quant_dirs) <- samples

reference_gtf <- rtracklayer::import("/path/to/chm13v2.0_RefSeq_Liftoff_v5.2_SMN1masked_adjusted.gff3")
reference_gtf <- reference_gtf %>% dplyr::as_tibble()
reference_gtf <- reference_gtf %>% fill(gene_name, .direction = "down")
reference_gtf <- reference_gtf %>% dplyr::rename(ext_gene = gene_name)
reference_gtf <- reference_gtf %>% mutate(ens_gene = NA)
reference_gtf <- reference_gtf %>% mutate(ens_gene=ifelse(type == "gene", ID, NA))
reference_gtf <- reference_gtf %>% fill(ens_gene, .direction = "down")
reference_gtf <- reference_gtf %>% mutate(target_id=ifelse(type == "transcript", ID, NA))
reference_gtf <- reference_gtf %>% fill(target_id, .direction = "down") %>% mutate(target_id = ifelse(type!="gene", target_id, NA))

info <- reference_gtf[reference_gtf$type=="transcript",] %>% 
  select(ens_gene,ext_gene,target_id)

#Create Sleuth object for analysis (this normalizes the data)
so_transcript <- sleuth_prep(sample_groups, 
                  target_mapping = info, 
                  read_bootstrap_tpm = TRUE,
                  extra_bootstrap_summary = TRUE,
                  transformation_function = function(x) log2(x + 0.5))

smn2_transcripts <- unique(reference_gtf$target_id[reference_gtf$ens_gene=="SMN2"])
smn2_transcripts <- smn2_transcripts[2:length(smn2_transcripts)]

#Extract the normalized data
data <- so_transcript[["obs_norm"]]
SMN2_data <- data[data$target_id %in% smn2_transcripts,]

#Statistical tests between conditions
design <- ~ ID + condition

so_transcript <- sleuth_fit(so_transcript,formula=design)

models(so_transcript)

oe <- sleuth_wt(so_transcript, 
                which_beta = 'conditionR')

sleuth_results_oe_R <- sleuth_results(oe, 
                                    test = 'conditionR', 
                                    show_all = TRUE)

sleuth_results_SMN2 <-  sleuth_results_oe[sleuth_results_oe$ens_gene=="SMN2",]

#Heatmap
smn2_transcripts_exp <- SMN2_data %>% filter(tpm != 0) %>% select(target_id)
smn2_transcripts_exp <- unique(smn2_transcripts_exp$target_id)
plot_transcript_heatmap(oe, 
                        transcripts = smn2_transcripts_exp)

##To analyze differences between SMA types
sample_groups_unt <- sample_groups %>% filter(condition=="C")

so_refseq_unt_SMAsamples <- sleuth_prep(sample_groups_unt, 
                             target_mapping = info, 
                             read_bootstrap_tpm = TRUE,
                             extra_bootstrap_summary = TRUE,
                             transformation_function = function(x) log2(x + 0.5))

sample_groups_unt$SMAtype <- factor(sample_groups_unt$SMAtype)
sample_groups_unt$SMAtype <- relevel(sample_groups_unt$SMAtype, ref = "2a")
levels(sample_groups_unt$SMAtype)
so_refseq_unt_SMAsamples_2abase <- sleuth_prep(sample_groups_unt, 
                                        target_mapping = info, 
                                        read_bootstrap_tpm = TRUE,
                                        extra_bootstrap_summary = TRUE,
                                        transformation_function = function(x) log2(x + 0.5))

sample_groups_unt <- sample_groups %>% filter(condition=="C")
sample_groups_unt$SMAtype <- factor(sample_groups_unt$SMAtype)
sample_groups_unt$SMAtype <- relevel(sample_groups_unt$SMAtype, ref = "2b")
levels(sample_groups_unt$SMAtype)
so_refseq_unt_SMAsamples_2bbase <- sleuth_prep(sample_groups_unt, 
                                               target_mapping = info, 
                                               read_bootstrap_tpm = TRUE,
                                               extra_bootstrap_summary = TRUE,
                                               transformation_function = function(x) log2(x + 0.5))

design <- ~ SMAtype
so_refseq_unt_SMAsamples <- sleuth_fit(so_refseq_unt_SMAsamples,formula=design)
so_refseq_unt_SMAsamples_2abase <- sleuth_fit(so_refseq_unt_SMAsamples_2abase,formula=design)
so_refseq_unt_SMAsamples_2bbase <- sleuth_fit(so_refseq_unt_SMAsamples_2bbase,formula=design)

oe_refseq_unt_SMAsamples_2a <- sleuth_wt(so_refseq_unt_SMAsamples, 
                                         which_beta = 'SMAtype2a')
sleuth_results_oe_2a <- sleuth_results(oe_refseq_unt_SMAsamples_2a,
                                       test="SMAtype2a")

oe_refseq_unt_SMAsamples_2b <- sleuth_wt(so_refseq_unt_SMAsamples, 
                                         which_beta = 'SMAtype2b')
sleuth_results_oe_2b <- sleuth_results(oe_refseq_unt_SMAsamples_2b,
                                       test="SMAtype2b")

oe_refseq_unt_SMAsamples_3a <- sleuth_wt(so_refseq_unt_SMAsamples, 
                                         which_beta = 'SMAtype3a')
sleuth_results_oe_3a <- sleuth_results(oe_refseq_unt_SMAsamples_3a,
                                       test="SMAtype3a")

oe_refseq_unt_SMAsamples_2a_2b <- sleuth_wt(so_refseq_unt_SMAsamples_2abase, 
                                         which_beta = 'SMAtype2b')
sleuth_results_oe_2a_2b <- sleuth_results(oe_refseq_unt_SMAsamples_2a_2b,
                                       test="SMAtype2b")

oe_refseq_unt_SMAsamples_2a_3a <- sleuth_wt(so_refseq_unt_SMAsamples_2abase, 
                                         which_beta = 'SMAtype3a')
sleuth_results_oe__2a_3a <- sleuth_results(oe_refseq_unt_SMAsamples_2a_3a,
                                       test="SMAtype3a")

oe_refseq_unt_SMAsamples_2b_3a <- sleuth_wt(so_refseq_unt_SMAsamples_2bbase, 
                                            which_beta = 'SMAtype3a')
sleuth_results_oe__2b_3a <- sleuth_results(oe_refseq_unt_SMAsamples_2b_3a,
                                           test="SMAtype3a")

##To analyze differences between SMA types when treated with risdiplam
sample_groups_ris <- sample_groups %>% filter(condition=="R")

so_refseq_ris_SMAtotal <- sleuth_prep(sample_groups_ris, 
                                        target_mapping = info, 
                                        read_bootstrap_tpm = TRUE,
                                        extra_bootstrap_summary = TRUE,
                                        transformation_function = function(x) log2(x + 0.5))

design <- ~ SMAtype
so_refseq_ris_SMAtotal <- sleuth_fit(so_refseq_ris_SMAtotal,formula=design)

oe_refseq_ris_SMAtotal_3a <- sleuth_wt(so_refseq_ris_SMAtotal, 
                                         which_beta = 'SMAtype3a')
oe_refseq_ris_SMAtotal_3a_results <- sleuth_results(oe_refseq_ris_SMAtotal_3a,
                                       test="SMAtype3a")
