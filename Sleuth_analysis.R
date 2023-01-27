#### Waterhemp Sleuth analysis #####

library(tidyverse)
library(sleuth)

# load metadata

metadata <- read_csv("waterhemp_metadata.csv")

# sample id

sample_id <- metadata$sample


# get kallisto results

kal_dirs <- file.path("../kallisto/waterhemp/kallisto_out",  sample_id)

# join pathway to metadata

metadata <- metadata %>% 
  mutate(path = kal_dirs)

# get gene to transcript association

txi2gene <- read.table("../trinity_assembly/trinity_waterhemp/trinity_out_dir/Trinity-GG.fasta.gene_trans_map")
colnames(txi2gene) <- c("gene_id","target_id")
txi2gene <- txi2gene[,c('target_id','gene_id')]

# get annotation 

blastp_annot <- read_csv("../trinity_assembly/trinity_waterhemp/annotation/blastP_annotation.csv")
blastp_annot_gene <- blastp_annot %>% 
  distinct(gene, .keep_all = T) %>% 
  select(-transcript, -protein)
# Sleuth prep

so <- sleuth_prep(metadata, target_mapping = txi2gene,
                  aggregation_column = 'gene_id', 
                  extra_bootstrap_summary = TRUE)


# fit models

so <- sleuth_fit(so, ~tissue, 'reduced')
so <- sleuth_fit(so, ~ gender + tissue, 'full')

so <- sleuth_lrt(so, 'reduced', 'full')

#gene level
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
gene_level_results <- dplyr::filter(sleuth_table_gene, qval <= 0.05)


head(gene_level_results)

gene_level_results <- tibble(
  gene_level_results %>%
  rename(gene = target_id) %>% 
  left_join(blastp_annot_gene)
)

View(gene_level_results %>% 
       filter(!is.na(uniprot)))

#transcript level
sleuth_table_tx <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
tx_results <- dplyr::filter(sleuth_table_tx, qval <= 0.05)
head(tx_results)
tx_results <-tibble(
  tx_results %>%
    rename(transcript = target_id, gene = gene_id) %>% 
    left_join(blastp_annot)
)

tx_results <- tx_results %>% 
  filter(!is.na(uniprot))

sleuth_live(so)
