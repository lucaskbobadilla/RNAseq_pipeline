# EDGER Waterhemp analysis ####################


# packages ===================================
library(tidyverse)
library(tximport)
library(edgeR)
library(PCAtools)

# load metadata ========================

sample_info <- read_csv("waterhemp/data/waterhemp_metadata.csv")

sample_info <- data.frame(sample_info)

rownames(sample_info) <- sample_info$sample

sample_info <- sample_info[-1]

sample_info <- sample_info %>% 
  mutate(groups = factor(paste(gender,tissue,sep = ".")))

## Import Transcriptome annotation ==========================================


# blastP
annot_blastp <- read_csv("../trinity_assembly/trinity_waterhemp/annotation/blastP_annotation.csv")

annot_blastp <- annot_blastp %>%
  group_by(gene) %>% 
  slice_max(identity) %>%
  distinct(gene, .keep_all = T) %>% 
  dplyr::select(gene, uniprot, identity, blastp_eval = evalue,blastP_annot = name, lineage,genus) %>% 
  ungroup()


# Pfam
annot_pfam <- read_csv("../trinity_assembly/trinity_waterhemp/annotation/pfam_annotation.csv")

annot_pfam <- annot_pfam %>%
  group_by(gene) %>% 
  slice_min(evalue) %>%
  distinct(gene, .keep_all = T) %>% 
  dplyr::select(gene, pfam, pfam_symbol = symbol,Pfam_annot = name, pfam_eval = evalue) %>% 
  ungroup()

waterhemp_annot <- annot_blastp %>% 
  full_join(annot_pfam, by = "gene")


# load gene location

gff <- ape::read.gff("../trinity_assembly/trinity_waterhemp/trinity_out_dir/trinity_gene_Male_WH.gff")

gene_location <- gff |> 
  filter(type == "mRNA") |> 
  separate(attributes, into = c('ID', 'Name','Parent', 'Dir', 'Coverage', "Identity", "Matches",
                                'Indels', "Unknowns"), sep = ";") |> 
  mutate_all(sub, pattern = ".*=", replacement = '') |> 
  dplyr::select(Name, scaffold = seqid, start, end, Coverage:Unknowns) |> 
  tibble() |> 
  arrange(Name, desc(Coverage), desc(Identity))

gene_location <- gene_location %>% 
  separate(scaffold, into = c('scaffold', 'arrow'), sep = "\\|") %>% 
  dplyr::select(-arrow)

write_csv(gene_location, "./waterhemp/data/gene_WH_location.csv")


View(gene_location |> filter(Name == 'TRINITY_GG_10910_c469_g1'))

waterhemp_annot <- waterhemp_annot %>% 
  left_join(gene_location |> rename(gene = Name)) %>% 
  mutate(fac = paste(scaffold,gene, sep = ".")) %>% 
  distinct(fac, .keep_all = T) %>% 
  dplyr::select(-fac)
## Import count files - gene level ==========================================

# get transcript to gene information
tx2gene <- read.table("../trinity_assembly/trinity_waterhemp/trinity_out_dir/Trinity-GG.fasta.gene_trans_map")
colnames(tx2gene) <- c("gene_id", "transcript_id")
tx2gene <- tx2gene[c("transcript_id", "gene_id")]

# get file pathway
files <- file.path("../kallisto/waterhemp/kallisto_out",  rownames(sample_info), 
                   "abundance.h5")
names(files) <- rownames(sample_info)

# import read quantification
txi.rsem <- tximport(files, type = "kallisto", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene)


### get count data ==============================
#counts <- read.table("waterhemp/data/raw_counts.txt", sep = "\t") # for raw counts
counts <- read.table("waterhemp/data/annotated_counts.txt") # for annotated plant genes


waterhemp_annot_genes <-  waterhemp_annot %>% 
  distinct(gene,.keep_all = T)

counts <- data.frame(counts) %>% filter(rownames(counts) %in% waterhemp_annot$gene) |> as.matrix()


counts <- data.frame(txi.rsem$counts)
# counts[1:10,1:10]
# write.table(counts, "waterhemp/data/raw_counts.txt", 
#              sep = "\t", quote = F, 
#              col.names = T, row.names = T)


#counts <- subset(counts, rownames(counts) %in% waterhemp_annot$gene)

rm(txi.rsem) # remove tximport file to save RAM



##  Create DGE list and normalization ================================

dgList <- DGEList(counts=counts, genes=rownames(counts), group = sample_info$groups)
?DGEList
colnames(counts)
### Make a design ==============
design <- model.matrix(~0+groups, sample_info)
colnames(design) <- levels(sample_info$groups)
### Make contrasts ========================================

levels(sample_info$groups)

contrasts <- makeContrasts(
  # Comparison across genders
  MaleFlowerVsFemaleFlower = Male.Flower-Female.Flower,
  MaleFMvsFemaleFM = Male.Flower_meristem -Female.Flower_meristem,
  MaleSAMvsFemaleSAM = Male.SA_meristem - Female.SA_meristem,
  MalevsFemale = (Male.Flower + Male.Flower_meristem + Male.SA_meristem)-(Female.Flower + Female.Flower_meristem + Female.SA_meristem),
  # Comparison across tissues within each gender
  Male_SAMvsFM = Male.SA_meristem - Male.Flower_meristem,
  Male_SAMvsFlower = Male.SA_meristem - Male.Flower,
  Male_FMvsFlower = Male.Flower_meristem - Male.Flower,
  Female_SAMvsFM = Female.SA_meristem - Female.Flower_meristem,
  Female_SAMvsFlower = Female.SA_meristem - Female.Flower,
  Female_FMvsFlower = Female.Flower_meristem - Female.Flower, 
                           levels = design)


### Filter low expressed genes =======================

#get counts per million and remove counts that are less than 1
#countsPerMillion <- cpm(dgList) #get counts
#countCheck <- countsPerMillion > 1 # create a conditional matrix for values > 1
#keep <- which(rowSums(countCheck) >= 1) #filter values
keep <- filterByExpr(dgList)
table(keep)
#dgList <- dgList[keep,]
dgList <- dgList[keep, , keep.lib.sizes=FALSE]



### normalization =======================================

dgList <- calcNormFactors(dgList, method="TMM")
dgList$samples
dgList$counts[1:5,1:5]

?calcNormFactors
### Dispersion =============================================

# apply robust method for QL model
dgList <- estimateDisp(dgList)
plotBCV(dgList)


### Data exploration MDS plot =================================

pdf("waterhemp/results/MDS_plot.pdf")
points <- c(15,16,17,18,1,2)
colors <- c("blue", "darkgreen", "red", "purple","black","brown")
plotMDS(dgList, col=colors[sample_info$groups], pch=points[sample_info$groups])
legend("topleft", legend=levels(sample_info$groups), pch=points, col=colors, ncol=2,cex=0.8)
dev.off()

#### DEG analysis ========================

# fit QL model
fit <- glmQLFit(dgList, design, robust = TRUE)





# make QL tests ----------------------------

QLtests <- list()
for (i in colnames(contrasts)) {
  qlf <- glmQLFTest(fit,contrast = contrasts[,i])
  QLtests[[i]] <- qlf
}

# qlf_flower <- glmQLFTest(fit,contrast = contrasts[,'MaleFlowerVsFemaleFlower'])
# qlf_FM <- glmQLFTest(fit,contrast = contrasts[,'MaleFMvsFemaleFM'])
# qlf_SAM <- glmQLFTest(fit,contrast = contrasts[,'MaleSAMvsFemaleSAM'])
# qlf_MaleVsFemale <- glmQLFTest(fit,contrast = contrasts[,'MalevsFemale'])

# DE tables ------------


DE_tables <- list()
for (i in colnames(contrasts)) {
  test <- QLtests[[i]]
  results <- topTags(test, n = nrow(test$table), adjust.method = "fdr")
  results <- tibble(data.frame(results$table)) %>% 
    mutate(
      Expression = case_when(logFC >= log(2) & FDR <= 0.05 ~ "Up-regulated",
                             logFC <= -log(2) & FDR <= 0.05 ~ "Down-regulated",
                             TRUE ~ "Unchanged")
    ) %>% 
    mutate(
      Significance = case_when(
        abs(logFC) >= log(2) & FDR <= 0.05 & FDR > 0.01 ~ "FDR 0.05", 
        abs(logFC) >= log(2) & FDR <= 0.01 & FDR > 0.001 ~ "FDR 0.01",
        abs(logFC) >= log(2) & FDR <= 0.001 ~ "FDR 0.001", 
        TRUE ~ "Unchanged")
    ) %>% 
    dplyr::rename(gene = genes) %>% 
    left_join(waterhemp_annot, by ="gene")
  DE_tables[[i]] <- results
}


# get all DE genes
DE_genes_all <- lapply(DE_tables, function(df){
  df %>% 
    filter(Expression != "Unchanged") %>%
    distinct(gene,.keep_all = T)
})


DE_genes_all$MaleFMvsFemaleFM |> View()


names(DE_genes_all)
# number of DEs
lapply(DE_tables, function(df){
  df %>%
    distinct(gene,.keep_all = T) |> 
    group_by(Expression) %>% 
    summarise(n = n())
})

# get number of DE genes
n_DE <- data.frame(sapply(DE_genes_all, function(df){
  nrow(df %>%
         distinct(gene,.keep_all = T))
}))

names(n_DE) <- "DE_genes"

n_DE <- cbind(n_DE, comparison = c(rep("Across gender",4), rep("Within gender", 6)))
n_DE <- n_DE %>% rownames_to_column("Contrast")

ggplot(n_DE[-4,], aes(y = DE_genes, x = Contrast, fill = comparison)) + geom_col(position = "dodge") +
  facet_wrap(comparison ~ ., scales = "free") +
  labs(x = "Comparison", y = "DE genes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label=DE_genes), position=position_dodge(width=0.9), vjust=-0.25) 

ggsave("waterhemp/results/n_DE_all_comparisons.pdf", dpi = 300, width = 30, height = 20, units = "cm")  


### Venn diagram ------------

if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)
# Across gender comparisons
Venn_genes_gender <- list(
  `Mature Flower` = DE_genes_all$MaleFlowerVsFemaleFlower$gene,
  SAM =  DE_genes_all$MaleSAMvsFemaleSAM$gene,
  `Floral meristem` =  DE_genes_all$MaleFMvsFemaleFM$gene
)

gender_common <- intersect(intersect(Venn_genes_gender$`Mature Flower`,Venn_genes_gender$`Floral meristem`),Venn_genes_gender$SAM)



waterhemp_venn_across_gender <- ggvenn(
  Venn_genes_gender, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
) + labs(title = "Across genders comparison") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("waterhemp/results/Across_gender_venn.pdf", width = 20, height = 20, units = "cm", dpi = 300)


# within gender comparisons


# What genes are unique for each gender across each tissue

# SAM vs FM

Venn_genes_SAMvsFM <- list(
  `Female` = DE_genes_all$Female_SAMvsFM$gene,
  `Male` = DE_genes_all$Male_SAMvsFM$gene
) 


Venn_SAMvsFM <- ggvenn(
  Venn_genes_SAMvsFM, 
  fill_color = c("#F9502C", "#2C9FF9"),
  stroke_size = 0.5, set_name_size = 5
) + labs(title = "SAM vs Floral meristem") +
  theme(plot.title = element_text(hjust = 0.5))



SAMvsFM_female_unique <- setdiff(Venn_genes_SAMvsFM$Female, Venn_genes_SAMvsFM$Male)
SAMvsFM_male_unique <- setdiff(Venn_genes_SAMvsFM$Male, Venn_genes_SAMvsFM$Female)

SAMvsFM_female_unique <- DE_genes_all$Female_SAMvsFM %>% 
  filter(gene %in% SAMvsFM_female_unique)

SAMvsFM_male_unique <- DE_genes_all$Male_SAMvsFM %>% 
  filter(gene %in% SAMvsFM_male_unique)


# SAM vs Mature flower

Venn_genes_SAMvsFlower <- list(
  `Female` = DE_genes_all$Female_SAMvsFlower$gene,
  `Male` = DE_genes_all$Male_SAMvsFlower$gene
) 

Venn_SAMvsFlower <- ggvenn(
  Venn_genes_SAMvsFlower, 
  fill_color = c("#F9502C", "#2C9FF9"),
  stroke_size = 0.5, set_name_size = 5
) + labs(title = "SAM vs Mature Flower") +
  theme(plot.title = element_text(hjust = 0.5))

SAMvsFlower_female_unique <- setdiff(Venn_genes_SAMvsFlower$Female, Venn_genes_SAMvsFlower$Male)
SAMvsFlower_male_unique <- setdiff(Venn_genes_SAMvsFlower$Male, Venn_genes_SAMvsFlower$Female)


SAMvsFlower_female_unique <- DE_genes_all$Female_SAMvsFlower %>% 
  filter(gene %in% SAMvsFlower_female_unique)

SAMvsFlower_male_unique <- DE_genes_all$Male_SAMvsFlower %>% 
  filter(gene %in% SAMvsFlower_male_unique)


# FM vs Mature flower

Venn_genes_FMvsFlower <- list(
  `Female` = DE_genes_all$Female_FMvsFlower$gene,
  `Male` = DE_genes_all$Male_FMvsFlower$gene
) 

Venn_FMvsFlower <- ggvenn(
  Venn_genes_FMvsFlower, 
  fill_color = c("#F9502C", "#2C9FF9"),
  stroke_size = 0.5, set_name_size = 5
) + labs(title = "Floral meristem vs Mature Flower") +
  theme(plot.title = element_text(hjust = 0.5))

FMvsFlower_female_unique <- setdiff(Venn_genes_FMvsFlower$Female, Venn_genes_FMvsFlower$Male)
FMvsFlower_male_unique <- setdiff(Venn_genes_FMvsFlower$Male, Venn_genes_FMvsFlower$Female)
length(FMvsFlower_male_unique)

FMvsFlower_female_unique <- DE_genes_all$Female_FMvsFlower %>% 
  filter(gene %in% FMvsFlower_female_unique) %>% 
  distinct(gene, .keep_all = T)

FMvsFlower_male_unique <- DE_genes_all$Male_FMvsFlower %>% 
  filter(gene %in% FMvsFlower_male_unique) %>% 
  distinct(gene, .keep_all = T)


# get venn diagram joint 
library(patchwork)
all_Ven_waterhemp <- Venn_SAMvsFM + Venn_SAMvsFlower + Venn_FMvsFlower

ggsave(plot = all_Ven_waterhemp, filename = "waterhemp/results/waterhemp_venn_all.pdf", width = 35, height = 20, units = "cm", dpi = 300)

# SAVE LIST OF DE ---------------

setwd("waterhemp/results/DE_genes_all_comparisons")
sapply(names(DE_tables), 
       function (x) write.table(DE_tables[[x]], 
                                file=paste0(x, "_ALLgenes_EDGER.txt"),
                                row.names = F, quote = F, sep = "\t"))

setwd("../../..")


files <- dir("waterhemp/results/DE_genes_all_comparisons", pattern = "*.txt",full.names = T)

df_names <- str_remove(dir("waterhemp/results/DE_genes_all_comparisons", pattern = "*.txt"),
                       pattern = "_ALLgenes_EDGER.txt")



DE_all_comparisons_WH <- DE_tables |> 
  bind_rows(.id = "Comparison")

DE_all_comparisons_WH |> filter(uniprot == "SUP_ARATH") |> View()
DE_all_comparisons_WH <- DE_all_comparisons_WH |> 
  dplyr::select(gene,uniprot, blastP_annot, Pfam_annot, Comparison,
                logFC,Expression, Significance, scaffold, start, end, Identity) |> 
  mutate(Significance = str_remove(Significance, "FDR ")) |> 
  mutate(group = paste(round(logFC,2),Expression,Significance, sep = "_")) |> 
  dplyr::select(-logFC,-Expression, - Significance) |>
  pivot_wider(names_from = Comparison, values_from = group)


write_csv(DE_all_comparisons_WH, "waterhemp/results/DE_genes_all_comparisons/DE_all_comparisons_WH.csv")


# MSY ----------

WH_MSY <- DE_all_comparisons_WH |> 
  filter(scaffold %in% c('tig00100752', 'tig00000455', 'tig00000336',
                         'tig00000340', 'tig00000542', 'tig00000546', 
                         'tig00000298', 'tig00001274', 'tig00000740', 
                         'tig00000708', 'tig00003161', 'tig00000595', 
                         'tig00000080', 'tig00004323'))




write_csv(WH_MSY, "waterhemp/results/DE_genes_all_comparisons/WH_MSY.csv")

