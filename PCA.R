### PCA waterhemp ############


library(tidyverse)
library(pcaExplorer)
library(DESeq2)


# load metadata ========================

sample_info <- read_csv("waterhemp/data/waterhemp_metadata.csv")

sample_info <- data.frame(sample_info)

rownames(sample_info) <- sample_info$sample

sample_info <- sample_info[-1]

sample_info <- sample_info %>% 
  mutate(groups = factor(paste(gender,tissue,sep = ".")))


col_data <- sample_info %>% 
  mutate(gender = factor(gender),
         tissue = factor(tissue))
# load counts

counts <- read.table("../DE_analysis/waterhemp/data/raw_counts.txt", sep = "\t")
counts[1:5,1:5]
dim(counts)
t(gene_count[1:5,1:5])
# load annotation

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

waterhemp_annot_only_plants <-  waterhemp_annot %>% 
  filter(str_detect(lineage,"Viridiplantae"))

annot <- waterhemp_annot%>% 
  dplyr::select(gene, gene_name = uniprot) 

annot <- data.frame(annot)
rownames(annot) <- annot$gene
annot <- annot[-1]
# %>% 
#   separate(gene, into = c("r1","r2"),sep = "_GG_",remove = F) %>% 
#   mutate(gene_name = paste(r2, uniprot, sep = "_")) %>% 
#   select(-r1,-r2, -uniprot)

# filter only annotated genes

count_plants_only <- data.frame(counts) %>% filter(rownames(counts) %in% waterhemp_annot_only_plants$gene)
count_annot_only <-  data.frame(counts) %>% filter(rownames(counts) %in% waterhemp_annot$gene)
write.table(count_plants_only, "../../data/annotated_counts.txt",
            sep = "\t", quote = F,
            col.names = T, row.names = T)
#create DESeq2 object

dds <- DESeqDataSetFromMatrix(countData = round(count_annot_only),
                              colData = col_data,
                              design= ~ gender + tissue)




ddt <- vst(dds)

# run pca explorer

pcaExplorer(dds = dds, annotation = annot)


# create plot for gene expression

WH_genes <- read_csv("waterhemp/results/WH_expression_genes.csv")





WH_genes <- WH_genes %>% 
  dplyr::rename(sample = `...1`) %>% 
  mutate(tissue = case_when(tissue == "Flower" ~ "Mature flower",
                            tissue == "Flower_meristem" ~ "Floral meristem",
                            tissue == "SA_meristem" ~ "SAM"),
         tissue = factor(tissue, levels = c("SAM", "Floral meristem", "Mature flower"))) 

WH_genes <- WH_genes |> 
  full_join(ORR24 |> mutate(Gene = "ORR24")) |> 
  full_join(SUP |> mutate(Gene = "SUP"))
factor(WH_genes$Gene)
WH_genes$tissue <- factor(WH_genes$tissue, levels = c("SAM", "Floral meristem", "Mature flower"))
WH_genes$Gene <- factor(WH_genes$Gene, levels = c("MADS18", "ORR24", "SUP", "LOB31", "PISTILLATA PMADS2", "MADS-CMB2", "MADS2", "MYB35"))

WH_plot_genes <- WH_genes |> 
  ggplot(aes(y = count, x = tissue, fill = gender)) +
  geom_boxplot() +
  facet_wrap(Gene ~ ., scales = "free", ncol = 4) +
  labs(x = "tissue type", y = "Log-normalized expression count", title = expression(italic("Amaranthus tuberculatus"))) + 
  theme_light() +
  theme(legend.position = "bottom") +
  theme(strip.text.x = element_text(face = "bold", colour = "black"))

WH_plot_genes +
  theme(text = element_text(size = 28),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("waterhemp/results/waterhemp_key_genes_presentation.png", dpi = 300, height = 12, width = 24, units = "in")

# MADSplot

MADS18 <- read_csv("waterhemp/results/MADS18_counts.csv")

MADS18 <- MADS18 %>% 
  rename(sample = `...1`) %>% 
  mutate(tissue = case_when(tissue == "Flower" ~ "Mature flower",
                            tissue == "Flower_meristem" ~ "Floral meristem",
                            tissue == "SA_meristem" ~ "SAM"))

MADS18 %>% 
  group_by(tissue,gender) %>% 
  summarise(avg = mean(count),
            sd = sd(count)) %>% 
  ggplot(aes(y = avg, x = tissue, fill = gender)) +
  geom_col(position = position_dodge(.9), width = 0.5) +
  theme_classic() +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2,
                position=position_dodge(.9))  +
  theme(legend.position = "bottom") +
  labs(x = "tissue type", y = "Log-normalized expression count", title = "MADS-box 18 - Waterhemp")

ggsave("waterhemp/results/MADS18.pdf", width = 20, height = 15, units = "cm", dpi = 300)


#ORR24

ORR24 <- read_csv("waterhemp/results/ORR24_waterhemp.csv")

ORR24 <- ORR24 %>% 
  dplyr::rename(sample = `...1`) %>% 
  mutate(tissue = case_when(tissue == "Flower" ~ "Mature flower",
                            tissue == "Flower_meristem" ~ "Floral meristem",
                            tissue == "SA_meristem" ~ "SAM"),
         tissue = factor(tissue, levels = c("SAM", "Floral meristem", "Mature flower"))) 

ORR24 %>% 
  group_by(tissue,gender) %>% 
  summarise(avg = mean(count),
            sd = sd(count)) %>% 
  ggplot(aes(y = avg, x = tissue, fill = gender)) +
  geom_col(position = position_dodge(.9), width = 0.5) +
  theme_classic() +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2,
                position=position_dodge(.9))  +
  theme(legend.position = "bottom") +
  labs(x = "tissue type", y = "Log-normalized expression count", title = "ORR24 - Waterhemp")

ggsave("waterhemp/results/ORR24.pdf", width = 20, height = 15, units = "cm", dpi = 300)

#SUP

SUP <- read_csv("waterhemp/results/SUP_expression_WH.csv")

SUP <- SUP %>% 
  dplyr::rename(sample = `...1`) %>% 
  mutate(tissue = case_when(tissue == "Flower" ~ "Mature flower",
                            tissue == "Flower_meristem" ~ "Floral meristem",
                            tissue == "SA_meristem" ~ "SAM"))

SUP %>% 
  group_by(tissue,gender) %>% 
  summarise(avg = mean(count),
            sd = sd(count)) %>% 
  ggplot(aes(y = avg, x = tissue, fill = gender)) +
  geom_col(position = position_dodge(.9), width = 0.5) +
  theme_classic() +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2,
                position=position_dodge(.9))  +
  theme(legend.position = "bottom") +
  labs(x = "tissue type", y = "Log-normalized expression count", title = "SUP - Waterhemp")

ggsave("waterhemp/results/SUP.png", width = 20, height = 15, units = "cm", dpi = 300)



# plots

pc12 <- pcaplot(ddt,intgroup = c("gender","tissue"), pcX = 1, pcY = 2, 
                title = "PC1 vs PC2", text_labels = FALSE) +
  theme_classic() + 
  scale_colour_discrete(name = "Gender/tissue: ", 
                        labels = c("Female/Mature Flower","Female/Floral meristem", "Female/SAM",
                                   "Male/Mature Flower", "Male/Floral meristem", "Male/SAM")) +
  theme(legend.position = "bottom") 

pc13 <- pcaplot(ddt,intgroup = c("gender","tissue"), pcX = 1, pcY = 3, 
                title = "PC1 vs PC3", text_labels = FALSE) +
  theme_classic() + 
  scale_colour_discrete(name = "Gender/tissue: ", 
                        labels = c("Female/Mature Flower","Female/Floral meristem", "Female/SAM",
                                   "Male/Mature Flower", "Male/Floral meristem", "Male/SAM")) +
  theme(legend.position = "bottom") 

pc23 <- pcaplot(ddt,intgroup = c("gender","tissue"), pcX = 2, pcY = 3, 
                title = "PC2 vs PC3", text_labels = FALSE) +
  theme_classic() + 
  scale_colour_discrete(name = "Gender/tissue: ", 
                        labels = c("Female/Mature Flower","Female/Floral meristem", "Female/SAM",
                                   "Male/Mature Flower", "Male/Floral meristem", "Male/SAM")) +
  theme(legend.position = "bottom")

library(patchwork)

PCs_plots <- pc12 + pc13 + pc23 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

#scree
pcaobj_ddt <- prcomp(t(assay(ddt)))
scree <- pcascree(pcaobj_ddt,type="pev",
                  title="Proportion of explained proportion of variance") + scale_y_continuous(labels = scales::percent)

final_waterhemp_PC_plots <- PCs_plots / scree
ggsave("waterhemp/results/PCA/waterhemp_PCA_final_plot.pdf", width = 15, height = 12, units = "in", dpi = 300)

# correlate PCs
res_pcairway <- correlatePCs(pcaobj_ddt,colData(ddt))
res_pcairway <- res_pcairway[,1:3]

pdf("waterhemp/results/PCA/PC_corretions.pdf", width = 12, height = 8)
par(mfrow=c(2,2))
corPC1 <- plotPCcorrs(res_pcairway)
corPC2 <- plotPCcorrs(res_pcairway,pc = 2)
corPC3 <- plotPCcorrs(res_pcairway,pc = 3)
corPC4 <- plotPCcorrs(res_pcairway,pc = 4)
dev.off()

# extract the table of the genes with high loadings
load_PC1 <- hi_loadings(pcaobj_ddt,topN = 30,exprTable=assay(ddt),whichpc = 1)
load_PC1_annot <- annot_blastp %>% 
  filter(gene %in% rownames(load_PC1)) %>% 
  select(gene, uniprot, blastP_annot) %>% 
  left_join(data.frame(load_PC1) %>% rownames_to_column("gene")) 

load_PC2 <- hi_loadings(pcaobj_ddt,topN = 30,exprTable=assay(ddt),whichpc = 2)
load_PC2_annot <- annot_blastp %>% 
  filter(gene %in% rownames(load_PC2)) %>% 
  select(gene, uniprot, blastP_annot) %>% 
  left_join(data.frame(load_PC2) %>% rownames_to_column("gene")) 

load_PC3 <- hi_loadings(pcaobj_ddt,topN = 30,exprTable=assay(ddt),whichpc = 3)
load_PC3_annot <- annot_blastp %>% 
  filter(gene %in% rownames(load_PC3)) %>% 
  select(gene, uniprot, blastP_annot) %>% 
  left_join(data.frame(load_PC3) %>% rownames_to_column("gene")) 

load_PC4 <- hi_loadings(pcaobj_ddt,topN = 30,exprTable=assay(ddt),whichpc = 4)
load_PC4_annot <- annot_blastp %>% 
  filter(gene %in% rownames(load_PC4)) %>% 
  select(gene, uniprot, blastP_annot) %>% 
  left_join(data.frame(load_PC4) %>% rownames_to_column("gene")) 

pdf("waterhemp/results/PCA/waterhemp_PC_loadings.pdf", width = 12, height = 8)
par(mar = c(8, 4, 4, 4), mfrow=c(2,2))
hi_loadings(pcaobj_ddt,topN = 30,annotation = annot,whichpc = 1)
hi_loadings(pcaobj_ddt,topN = 30,annotation = annot,whichpc = 2)
hi_loadings(pcaobj_ddt,topN = 30,annotation = annot,whichpc = 3)
hi_loadings(pcaobj_ddt,topN = 30,annotation = annot,whichpc = 4)
dev.off()


