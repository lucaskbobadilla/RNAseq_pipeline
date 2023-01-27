# GOTERM ANALYSIS WATERHEMP FLOWER
library(tidyverse)
library(topGO)


# load data
#DE_genes <- read_csv("waterhemp/results/Male_vs_Female/sign_Male_vs_Female.csv")
DE_genes <- SAMvsFlower_male_unique
# Load GOterm annot

# GOterms --------------

GOterms_annot <- read_csv("../trinity_assembly/trinity_waterhemp/annotation/GOterm_blastP.csv")
GOterms_annot_pfam <- read_csv("../trinity_assembly/trinity_waterhemp/annotation/GOterm_Pfam.csv")
GOterms_annot_blastX <- read_csv("../trinity_assembly/trinity_waterhemp/annotation/GOterm_blastX.csv")

# organize
GOterms_annot_pfam <- GOterms_annot_pfam %>% 
  group_by(gene,ontology) %>% 
  #distinct(gene, .keep_all = T) %>% 
  dplyr::select(gene, go, ontology, GOname= name) %>% 
  mutate(join_fac = paste(gene,go, sep = "-")) %>% 
  distinct(join_fac, .keep_all = T)
 
 
GOterms_annot <- GOterms_annot %>% 
  group_by(gene,ontology) %>% 
  #distinct(gene, .keep_all = T) %>% 
  dplyr::select(gene, go, ontology, GOname= name) %>% 
  mutate(join_fac = paste(gene,go, sep = "-")) %>% 
  distinct(join_fac, .keep_all = T) 

GOterms_annot_blastX <- GOterms_annot_blastX %>% 
  group_by(gene,ontology) %>% 
  #distinct(gene, .keep_all = T) %>% 
  dplyr::select(gene, go, ontology, GOname= name) %>% 
  mutate(join_fac = paste(gene,go, sep = "-")) %>% 
  distinct(join_fac, .keep_all = T) 

GOterms_annot <- GOterms_annot %>%
  full_join(GOterms_annot_blastX) %>% 
  full_join(GOterms_annot_pfam) %>% 
  distinct(join_fac, .keep_all = T)

GOterms_annot <- GOterms_annot %>%
  filter(gene %in% DE_genes$genes) %>%
  ungroup()

 
# # prepare GOterm database
# 
# geneID2GO <- unlist(as.list(
#   data.frame(GOterms_annot %>% 
#                arrange(gene) %>% 
#                dplyr::select(gene, go) %>% 
#                pivot_wider(names_from = gene, values_from = go))
# ),recursive=FALSE)
# 
# # check the list
# str(head(geneID2GO))


# Importing Trinotate list names --------------------
Trinotate_go <- read_delim("~/lucas/gender_rnaseq/trinity_assembly/trinity_waterhemp/annotation/Trinotate.xls.gene_ontology", 
                           delim = "\t", escape_double = FALSE, 
                           col_names = FALSE, trim_ws = TRUE)

Trinotate_go <- Trinotate_go %>% mutate(X1 = gsub("_i.*", "", X1)) %>% 
  mutate(n = nchar(X2)) %>% 
  group_by(X1) %>% 
  arrange(n, .by_group = T) %>% 
  distinct(X1, .keep_all = T)



head(Trinotate_go)
Trinotate_go <- data.frame(Trinotate_go)
rownames(Trinotate_go) <- Trinotate_go$X1


# convert to list
geneID2GO <- list()
for (i in Trinotate_go$X2) {
  GO <- unlist(strsplit(i,","))
  geneID2GO[[length(geneID2GO) + 1]] <- GO
}

# change list names
names(geneID2GO) <- Trinotate_go$X1

# check the list
str(head(geneID2GO))


# get gene names
geneNames <- names(geneID2GO)

# get genes of interest --------------
gene_of_interest <- (DE_genes %>% filter(Expression == "Up-regulated"))$gene
gene_of_interest <- (SAMvsFlower_male_unique$gene)
# get final gene list
geneList <- factor(as.integer(geneNames %in% gene_of_interest))
names(geneList) <- geneNames
str(geneList)


# start GOterm analysis - Biological function -----------------

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO,
              description = "Enrichment analysis for waterhemp - male vs Female")

# define test using the classic algorithm with ks 
classic_ks_result <- runTest(GOdata, algorithm='classic', statistic='fisher')
# define test using the weight01 algorithm with ks
weight_ks_result <- runTest(GOdata, algorithm='weight01', statistic='fisher') 

# generate a table of results:
# We can use the GenTable function to generate a summary table with the results 
# from tests applied to the topGOdata object.

allGO <- usedGO(GOdata)
all_res <- GenTable(GOdata, weightks= weight_ks_result,
                    classic_ks_result = classic_ks_result, orderBy='weightks', topNodes=length(allGO))

glimpse(all_res)

all_res <- all_res %>% mutate(weightks = as.numeric(str_remove_all(weightks, pattern = "< ")),
                              classic_ks_result = as.numeric(str_remove_all(classic_ks_result, pattern = "< ")))
#performing BH correction on our p values
p.adj_weight <- p.adjust(all_res$weightks,method="fdr")
p.adj_classic <- p.adjust(all_res$classic_ks_result, method="fdr")

all_res <- tibble(all_res) %>% 
  dplyr::rename(weightks_pval = weightks, 
                classic_rank = `Rank in classic_ks_result`,
                classic_ks_pval =classic_ks_result) %>% 
   mutate(p.adj_weight = p.adj_weight,
          p.adj_classic = p.adj_classic)

# Conditional enrichment
cond_enrich <- all_res %>% 
  filter(weightks_pval <= 0.05)

# Over-representation enrichment

overRep_enrich <- all_res %>% 
  filter(classic_ks_pval <= 0.05)

View(all_res %>% 
       filter(p.adj_weight <= 0.05 | p.adj_classic <= 0.05))

pdf(file='waterhemp/results/flower_comp/topGOPlot_classic_top5.pdf', height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(GOdata, score(classic_ks_result), firstSigNodes = 30, useInfo = 'all')
dev.off()

write_csv(all_res, "waterhemp/results/Male_vs_Female/GOenrichment_analysis.csv")

glimpse(all_res)

# top GO terms
all_res %>%
  mutate(`Conditional enrichment` = as.numeric(weightks_pval),
         `Over-representation enrichment` = as.numeric(classic_ks_pval)) %>% 
  filter(`Conditional enrichment` <= 0.05 & `Over-representation enrichment` <= 0.05) %>%
  head(n = 30) %>% 
  ggplot(aes(y = Term, x = Annotated, colour = `Conditional enrichment`, size = `Over-representation enrichment`)) +
  geom_point() +
  theme_light() +
  theme(legend.position = "right") +
  labs(title = "Enriched GO-terms on both tests - Male vs Female Waterhemp")

sign_GOterms <- overRep_enrich %>% 
  full_join(cond_enrich)

both_GOterms <- all_res %>%
  mutate(`Conditional enrichment` = as.numeric(weightks_pval),
         `Over-representation enrichment` = as.numeric(classic_ks_pval)) %>% 
  filter(`Conditional enrichment` <= 0.05 & `Over-representation enrichment` <= 0.05)

write_csv(both_GOterms,
          "waterhemp/results/tissue_per_gender_comp/GOterms_SAMvsFlower_male_unique.csv")

# get genes in significant terms -------------------------

mygenes = genesInTerm(GOdata, both_GOterms$GO.ID)
glimpse(tibble(plyr::ldply(mygenes, rbind)))

View(DE_genes_enriched %>% 
       distinct(gene, .keep_all = T) %>% 
       group_by(scaffold) %>% 
       summarise(n = n()) %>% 
       arrange(desc(n)))


DE_genes_enriched <- tibble(plyr::ldply(mygenes, rbind)) %>% 
  pivot_longer(cols = `1`:`4480`,names_to = "rem", values_to = "gene") %>% 
  filter(gene %in% gene_of_interest) %>% 
  dplyr::select(-rem) %>% 
  left_join(all_res %>% dplyr::rename(.id = GO.ID) %>% dplyr::select(.id:Term)) %>% 
  left_join(waterhemp_annot) %>% 
  arrange(gene) 




DE_genes_enriched %>% 
  distinct(gene, .keep_all = T) %>% 
  View()

DE_genes_enriched_set <- DE_genes_enriched %>% 
  filter(.id %in% c(
    'GO:0009908'
  )) %>% 
  distinct(gene, .keep_all = T) %>%
  dplyr::select(.id:blastP_annot, scaffold, start, end)

DE_genes_enriched_set <- DE_genes_enriched_set %>% 
  left_join(DE_genes %>% dplyr::select(gene,logFC, FDR, Expression))


write_csv(DE_genes_enriched, "waterhemp/results/tissue_per_gender_comp/Genes_SAMvsFlower_male_unique.csv")
# OLD -------------------
# select meaningful terms


GOterms <- c("GO:0010197", "GO:0051301", "GO:0009567", "GO:0051865", "GO:1902182",
           "GO:0010227", "GO:0009909", "GO:0010074", "GO:0009559", "GO:0090700",
           "GO:0009566", "GO:0045596", "GO:2000242","GO:0009690")

all_res %>%
  mutate(`Conditional enrichment` = as.numeric(weightks_pval),
         `Over-representation enrichment` = as.numeric(classic_ks_pval)) %>% 
  filter(GO.ID %in% GOterms) %>% 
  ggplot(aes(y = Term, x = Annotated, colour = `Conditional enrichment`, size = `Over-representation enrichment`)) +
  geom_point() +
  theme_light() +
  theme(legend.position = "right") +
  labs(title = "Enriched GO-terms related to reproduction- Male vs Female Waterhemp")

sign_terms_annot <- GOterms_annot %>% 
  filter(go %in% sign_GOterms$GO.ID) %>% 
  filter(gene %in% DE_genes$genes) %>% 
  left_join(waterhemp_annot_only_plants %>% dplyr::select(gene,uniprot, blastP_annot), by = "gene") %>% 
  left_join(DE_genes %>% dplyr::rename(gene = genes)) %>% 
  filter(!str_detect(blastP_annot, "Retrovirus"))

write_csv(sign_terms_annot, "waterhemp/results/Male_vs_Female/sign_GOterms_annot.csv")

reproduction_terms <- GOterms_annot %>% 
  filter(go %in% GOterms) %>% 
  filter(gene %in% DE_genes$genes) %>% 
    left_join(waterhemp_annot_only_plants %>% dplyr::select(gene,uniprot, blastP_annot), by = "gene") %>% 
  left_join(DE_genes %>% dplyr::rename(gene = genes))

write_csv(reproduction_terms, "waterhemp/results/Male_vs_Female/genes_annot_Goterms_reproduction.csv")
