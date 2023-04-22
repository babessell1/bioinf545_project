library(tidyverse)
library(biomaRt)
library(gprofiler2)
library(ggplot2)
library(dplyr)
setwd('/Users/xieshuyi/Desktop/Callipyge/bioinf545_project/DGE/DESeq_table')
CC_NC_no_dup <- read.table("CC_vs_NC_no_dup_cow", header = TRUE)
DEG_CC_NC_no_dup <-subset(CC_NC_no_dup, pvalue < 0.05 & abs(log2FoldChange) > 1)
ensembl<-  useMart("ensembl", dataset="btaurus_gene_ensembl")
gene_set_CC_NC_no_dup <- row.names(DEG_CC_NC_no_dup)
conversion_table <- getBM(attributes=c("uniprot_gn_symbol", "ensembl_gene_id"), filters = "uniprot_gn_symbol", values = gene_set_CC_NC_no_dup, mart= ensembl)
uniprot_gn_symbol<-row.names(DEG_CC_NC_no_dup)
DEG_CC_NC_no_dup<-cbind(uniprot_gn_symbol,DEG_CC_NC_no_dup)
Merge_CC_NC_no_dup<- merge(DEG_CC_NC_no_dup,conversion_table,by="uniprot_gn_symbol")
CC_NC_no_dup_down <- filter(Merge_CC_NC_no_dup,log2FoldChange < -1)
CC_NC_no_dup_up <- filter(Merge_CC_NC_no_dup,log2FoldChange > 1)
CC_NC_no_dup_down  <- CC_NC_no_dup_down [order(CC_NC_no_dup_down$pvalue),]
CC_NC_no_dup_up  <- CC_NC_no_dup_up[order(CC_NC_no_dup_up$pvalue),]
gp_CC_NC_no_dup<- gost(query = list("upregulated" = CC_NC_no_dup_up$ensembl_gene_id,
                                    "downregulated" = CC_NC_no_dup_down$ensembl_gene_id), 
                       organism = "btaurus", ordered_query = TRUE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                       measure_underrepresentation = FALSE, evcodes = TRUE, 
                       user_threshold = 0.05, correction_method = "g_SCS", 
                       domain_scope = "annotated", custom_bg = NULL, 
                       numeric_ns = "", sources = "GO", as_short_link = FALSE)
p<-gostplot(gp_CC_NC_no_dup, capped = TRUE, interactive = FALSE)
pp <- publish_gostplot(p, highlight_terms = gp_CC_NC_no_dup$result[c(4,16,60,64),],
                       width = 0.00000002, height = 10, filename = NULL)
gp_mod = gp_CC_NC_no_dup$result[,c("query", "source", "term_id",
                            "term_name", "p_value", "query_size", 
                            "intersection_size", "term_size", 
                            "effective_domain_size", "intersection")]
gp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)
gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                  "query_size", "Count", "term_size", "effective_domain_size", 
                  "geneID", "GeneRatio", "BgRatio")
gp_mod$geneID = gsub(",", "/", gp_mod$geneID)
row.names(gp_mod) = gp_mod$ID
gp_mod_enrich  = new("enrichResult", result = gp_mod)
barplot(gp_mod_enrich, showCategory = 10, font.size = 8) + 
  ggplot2::facet_grid(~Cluster) +
  ggplot2::ylab("Intersection size")
gem <- gp_CC_NC_no_dup$result[,c("query", "term_id", "term_name", "p_value", "intersection")]
colnames(gem) <- c("query", "term_id", "Description", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem %>% group_by(query) %>%
  group_walk(~
               write.table(data.frame(.x[,c("term_id", "Description", "p.Val", "FDR", "Genes")]), 
                           file = paste0("/Users/xieshuyi/Desktop/Callipyge/bioinf545_project/pathway/gProfiler_", unique(.y$query), "_gem.txt"),
                           sep = "\t", quote = F, row.names = F))