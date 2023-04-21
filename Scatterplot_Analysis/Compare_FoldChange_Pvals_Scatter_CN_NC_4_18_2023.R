#' The purpose of this file is to merge the datasets provided by DESeq2 and
#' both generate scatterplots of log-2 fold change for each gene across species
#' alignment and across deduplication status, and scatterplots of Wald test
#' statistics in the same conditions. Additionally, we identify genes with the
#' five highest Cook's D in our various datasets.
#'
#' This file in particular looks at the CN vs. NC comparison.
#'
#' Author: Hasan Abu-Amara
#' Date: 4-18-2023

library(ggplot2)
library(qqman)
library(data.table)
library(tidyverse)
library(dplyr)
library(broom)
library(ggrepel)

workspace_chromebook = FALSE
results_dir = ifelse(workspace_chromebook, "/mnt/chromeos/removable/HASAN LEXAR/BIOSTAT 646/Project/DESeq2 Katarina 4-11-2023/", "H:/BIOSTAT 646/Project/DESeq2 Katarina 4-11-2023/")
cc_nc_nodup_sheep = fread(paste0(results_dir, "CC_vs_NC_no_dup_sheep"))
cc_nc_nodup_cow = fread(paste0(results_dir, "CC_vs_NC_no_dup_cow"))
cc_nc_dup_sheep = fread(paste0(results_dir, "CC_vs_NC_with_dup_sheep"))
cc_nc_dup_cow = fread(paste0(results_dir, "CC_vs_NC_with_dup_cow"))
nn_nc_nodup_sheep = fread(paste0(results_dir, "NN_vs_NC_no_dup_sheep"))
nn_nc_nodup_cow = fread(paste0(results_dir, "NN_vs_NC_no_dup_cow"))
nn_nc_dup_sheep = fread(paste0(results_dir, "NN_vs_NC_with_dup_sheep"))
nn_nc_dup_cow = fread(paste0(results_dir, "NN_vs_NC_with_dup_cow"))
cn_nc_nodup_sheep = fread(paste0(results_dir, "CN_vs_NC_no_dup_sheep"))
cn_nc_nodup_cow = fread(paste0(results_dir, "CN_vs_NC_no_dup_cow"))
cn_nc_dup_sheep = fread(paste0(results_dir, "CN_vs_NC_with_dup_sheep"))
cn_nc_dup_cow = fread(paste0(results_dir, "CN_vs_NC_with_dup_cow"))

# Okay. That is all of the datasets loaded.
# Let's make some scatterplots.

cn_nc_nodup_cow_sheep = merge(cn_nc_nodup_cow, cn_nc_nodup_sheep, by = "V1", suffixes = c(".cow", ".sheep"))
ggplot(cn_nc_nodup_cow_sheep, aes(x = stat.cow, y = stat.sheep)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic in Cow Alignment", y = "Test Statistic in Sheep Alignment", title = "CN Vs. NC Alignment Test Statistic Comparison Across Species - No Duplicates") +
  theme_classic()
ggplot(cn_nc_nodup_cow_sheep, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "CN Vs. NC Alignment log2-Fold-Change Comparison Across Species - No Duplicates") +
  theme_classic()

cor(cn_nc_nodup_cow_sheep$stat.cow, cn_nc_nodup_cow_sheep$stat.sheep)
cor(cn_nc_nodup_cow_sheep$log2FoldChange.cow, cn_nc_nodup_cow_sheep$log2FoldChange.sheep)

cn_nc_dup_cow_sheep = merge(cn_nc_dup_cow, cn_nc_dup_sheep, by = "V1", suffixes = c(".cow", ".sheep"))
ggplot(cn_nc_dup_cow_sheep, aes(x = stat.cow, y = stat.sheep)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic in Cow Alignment", y = "Test Statistic in Sheep Alignment", title = "CN Vs. NC Alignment Test Statistic Comparison Across Species - With Duplicates") +
  theme_classic()
ggplot(cn_nc_dup_cow_sheep, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "CN Vs. NC Alignment log2-Fold-Change Comparison Across Species - With Duplicates") +
  theme_classic()

cor(cn_nc_dup_cow_sheep$stat.cow, cn_nc_dup_cow_sheep$stat.sheep)
cor(cn_nc_dup_cow_sheep$log2FoldChange.cow, cn_nc_dup_cow_sheep$log2FoldChange.sheep)

cn_nc_dup_no_dup_cow = merge(cn_nc_dup_cow, cn_nc_nodup_cow, by = "V1", suffixes = c(".dup", ".nodup"))
cn_nc_dup_no_dup_sheep = merge(cn_nc_dup_sheep, cn_nc_nodup_sheep, by = "V1", suffixes = c(".dup", ".nodup"))
cor(cn_nc_dup_no_dup_cow$log2FoldChange.dup, cn_nc_dup_no_dup_cow$log2FoldChange.nodup)
cor(cn_nc_dup_no_dup_sheep$log2FoldChange.dup, cn_nc_dup_no_dup_sheep$log2FoldChange.nodup)

ggplot(cn_nc_dup_no_dup_cow, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "CN Vs. NC Duplicates log2-Fold-Change Comparison - Cows") +
  theme_classic()
ggplot(cn_nc_dup_no_dup_cow, aes(x = stat.dup, y = stat.nodup)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Stats with Duplicates", y = "Test Stats without Duplicates", title = "CN Vs. NC Duplicates Test Statistics Comparison - Cows") +
  theme_classic()
ggplot(cn_nc_dup_no_dup_sheep, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "CN Vs. NC Duplicates log2-Fold-Change Comparison - Sheep") +
  theme_classic()
ggplot(cn_nc_dup_no_dup_sheep, aes(x = stat.dup, y = stat.nodup)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Stats with Duplicates", y = "Test Stats without Duplicates", title = "CN Vs. NC Duplicates Test Statistics Comparison - Sheep") +
  theme_classic()

# Let's label some outliers. (Well, high leverage points).
# We'll first get some diagnostics so we can objectively choose points to look
# at.

# Regress log2-Fold-Change for all regression models
cn_nc_cow_sheep_dup_model = lm(log2FoldChange.sheep ~ log2FoldChange.cow, data = cn_nc_dup_cow_sheep)
cn_nc_cow_sheep_nodup_model = lm(log2FoldChange.sheep ~ log2FoldChange.cow, data = cn_nc_nodup_cow_sheep)
cn_nc_dupstat_cow_model = lm(log2FoldChange.dup ~ log2FoldChange.nodup, data = cn_nc_dup_no_dup_cow)
cn_nc_dupstat_sheep_model = lm(log2FoldChange.dup ~ log2FoldChange.nodup, data = cn_nc_dup_no_dup_sheep)

summary(cn_nc_cow_sheep_dup_model)
summary(cn_nc_cow_sheep_nodup_model)
summary(cn_nc_dupstat_cow_model)
summary(cn_nc_dupstat_sheep_model)

x11()
par(mfrow = c(2, 2))
plot(cn_nc_cow_sheep_dup_model)
plot(cn_nc_cow_sheep_nodup_model)
plot(cn_nc_dupstat_cow_model)
plot(cn_nc_dupstat_sheep_model)
dev.off()

# There appear to be some problems that are deviations from the assumptions of
# linear regression. That's fine -  we weren't interested in an actual
# regression model and there is no reason to expect a linear relationship
# between the log2-fold-change in cows and in sheep. We want to find influential
# points.

# What we want are the high influence points. We can look at this in multiple
# ways: Cook's distance and DFBETAs. Let's look at these.
par(mfrow = c(1, 1))
plot(cn_nc_cow_sheep_dup_model, 4, id.n = 4) # 1960, 3215, 8017, 9180
plot(cn_nc_cow_sheep_dup_model, 5)
plot(cn_nc_cow_sheep_nodup_model, 4, id.n = 4) # 2131, 3674, 3768, 6875
plot(cn_nc_cow_sheep_nodup_model, 5)
plot(cn_nc_dupstat_cow_model, 4) # 1020, 1604, 2591
plot(cn_nc_dupstat_cow_model, 5)
plot(cn_nc_dupstat_sheep_model, 4) # 1578, 2571, 2619
plot(cn_nc_dupstat_sheep_model, 5)

cookd_cn_nc_cow_sheep_dup_model = as.data.frame(cooks.distance(cn_nc_cow_sheep_dup_model))
cookd_cn_nc_cow_sheep_nodup_model = as.data.frame(cooks.distance(cn_nc_cow_sheep_nodup_model))
cookd_cn_nc_dupstat_cow_model = as.data.frame(cooks.distance(cn_nc_dupstat_cow_model))
cookd_cn_nc_dupstat_sheep_model = as.data.frame(cooks.distance(cn_nc_dupstat_sheep_model))

# Augment the data with the diagnostics for easier reference.
diag_metrics_cn_nc_cow_sheep_dup_model = augment(cn_nc_cow_sheep_dup_model)
diag_metrics_cn_nc_cow_sheep_nodup_model = augment(cn_nc_cow_sheep_nodup_model)
diag_metrics_cn_nc_dupstat_cow_model = augment(cn_nc_dupstat_cow_model)
diag_metrics_cn_nc_dupstat_sheep_model = augment(cn_nc_dupstat_sheep_model)

# Let's add in the Cook's distance.
cn_nc_dup_cow_sheep_diag_info = cbind(cn_nc_dup_cow_sheep, diag_metrics_cn_nc_cow_sheep_dup_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
cn_nc_nodup_cow_sheep_diag_info = cbind(cn_nc_nodup_cow_sheep, diag_metrics_cn_nc_cow_sheep_nodup_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
cn_nc_dup_no_dup_cow_diag_info = cbind(cn_nc_dup_no_dup_cow, diag_metrics_cn_nc_dupstat_cow_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
cn_nc_dup_no_dup_sheep_diag_info = cbind(cn_nc_dup_no_dup_sheep, diag_metrics_cn_nc_dupstat_sheep_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])

# Remake the plots.
cn_nc_dup_cow_sheep_diag_info = cn_nc_dup_cow_sheep_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
cn_nc_nodup_cow_sheep_diag_info = cn_nc_nodup_cow_sheep_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
cn_nc_dup_no_dup_cow_diag_info = cn_nc_dup_no_dup_cow_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
cn_nc_dup_no_dup_sheep_diag_info = cn_nc_dup_no_dup_sheep_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))

ggplot(cn_nc_dup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "CN Vs. NC Alignment log2-Fold-Change Comparison Across Species - With Duplicates") +
  theme_classic()
ggplot(cn_nc_nodup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "CN Vs. NC Alignment log2-Fold-Change Comparison Across Species - No Duplicates") +
  theme_classic()
ggplot(cn_nc_dup_no_dup_cow_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "CN Vs. NC Duplicates log2-Fold-Change Comparison - Cows") +
  theme_classic()
ggplot(cn_nc_dup_no_dup_sheep_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "CN Vs. NC Duplicates log2-Fold-Change Comparison - Sheep") +
  theme_classic()

# Regress the test statistics for all regression models
cn_nc_cow_sheep_dup_stat_model = lm(stat.sheep ~ stat.cow, data = cn_nc_dup_cow_sheep)
cn_nc_cow_sheep_nodup_stat_model = lm(stat.sheep ~ stat.cow, data = cn_nc_nodup_cow_sheep)
cn_nc_dupstat_cow_stat_model = lm(stat.dup ~ stat.nodup, data = cn_nc_dup_no_dup_cow)
cn_nc_dupstat_sheep_stat_model = lm(stat.dup ~ stat.nodup, data = cn_nc_dup_no_dup_sheep)

x11()
par(mfrow = c(2, 2))
plot(cn_nc_cow_sheep_dup_stat_model)
plot(cn_nc_cow_sheep_nodup_stat_model)
plot(cn_nc_dupstat_cow_stat_model)
plot(cn_nc_dupstat_sheep_stat_model)
dev.off()

diag_metrics_cn_nc_cow_sheep_dup_stat_model = augment(cn_nc_cow_sheep_dup_stat_model)
diag_metrics_cn_nc_cow_sheep_nodup_stat_model = augment(cn_nc_cow_sheep_nodup_stat_model)
diag_metrics_cn_nc_dup_stat_cow_stat_model = augment(cn_nc_dupstat_cow_stat_model)
diag_metrics_cn_nc_dup_stat_sheep_stat_model = augment(cn_nc_dupstat_sheep_stat_model)

cn_nc_dup_cow_sheep_stat_diag_info = cbind(cn_nc_dup_cow_sheep, diag_metrics_cn_nc_cow_sheep_dup_stat_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
cn_nc_nodup_cow_sheep_stat_diag_info = cbind(cn_nc_nodup_cow_sheep, diag_metrics_cn_nc_cow_sheep_nodup_stat_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
cn_nc_dup_no_dup_cow_stat_diag_info = cbind(cn_nc_dup_no_dup_cow, diag_metrics_cn_nc_dup_stat_cow_stat_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
cn_nc_dup_no_dup_sheep_stat_diag_info = cbind(cn_nc_dup_no_dup_sheep, diag_metrics_cn_nc_dup_stat_sheep_stat_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])

cn_nc_dup_cow_sheep_stat_diag_info = cn_nc_dup_cow_sheep_stat_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
cn_nc_nodup_cow_sheep_stat_diag_info = cn_nc_nodup_cow_sheep_stat_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
cn_nc_dup_no_dup_cow_stat_diag_info = cn_nc_dup_no_dup_cow_stat_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
cn_nc_dup_no_dup_sheep_stat_diag_info = cn_nc_dup_no_dup_sheep_stat_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))

ggplot(cn_nc_dup_cow_sheep_stat_diag_info, aes(x = stat.cow, y = stat.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic in Cow Alignment", y = "Test Statistic in Sheep Alignment", title = "CN Vs. NC Alignment Test Statistic Comparison Across Species - With Duplicates") +
  theme_classic()
ggplot(cn_nc_nodup_cow_sheep_stat_diag_info, aes(x = stat.cow, y = stat.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic in Cow Alignment", y = "Test Statistic in Sheep Alignment", title = "CN Vs. NC Alignment Test Statistic Comparison Across Species - No Duplicates") +
  theme_classic()
ggplot(cn_nc_dup_no_dup_cow_stat_diag_info, aes(x = stat.dup, y = stat.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic with Duplicates", y = "Test Statistic without Duplicates", title = "CN Vs. NC Duplicates Test Statistic Comparison - Cows") +
  theme_classic()
ggplot(cn_nc_dup_no_dup_sheep_stat_diag_info, aes(x = stat.dup, y = stat.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic with Duplicates", y = "Test Statistic without Duplicates", title = "CN Vs. NC Duplicates Test Statistic Comparison - Sheep") +
  theme_classic()

# 4-16-2023
# Let's also make some scatterplots of the genes that appeared in the paper.
fdr_sig_genes_paper = c("TEAD1",
                        "AEBP2",
                        "KHDRBS1",
                        "CUX1",
                        "NCOA1",
                        "PUM1",
                        "ASH1L",
                        "SPAST",
                        "ZRANB2",
                        "NUFIP2",
                        "KPNB1",
                        "HNRNPR",
                        "ZNF207",
                        "KPNA3",
                        "MEF2A",
                        "HNRNPD",
                        "CTNNB1",
                        "ZFP106",
                        "SYNCRIP",
                        "USP34",
                        "RPS6KA3",
                        "ETF1",
                        "DHX15",
                        "CAMK2D",
                        "FXR1",
                        "ZBTB44",
                        "PPP3CA",
                        "HNRNPA2B1",
                        "RBM39",
                        "PSIP1",
                        "ARID2",
                        "TRIM55",
                        "CUL2",
                        "PLCB1")

# Now, make an indicator for these so we can look at these again.
cn_nc_dup_cow_sheep_diag_info$is_claimed_fdr_sig_paper = ifelse(cn_nc_dup_cow_sheep_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
cn_nc_nodup_cow_sheep_diag_info$is_claimed_fdr_sig_paper = ifelse(cn_nc_nodup_cow_sheep_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
cn_nc_dup_no_dup_cow_diag_info$is_claimed_fdr_sig_paper = ifelse(cn_nc_dup_no_dup_cow_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
cn_nc_dup_no_dup_sheep_diag_info$is_claimed_fdr_sig_paper = ifelse(cn_nc_dup_no_dup_sheep_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)

cn_nc_dup_cow_sheep_diag_info$is_claimed_fdr_sig_paper = as.factor(cn_nc_dup_cow_sheep_diag_info$is_claimed_fdr_sig_paper)
cn_nc_nodup_cow_sheep_diag_info$is_claimed_fdr_sig_paper = as.factor(cn_nc_nodup_cow_sheep_diag_info$is_claimed_fdr_sig_paper)
cn_nc_dup_no_dup_cow_diag_info$is_claimed_fdr_sig_paper = as.factor(cn_nc_dup_no_dup_cow_diag_info$is_claimed_fdr_sig_paper)
cn_nc_dup_no_dup_sheep_diag_info$is_claimed_fdr_sig_paper = as.factor(cn_nc_dup_no_dup_sheep_diag_info$is_claimed_fdr_sig_paper)

cn_nc_dup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper = ifelse(cn_nc_dup_cow_sheep_stat_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
cn_nc_nodup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper = ifelse(cn_nc_nodup_cow_sheep_stat_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
cn_nc_dup_no_dup_cow_stat_diag_info$is_claimed_fdr_sig_paper = ifelse(cn_nc_dup_no_dup_cow_stat_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
cn_nc_dup_no_dup_sheep_stat_diag_info$is_claimed_fdr_sig_paper = ifelse(cn_nc_dup_no_dup_sheep_stat_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)

cn_nc_dup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper = as.factor(cn_nc_dup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper)
cn_nc_nodup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper = as.factor(cn_nc_nodup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper)
cn_nc_dup_no_dup_cow_stat_diag_info$is_claimed_fdr_sig_paper = as.factor(cn_nc_dup_no_dup_cow_stat_diag_info$is_claimed_fdr_sig_paper)
cn_nc_dup_no_dup_sheep_stat_diag_info$is_claimed_fdr_sig_paper = as.factor(cn_nc_dup_no_dup_sheep_stat_diag_info$is_claimed_fdr_sig_paper)

table(cn_nc_dup_cow_sheep_diag_info$is_claimed_fdr_sig_paper, useNA = "ifany")
table(cn_nc_nodup_cow_sheep_diag_info$is_claimed_fdr_sig_paper, useNA = "ifany")

cn_nc_dup_cow_sheep_diag_info_claimed_sig = subset(cn_nc_dup_cow_sheep_diag_info, is_claimed_fdr_sig_paper == 1)
cn_nc_nodup_cow_sheep_diag_info_claimed_sig = subset(cn_nc_nodup_cow_sheep_diag_info, is_claimed_fdr_sig_paper == 1)
cn_nc_dup_no_dup_cow_diag_info_claimed_sig = subset(cn_nc_dup_no_dup_cow_diag_info, is_claimed_fdr_sig_paper == 1)
cn_nc_dup_no_dup_sheep_diag_info_claimed_sig = subset(cn_nc_dup_no_dup_sheep_diag_info, is_claimed_fdr_sig_paper == 1)

cn_nc_dup_cow_sheep_stat_diag_info_claimed_sig = subset(cn_nc_dup_cow_sheep_stat_diag_info, is_claimed_fdr_sig_paper == 1)
cn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig = subset(cn_nc_nodup_cow_sheep_stat_diag_info, is_claimed_fdr_sig_paper == 1)
cn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig = subset(cn_nc_dup_no_dup_cow_stat_diag_info, is_claimed_fdr_sig_paper == 1)
cn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig = subset(cn_nc_dup_no_dup_sheep_stat_diag_info, is_claimed_fdr_sig_paper == 1)

cn_nc_dup_cow_sheep_diag_info_claimed_sig$gene_label = cn_nc_dup_cow_sheep_diag_info_claimed_sig$V1
cn_nc_nodup_cow_sheep_diag_info_claimed_sig$gene_label = cn_nc_nodup_cow_sheep_diag_info_claimed_sig$V1
cn_nc_dup_no_dup_cow_diag_info_claimed_sig$gene_label = cn_nc_dup_no_dup_cow_diag_info_claimed_sig$V1
cn_nc_dup_no_dup_sheep_diag_info_claimed_sig$gene_label = cn_nc_dup_no_dup_sheep_diag_info_claimed_sig$V1

cn_nc_dup_cow_sheep_stat_diag_info_claimed_sig$gene_label = cn_nc_dup_cow_sheep_stat_diag_info_claimed_sig$V1
cn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig$gene_label = cn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig$V1
cn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig$gene_label = cn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig$V1
cn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig$gene_label = cn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig$V1

claimed_sig_cn_nc_dup_cow_sheep_lm = lm(log2FoldChange.sheep ~ log2FoldChange.cow, data = cn_nc_dup_cow_sheep_diag_info_claimed_sig)
claimed_sig_cn_nc_nodup_cow_sheep_lm = lm(log2FoldChange.sheep ~ log2FoldChange.cow, data = cn_nc_nodup_cow_sheep_diag_info_claimed_sig)
claimed_sig_cn_nc_dup_no_dup_cow_lm = lm(log2FoldChange.nodup ~ log2FoldChange.dup, data = cn_nc_dup_no_dup_cow_diag_info_claimed_sig)
claimed_sig_cn_nc_dup_no_dup_sheep_lm = lm(log2FoldChange.nodup ~ log2FoldChange.dup, data = cn_nc_dup_no_dup_sheep_diag_info_claimed_sig)

claimed_sig_cn_nc_dup_cow_sheep_stat_lm = lm(stat.sheep ~ stat.cow, data = cn_nc_dup_cow_sheep_stat_diag_info_claimed_sig)
claimed_sig_cn_nc_nodup_cow_sheep_stat_lm = lm(stat.sheep ~ stat.cow, data = cn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig)
claimed_sig_cn_nc_dup_no_dup_cow_stat_lm = lm(stat.nodup ~ stat.dup, data = cn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig)
claimed_sig_cn_nc_dup_no_dup_sheep_stat_lm = lm(stat.nodup ~ stat.dup, data = cn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig)

claimed_sig_cn_nc_dup_cow_sheep_lm_diag_info = augment(claimed_sig_cn_nc_dup_cow_sheep_lm)
claimed_sig_cn_nc_nodup_cow_sheep_lm_diag_info = augment(claimed_sig_cn_nc_nodup_cow_sheep_lm)
claimed_sig_cn_nc_dup_no_dup_cow_lm_diag_info = augment(claimed_sig_cn_nc_dup_no_dup_cow_lm)
claimed_sig_cn_nc_dup_no_dup_sheep_lm_diag_info = augment(claimed_sig_cn_nc_dup_no_dup_sheep_lm)

claimed_sig_cn_nc_dup_cow_sheep_stat_lm_diag_info = augment(claimed_sig_cn_nc_dup_cow_sheep_stat_lm)
claimed_sig_cn_nc_nodup_cow_sheep_stat_lm_diag_info = augment(claimed_sig_cn_nc_nodup_cow_sheep_stat_lm)
claimed_sig_cn_nc_dup_no_dup_cow_stat_lm_diag_info = augment(claimed_sig_cn_nc_dup_no_dup_cow_stat_lm)
claimed_sig_cn_nc_dup_no_dup_sheep_stat_lm_diag_info = augment(claimed_sig_cn_nc_dup_no_dup_sheep_stat_lm)

colnames(claimed_sig_cn_nc_dup_cow_sheep_lm_diag_info) = paste0(colnames(claimed_sig_cn_nc_dup_cow_sheep_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_cn_nc_nodup_cow_sheep_lm_diag_info) = paste0(colnames(claimed_sig_cn_nc_nodup_cow_sheep_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_cn_nc_dup_no_dup_cow_lm_diag_info) = paste0(colnames(claimed_sig_cn_nc_dup_no_dup_cow_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_cn_nc_dup_no_dup_sheep_lm_diag_info) = paste0(colnames(claimed_sig_cn_nc_dup_no_dup_sheep_lm_diag_info), ".claimed_sig")

colnames(claimed_sig_cn_nc_dup_cow_sheep_stat_lm_diag_info) = paste0(colnames(claimed_sig_cn_nc_dup_cow_sheep_stat_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_cn_nc_nodup_cow_sheep_stat_lm_diag_info) = paste0(colnames(claimed_sig_cn_nc_nodup_cow_sheep_stat_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_cn_nc_dup_no_dup_cow_stat_lm_diag_info) = paste0(colnames(claimed_sig_cn_nc_dup_no_dup_cow_stat_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_cn_nc_dup_no_dup_sheep_stat_lm_diag_info) = paste0(colnames(claimed_sig_cn_nc_dup_no_dup_sheep_stat_lm_diag_info), ".claimed_sig")

cn_nc_dup_cow_sheep_diag_info_claimed_sig2 = cbind.data.frame(cn_nc_dup_cow_sheep_diag_info_claimed_sig, claimed_sig_cn_nc_dup_cow_sheep_lm_diag_info)
cn_nc_nodup_cow_sheep_diag_info_claimed_sig2 = cbind.data.frame(cn_nc_nodup_cow_sheep_diag_info_claimed_sig, claimed_sig_cn_nc_nodup_cow_sheep_lm_diag_info)
cn_nc_dup_no_dup_cow_diag_info_claimed_sig2 = cbind.data.frame(cn_nc_dup_no_dup_cow_diag_info_claimed_sig, claimed_sig_cn_nc_dup_no_dup_cow_lm_diag_info)
cn_nc_dup_no_dup_sheep_diag_info_claimed_sig2 = cbind.data.frame(cn_nc_dup_no_dup_sheep_diag_info_claimed_sig, claimed_sig_cn_nc_dup_no_dup_sheep_lm_diag_info)

cn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2 = cbind.data.frame(cn_nc_dup_cow_sheep_stat_diag_info_claimed_sig, claimed_sig_cn_nc_dup_cow_sheep_stat_lm_diag_info)
cn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2 = cbind.data.frame(cn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig, claimed_sig_cn_nc_nodup_cow_sheep_stat_lm_diag_info)
cn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2 = cbind.data.frame(cn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig, claimed_sig_cn_nc_dup_no_dup_cow_stat_lm_diag_info)
cn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2 = cbind.data.frame(cn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig, claimed_sig_cn_nc_dup_no_dup_sheep_stat_lm_diag_info)

cn_nc_dup_cow_sheep_diag_info_claimed_sig2$gene_label = ifelse(cn_nc_dup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig %in% cn_nc_dup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig[order(cn_nc_dup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], cn_nc_dup_cow_sheep_diag_info_claimed_sig2$V1, "")
cn_nc_nodup_cow_sheep_diag_info_claimed_sig2$gene_label = ifelse(cn_nc_nodup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig %in% cn_nc_nodup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig[order(cn_nc_nodup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], cn_nc_nodup_cow_sheep_diag_info_claimed_sig2$V1, "")
cn_nc_dup_no_dup_cow_diag_info_claimed_sig2$gene_label = ifelse(cn_nc_dup_no_dup_cow_diag_info_claimed_sig2$.cooksd.claimed_sig %in% cn_nc_dup_no_dup_cow_diag_info_claimed_sig2$.cooksd.claimed_sig[order(cn_nc_dup_no_dup_cow_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], cn_nc_dup_no_dup_cow_diag_info_claimed_sig2$V1, "")
cn_nc_dup_no_dup_sheep_diag_info_claimed_sig2$gene_label = ifelse(cn_nc_dup_no_dup_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig %in% cn_nc_dup_no_dup_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig[order(cn_nc_dup_no_dup_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], cn_nc_dup_no_dup_sheep_diag_info_claimed_sig2$V1, "")

cn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2$gene_label = ifelse(cn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig %in% cn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig[order(cn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], cn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2$V1, "")
cn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2$gene_label = ifelse(cn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig %in% cn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig[order(cn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], cn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2$V1, "")
cn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2$gene_label = ifelse(cn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2$.cooksd.claimed_sig %in% cn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2$.cooksd.claimed_sig[order(cn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], cn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2$V1, "")
cn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2$gene_label = ifelse(cn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig %in% cn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig[order(cn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], cn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2$V1, "")

# Make the new plots.
ggplot(cn_nc_dup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point(aes(color = is_claimed_fdr_sig_paper, shape = is_claimed_fdr_sig_paper)) +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, aes(linetype = is_claimed_fdr_sig_paper, group = is_claimed_fdr_sig_paper)) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  stat_ellipse() +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "CN Vs. NC Alignment log2-Fold-Change Comparison Across Species - With Duplicates") +
  theme_classic()
ggplot(cn_nc_dup_cow_sheep_diag_info_claimed_sig2, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "CN Vs. NC Alignment log2-Fold-Change Comparison Across Species - With Duplicates") +
  theme_classic()


ggplot(cn_nc_nodup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep, color = is_claimed_fdr_sig_paper, shape = is_claimed_fdr_sig_paper)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  stat_ellipse() +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "CN Vs. NC Alignment log2-Fold-Change Comparison Across Species - No Duplicates") +
  theme_classic()
ggplot(cn_nc_nodup_cow_sheep_diag_info_claimed_sig2, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "CN Vs. NC Alignment log2-Fold-Change Comparison Across Species - No Duplicates") +
  theme_classic()

ggplot(cn_nc_dup_no_dup_cow_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup, color = is_claimed_fdr_sig_paper, shape = is_claimed_fdr_sig_paper)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  stat_ellipse() +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "CN Vs. NC Duplicates log2-Fold-Change Comparison - Cows") +
  theme_classic()
ggplot(cn_nc_dup_no_dup_cow_diag_info_claimed_sig2, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "CN Vs. NC Duplicates log2-Fold-Change Comparison - Cows") +
  theme_classic()

ggplot(cn_nc_dup_no_dup_sheep_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup, color = is_claimed_fdr_sig_paper, shape = is_claimed_fdr_sig_paper)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  stat_ellipse() +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "CN Vs. NC Duplicates log2-Fold-Change Comparison - Sheep") +
  theme_classic()
ggplot(cn_nc_dup_no_dup_sheep_diag_info_claimed_sig2, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "CN Vs. NC Duplicates log2-Fold-Change Comparison - Sheep") +
  theme_classic()

"PDE4D" %in% cn_nc_dup_cow$V1
"PDE4D" %in% cn_nc_dup_sheep$V1

cn_nc_dup_sheep$V1[which(cn_nc_dup_sheep$sig == "sig")]
cn_nc_nodup_sheep$V1[which(cn_nc_nodup_sheep$sig == "sig")]

meta_analysis_df = fread("/mnt/chromeos/removable/HASAN LEXAR/BIOSTAT 646/Project/Other_Papers_Meta_Data.csv")
str(meta_analysis_df)

cn_nc_dup_sheep$V1[which(cn_nc_dup_sheep$sig == "sig" & cn_nc_dup_sheep$V1 %in% meta_analysis_df$Gene)]

cn_nc_nodup_sheep$V1[which(cn_nc_nodup_sheep$sig == "sig")]

cn_nc_dup_cow$V1[which(cn_nc_dup_cow$sig == "sig")]
# [1] "LOC100847413" "LOC101907941"

cn_nc_nodup_cow$V1[which(cn_nc_nodup_cow$sig == "sig")]

# Check something.
cn_nc_dup_sheep$is_claimed_fdr_sig = ifelse(cn_nc_dup_sheep$V1 %in% fdr_sig_genes_paper, 1, 0)
cn_nc_nodup_sheep$is_claimed_fdr_sig = ifelse(cn_nc_nodup_sheep$V1 %in% fdr_sig_genes_paper, 1, 0)
cn_nc_dup_cow$is_claimed_fdr_sig = ifelse(cn_nc_dup_cow$V1 %in% fdr_sig_genes_paper, 1, 0)
cn_nc_nodup_cow$is_claimed_fdr_sig = ifelse(cn_nc_nodup_cow$V1 %in% fdr_sig_genes_paper, 1, 0)

table(cn_nc_dup_sheep$is_claimed_fdr_sig, cn_nc_dup_sheep$sig)
table(cn_nc_nodup_sheep$is_claimed_fdr_sig, cn_nc_nodup_sheep$sig)
table(cn_nc_dup_cow$is_claimed_fdr_sig, cn_nc_dup_cow$sig)
table(cn_nc_nodup_cow$is_claimed_fdr_sig, cn_nc_nodup_cow$sig)
