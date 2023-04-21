#' The purpose of this file is to merge the datasets provided by DESeq2 and
#' both generate scatterplots of log-2 fold change for each gene across species
#' alignment and across deduplication status, and scatterplots of Wald test
#' statistics in the same conditions. Additionally, we identify genes with the
#' five highest Cook's D in our various datasets.
#'
#' This file in particular looks at the NN vs. NC comparison.
#'
#' Author: Hasan Abu-Amara
#' Date: 4-11-2023

library(ggplot2)
library(qqman)
library(data.table)
library(tidyverse)
library(dplyr)
library(broom)
library(ggrepel)
library(ggpubr)

workspace_chromebook = FALSE
results_dir = ifelse(workspace_chromebook, "/mnt/chromeos/removable/HASAN LEXAR/BIOSTAT 646/Project/DESeq2 Katarina 4-11-2023/", "H:/BIOSTAT 646/Project/DESeq2 Katarina 4-11-2023/")
cc_nc_nodup_sheep = fread(paste0(results_dir, "CC_vs_NC_no_dup_sheep"))
cc_nc_nodup_cow = fread(paste0(results_dir, "CC_vs_NC_no_dup_cow"))
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

nn_nc_nodup_cow_sheep = merge(nn_nc_nodup_cow, nn_nc_nodup_sheep, by = "V1", suffixes = c(".cow", ".sheep"))
ggplot(nn_nc_nodup_cow_sheep, aes(x = stat.cow, y = stat.sheep)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic in Cow Alignment", y = "Test Statistic in Sheep Alignment", title = "NN Vs. NC Alignment Test Statistic Comparison Across Species - No Duplicates") +
  theme_classic()
ggplot(nn_nc_nodup_cow_sheep, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "NN Vs. NC Alignment log2-Fold-Change Comparison Across Species - No Duplicates") +
  theme_classic()

cor(nn_nc_nodup_cow_sheep$stat.cow, nn_nc_nodup_cow_sheep$stat.sheep)
cor(nn_nc_nodup_cow_sheep$log2FoldChange.cow, nn_nc_nodup_cow_sheep$log2FoldChange.sheep)

nn_nc_dup_cow_sheep = merge(nn_nc_dup_cow, nn_nc_dup_sheep, by = "V1", suffixes = c(".cow", ".sheep"))
ggplot(nn_nc_dup_cow_sheep, aes(x = stat.cow, y = stat.sheep)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic in Cow Alignment", y = "Test Statistic in Sheep Alignment", title = "NN Vs. NC Alignment Test Statistic Comparison Across Species - With Duplicates") +
  theme_classic()
ggplot(nn_nc_dup_cow_sheep, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "NN Vs. NC Alignment log2-Fold-Change Comparison Across Species - With Duplicates") +
  theme_classic()

cor(nn_nc_dup_cow_sheep$stat.cow, nn_nc_dup_cow_sheep$stat.sheep)
cor(nn_nc_dup_cow_sheep$log2FoldChange.cow, nn_nc_dup_cow_sheep$log2FoldChange.sheep)

nn_nc_dup_no_dup_cow = merge(nn_nc_dup_cow, nn_nc_nodup_cow, by = "V1", suffixes = c(".dup", ".nodup"))
nn_nc_dup_no_dup_sheep = merge(nn_nc_dup_sheep, nn_nc_nodup_sheep, by = "V1", suffixes = c(".dup", ".nodup"))
cor(nn_nc_dup_no_dup_cow$log2FoldChange.dup, nn_nc_dup_no_dup_cow$log2FoldChange.nodup)
cor(nn_nc_dup_no_dup_sheep$log2FoldChange.dup, nn_nc_dup_no_dup_sheep$log2FoldChange.nodup)
cor(nn_nc_dup_no_dup_cow$stat.dup, nn_nc_dup_no_dup_cow$stat.nodup)
cor(nn_nc_dup_no_dup_sheep$stat.dup, nn_nc_dup_no_dup_sheep$stat.nodup)

ggplot(nn_nc_dup_no_dup_cow, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "NN Vs. NC Duplicates log2-Fold-Change Comparison - Cows") +
  theme_classic()
ggplot(nn_nc_dup_no_dup_cow, aes(x = stat.dup, y = stat.nodup)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Stats with Duplicates", y = "Test Stats without Duplicates", title = "NN Vs. NC Duplicates Test Statistics Comparison - Cows") +
  theme_classic()
ggplot(nn_nc_dup_no_dup_sheep, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "NN Vs. NC Duplicates log2-Fold-Change Comparison - Sheep") +
  theme_classic()
ggplot(nn_nc_dup_no_dup_sheep, aes(x = stat.dup, y = stat.nodup)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Stats with Duplicates", y = "Test Stats without Duplicates", title = "NN Vs. NC Duplicates Test Statistics Comparison - Sheep") +
  theme_classic()

# 4-12-2023: Let's label some outliers. (Well, high leverage points).
# We'll first get some diagnostics so we can objectively choose points to look
# at.

# Regress log2-Fold-Change for all regression models
nn_nc_cow_sheep_dup_model = lm(log2FoldChange.sheep ~ log2FoldChange.cow, data = nn_nc_dup_cow_sheep)
nn_nc_cow_sheep_nodup_model = lm(log2FoldChange.sheep ~ log2FoldChange.cow, data = nn_nc_nodup_cow_sheep)
nn_nc_dupstat_cow_model = lm(log2FoldChange.dup ~ log2FoldChange.nodup, data = nn_nc_dup_no_dup_cow)
nn_nc_dupstat_sheep_model = lm(log2FoldChange.dup ~ log2FoldChange.nodup, data = nn_nc_dup_no_dup_sheep)

summary(nn_nc_cow_sheep_dup_model)
summary(nn_nc_cow_sheep_nodup_model)
summary(nn_nc_dupstat_cow_model)
summary(nn_nc_dupstat_sheep_model)

x11()
par(mfrow = c(2, 2))
plot(nn_nc_cow_sheep_dup_model)
plot(nn_nc_cow_sheep_nodup_model)
plot(nn_nc_dupstat_cow_model)
plot(nn_nc_dupstat_sheep_model)
dev.off()

# There appear to be some problems that are deviations from the assumptions of
# linear regression. That's fine -  we weren't interested in an actual
# regression model and there is no reason to expect a linear relationship
# between the log2-fold-change in cows and in sheep. We want to find influential
# points.

# What we want are the high influence points. We can look at this in multiple
# ways: Cook's distance and DFBETAs. Let's look at these.
par(mfrow = c(1, 1))
plot(nn_nc_cow_sheep_dup_model, 4, id.n = 4) # 1960, 3215, 8017, 9180
plot(nn_nc_cow_sheep_dup_model, 5)
plot(nn_nc_cow_sheep_nodup_model, 4, id.n = 4) # 2131, 3674, 3768, 6875
plot(nn_nc_cow_sheep_nodup_model, 5)
plot(nn_nc_dupstat_cow_model, 4) # 1020, 1604, 2591
plot(nn_nc_dupstat_cow_model, 5)
plot(nn_nc_dupstat_sheep_model, 4) # 1578, 2571, 2619
plot(nn_nc_dupstat_sheep_model, 5)

# Look at these genes.
nn_nc_dup_cow_sheep[c(1960, 3215, 8017, 9180),]
nn_nc_nodup_cow_sheep[c(2131, 3674, 3768, 6875),]
nn_nc_dup_no_dup_cow[c(1020, 1604, 2591),]
nn_nc_dup_no_dup_sheep[c(1578, 2571, 2619),]

cookd_nn_nc_cow_sheep_dup_model = as.data.frame(cooks.distance(nn_nc_cow_sheep_dup_model))
cookd_nn_nc_cow_sheep_nodup_model = as.data.frame(cooks.distance(nn_nc_cow_sheep_nodup_model))
cookd_nn_nc_dupstat_cow_model = as.data.frame(cooks.distance(nn_nc_dupstat_cow_model))
cookd_nn_nc_dupstat_sheep_model = as.data.frame(cooks.distance(nn_nc_dupstat_sheep_model))

# Augment the data with the diagnostics for easier reference.
diag_metrics_nn_nc_cow_sheep_dup_model = augment(nn_nc_cow_sheep_dup_model)
diag_metrics_nn_nc_cow_sheep_nodup_model = augment(nn_nc_cow_sheep_nodup_model)
diag_metrics_nn_nc_dupstat_cow_model = augment(nn_nc_dupstat_cow_model)
diag_metrics_nn_nc_dupstat_sheep_model = augment(nn_nc_dupstat_sheep_model)

# Let's add in the Cook's distance.
nn_nc_dup_cow_sheep_diag_info = cbind(nn_nc_dup_cow_sheep, diag_metrics_nn_nc_cow_sheep_dup_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
nn_nc_nodup_cow_sheep_diag_info = cbind(nn_nc_nodup_cow_sheep, diag_metrics_nn_nc_cow_sheep_nodup_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
nn_nc_dup_no_dup_cow_diag_info = cbind(nn_nc_dup_no_dup_cow, diag_metrics_nn_nc_dupstat_cow_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
nn_nc_dup_no_dup_sheep_diag_info = cbind(nn_nc_dup_no_dup_sheep, diag_metrics_nn_nc_dupstat_sheep_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])

# Remake the plots.
nn_nc_dup_cow_sheep_diag_info = nn_nc_dup_cow_sheep_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
nn_nc_nodup_cow_sheep_diag_info = nn_nc_nodup_cow_sheep_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
nn_nc_dup_no_dup_cow_diag_info = nn_nc_dup_no_dup_cow_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
nn_nc_dup_no_dup_sheep_diag_info = nn_nc_dup_no_dup_sheep_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))

ggplot(nn_nc_dup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "NN Vs. NC Alignment log2-Fold-Change Comparison Across Species - With Duplicates") +
  theme_classic()
ggplot(nn_nc_nodup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "NN Vs. NC Alignment log2-Fold-Change Comparison Across Species - No Duplicates") +
  theme_classic()
ggplot(nn_nc_dup_no_dup_cow_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "NN Vs. NC Duplicates log2-Fold-Change Comparison - Cows") +
  theme_classic()
ggplot(nn_nc_dup_no_dup_sheep_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "NN Vs. NC Duplicates log2-Fold-Change Comparison - Sheep") +
  theme_classic()

# Regress the test statistics for all regression models
nn_nc_cow_sheep_dup_stat_model = lm(stat.sheep ~ stat.cow, data = nn_nc_dup_cow_sheep)
nn_nc_cow_sheep_nodup_stat_model = lm(stat.sheep ~ stat.cow, data = nn_nc_nodup_cow_sheep)
nn_nc_dupstat_cow_stat_model = lm(stat.dup ~ stat.nodup, data = nn_nc_dup_no_dup_cow)
nn_nc_dupstat_sheep_stat_model = lm(stat.dup ~ stat.nodup, data = nn_nc_dup_no_dup_sheep)

x11()
par(mfrow = c(2, 2))
plot(nn_nc_cow_sheep_dup_stat_model)
plot(nn_nc_cow_sheep_nodup_stat_model)
plot(nn_nc_dupstat_cow_stat_model)
plot(nn_nc_dupstat_sheep_stat_model)
dev.off()

diag_metrics_nn_nc_cow_sheep_dup_stat_model = augment(nn_nc_cow_sheep_dup_stat_model)
diag_metrics_nn_nc_cow_sheep_nodup_stat_model = augment(nn_nc_cow_sheep_nodup_stat_model)
diag_metrics_nn_nc_dup_stat_cow_stat_model = augment(nn_nc_dupstat_cow_stat_model)
diag_metrics_nn_nc_dup_stat_sheep_stat_model = augment(nn_nc_dupstat_sheep_stat_model)

nn_nc_dup_cow_sheep_stat_diag_info = cbind(nn_nc_dup_cow_sheep, diag_metrics_nn_nc_cow_sheep_dup_stat_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
nn_nc_nodup_cow_sheep_stat_diag_info = cbind(nn_nc_nodup_cow_sheep, diag_metrics_nn_nc_cow_sheep_nodup_stat_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
nn_nc_dup_no_dup_cow_stat_diag_info = cbind(nn_nc_dup_no_dup_cow, diag_metrics_nn_nc_dup_stat_cow_stat_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])
nn_nc_dup_no_dup_sheep_stat_diag_info = cbind(nn_nc_dup_no_dup_sheep, diag_metrics_nn_nc_dup_stat_sheep_stat_model[, c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid")])

nn_nc_dup_cow_sheep_stat_diag_info = nn_nc_dup_cow_sheep_stat_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
nn_nc_nodup_cow_sheep_stat_diag_info = nn_nc_nodup_cow_sheep_stat_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
nn_nc_dup_no_dup_cow_stat_diag_info = nn_nc_dup_no_dup_cow_stat_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))
nn_nc_dup_no_dup_sheep_stat_diag_info = nn_nc_dup_no_dup_sheep_stat_diag_info %>%
  mutate(gene_label = ifelse(.cooksd %in% .cooksd[order(.cooksd, decreasing = TRUE)[1:5]], paste0(V1, ", ", round(.cooksd, 3)), ""))

ggplot(nn_nc_dup_cow_sheep_stat_diag_info, aes(x = stat.cow, y = stat.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic in Cow Alignment", y = "Test Statistic in Sheep Alignment", title = "NN Vs. NC Alignment Test Statistic Comparison Across Species - With Duplicates") +
  theme_classic()
ggplot(nn_nc_nodup_cow_sheep_stat_diag_info, aes(x = stat.cow, y = stat.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic in Cow Alignment", y = "Test Statistic in Sheep Alignment", title = "NN Vs. NC Alignment Test Statistic Comparison Across Species - No Duplicates") +
  theme_classic()
ggplot(nn_nc_dup_no_dup_cow_stat_diag_info, aes(x = stat.dup, y = stat.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic with Duplicates", y = "Test Statistic without Duplicates", title = "NN Vs. NC Duplicates Test Statistic Comparison - Cows") +
  theme_classic()
ggplot(nn_nc_dup_no_dup_sheep_stat_diag_info, aes(x = stat.dup, y = stat.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Test Statistic with Duplicates", y = "Test Statistic without Duplicates", title = "NN Vs. NC Duplicates Test Statistic Comparison - Sheep") +
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
nn_nc_dup_cow_sheep_diag_info$is_claimed_fdr_sig_paper = ifelse(nn_nc_dup_cow_sheep_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
nn_nc_nodup_cow_sheep_diag_info$is_claimed_fdr_sig_paper = ifelse(nn_nc_nodup_cow_sheep_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
nn_nc_dup_no_dup_cow_diag_info$is_claimed_fdr_sig_paper = ifelse(nn_nc_dup_no_dup_cow_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
nn_nc_dup_no_dup_sheep_diag_info$is_claimed_fdr_sig_paper = ifelse(nn_nc_dup_no_dup_sheep_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)

nn_nc_dup_cow_sheep_diag_info$is_claimed_fdr_sig_paper = as.factor(nn_nc_dup_cow_sheep_diag_info$is_claimed_fdr_sig_paper)
nn_nc_nodup_cow_sheep_diag_info$is_claimed_fdr_sig_paper = as.factor(nn_nc_nodup_cow_sheep_diag_info$is_claimed_fdr_sig_paper)
nn_nc_dup_no_dup_cow_diag_info$is_claimed_fdr_sig_paper = as.factor(nn_nc_dup_no_dup_cow_diag_info$is_claimed_fdr_sig_paper)
nn_nc_dup_no_dup_sheep_diag_info$is_claimed_fdr_sig_paper = as.factor(nn_nc_dup_no_dup_sheep_diag_info$is_claimed_fdr_sig_paper)

nn_nc_dup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper = ifelse(nn_nc_dup_cow_sheep_stat_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
nn_nc_nodup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper = ifelse(nn_nc_nodup_cow_sheep_stat_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
nn_nc_dup_no_dup_cow_stat_diag_info$is_claimed_fdr_sig_paper = ifelse(nn_nc_dup_no_dup_cow_stat_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)
nn_nc_dup_no_dup_sheep_stat_diag_info$is_claimed_fdr_sig_paper = ifelse(nn_nc_dup_no_dup_sheep_stat_diag_info$V1 %in% fdr_sig_genes_paper, 1, 0)

nn_nc_dup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper = as.factor(nn_nc_dup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper)
nn_nc_nodup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper = as.factor(nn_nc_nodup_cow_sheep_stat_diag_info$is_claimed_fdr_sig_paper)
nn_nc_dup_no_dup_cow_stat_diag_info$is_claimed_fdr_sig_paper = as.factor(nn_nc_dup_no_dup_cow_stat_diag_info$is_claimed_fdr_sig_paper)
nn_nc_dup_no_dup_sheep_stat_diag_info$is_claimed_fdr_sig_paper = as.factor(nn_nc_dup_no_dup_sheep_stat_diag_info$is_claimed_fdr_sig_paper)

table(nn_nc_dup_cow_sheep_diag_info$is_claimed_fdr_sig_paper, useNA = "ifany")
table(nn_nc_nodup_cow_sheep_diag_info$is_claimed_fdr_sig_paper, useNA = "ifany")

nn_nc_dup_cow_sheep_diag_info_claimed_sig = subset(nn_nc_dup_cow_sheep_diag_info, is_claimed_fdr_sig_paper == 1)
nn_nc_nodup_cow_sheep_diag_info_claimed_sig = subset(nn_nc_nodup_cow_sheep_diag_info, is_claimed_fdr_sig_paper == 1)
nn_nc_dup_no_dup_cow_diag_info_claimed_sig = subset(nn_nc_dup_no_dup_cow_diag_info, is_claimed_fdr_sig_paper == 1)
nn_nc_dup_no_dup_sheep_diag_info_claimed_sig = subset(nn_nc_dup_no_dup_sheep_diag_info, is_claimed_fdr_sig_paper == 1)

nn_nc_dup_cow_sheep_stat_diag_info_claimed_sig = subset(nn_nc_dup_cow_sheep_stat_diag_info, is_claimed_fdr_sig_paper == 1)
nn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig = subset(nn_nc_nodup_cow_sheep_stat_diag_info, is_claimed_fdr_sig_paper == 1)
nn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig = subset(nn_nc_dup_no_dup_cow_stat_diag_info, is_claimed_fdr_sig_paper == 1)
nn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig = subset(nn_nc_dup_no_dup_sheep_stat_diag_info, is_claimed_fdr_sig_paper == 1)

nn_nc_dup_cow_sheep_diag_info_claimed_sig$gene_label = nn_nc_dup_cow_sheep_diag_info_claimed_sig$V1
nn_nc_nodup_cow_sheep_diag_info_claimed_sig$gene_label = nn_nc_nodup_cow_sheep_diag_info_claimed_sig$V1
nn_nc_dup_no_dup_cow_diag_info_claimed_sig$gene_label = nn_nc_dup_no_dup_cow_diag_info_claimed_sig$V1
nn_nc_dup_no_dup_sheep_diag_info_claimed_sig$gene_label = nn_nc_dup_no_dup_sheep_diag_info_claimed_sig$V1

nn_nc_dup_cow_sheep_stat_diag_info_claimed_sig$gene_label = nn_nc_dup_cow_sheep_stat_diag_info_claimed_sig$V1
nn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig$gene_label = nn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig$V1
nn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig$gene_label = nn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig$V1
nn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig$gene_label = nn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig$V1

claimed_sig_nn_nc_dup_cow_sheep_lm = lm(log2FoldChange.sheep ~ log2FoldChange.cow, data = nn_nc_dup_cow_sheep_diag_info_claimed_sig)
claimed_sig_nn_nc_nodup_cow_sheep_lm = lm(log2FoldChange.sheep ~ log2FoldChange.cow, data = nn_nc_nodup_cow_sheep_diag_info_claimed_sig)
claimed_sig_nn_nc_dup_no_dup_cow_lm = lm(log2FoldChange.nodup ~ log2FoldChange.dup, data = nn_nc_dup_no_dup_cow_diag_info_claimed_sig)
claimed_sig_nn_nc_dup_no_dup_sheep_lm = lm(log2FoldChange.nodup ~ log2FoldChange.dup, data = nn_nc_dup_no_dup_sheep_diag_info_claimed_sig)

claimed_sig_nn_nc_dup_cow_sheep_stat_lm = lm(stat.sheep ~ stat.cow, data = nn_nc_dup_cow_sheep_stat_diag_info_claimed_sig)
claimed_sig_nn_nc_nodup_cow_sheep_stat_lm = lm(stat.sheep ~ stat.cow, data = nn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig)
claimed_sig_nn_nc_dup_no_dup_cow_stat_lm = lm(stat.nodup ~ stat.dup, data = nn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig)
claimed_sig_nn_nc_dup_no_dup_sheep_stat_lm = lm(stat.nodup ~ stat.dup, data = nn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig)

claimed_sig_nn_nc_dup_cow_sheep_lm_diag_info = augment(claimed_sig_nn_nc_dup_cow_sheep_lm)
claimed_sig_nn_nc_nodup_cow_sheep_lm_diag_info = augment(claimed_sig_nn_nc_nodup_cow_sheep_lm)
claimed_sig_nn_nc_dup_no_dup_cow_lm_diag_info = augment(claimed_sig_nn_nc_dup_no_dup_cow_lm)
claimed_sig_nn_nc_dup_no_dup_sheep_lm_diag_info = augment(claimed_sig_nn_nc_dup_no_dup_sheep_lm)

claimed_sig_nn_nc_dup_cow_sheep_stat_lm_diag_info = augment(claimed_sig_nn_nc_dup_cow_sheep_stat_lm)
claimed_sig_nn_nc_nodup_cow_sheep_stat_lm_diag_info = augment(claimed_sig_nn_nc_nodup_cow_sheep_stat_lm)
claimed_sig_nn_nc_dup_no_dup_cow_stat_lm_diag_info = augment(claimed_sig_nn_nc_dup_no_dup_cow_stat_lm)
claimed_sig_nn_nc_dup_no_dup_sheep_stat_lm_diag_info = augment(claimed_sig_nn_nc_dup_no_dup_sheep_stat_lm)

colnames(claimed_sig_nn_nc_dup_cow_sheep_lm_diag_info) = paste0(colnames(claimed_sig_nn_nc_dup_cow_sheep_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_nn_nc_nodup_cow_sheep_lm_diag_info) = paste0(colnames(claimed_sig_nn_nc_nodup_cow_sheep_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_nn_nc_dup_no_dup_cow_lm_diag_info) = paste0(colnames(claimed_sig_nn_nc_dup_no_dup_cow_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_nn_nc_dup_no_dup_sheep_lm_diag_info) = paste0(colnames(claimed_sig_nn_nc_dup_no_dup_sheep_lm_diag_info), ".claimed_sig")

colnames(claimed_sig_nn_nc_dup_cow_sheep_stat_lm_diag_info) = paste0(colnames(claimed_sig_nn_nc_dup_cow_sheep_stat_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_nn_nc_nodup_cow_sheep_stat_lm_diag_info) = paste0(colnames(claimed_sig_nn_nc_nodup_cow_sheep_stat_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_nn_nc_dup_no_dup_cow_stat_lm_diag_info) = paste0(colnames(claimed_sig_nn_nc_dup_no_dup_cow_stat_lm_diag_info), ".claimed_sig")
colnames(claimed_sig_nn_nc_dup_no_dup_sheep_stat_lm_diag_info) = paste0(colnames(claimed_sig_nn_nc_dup_no_dup_sheep_stat_lm_diag_info), ".claimed_sig")

nn_nc_dup_cow_sheep_diag_info_claimed_sig2 = cbind.data.frame(nn_nc_dup_cow_sheep_diag_info_claimed_sig, claimed_sig_nn_nc_dup_cow_sheep_lm_diag_info)
nn_nc_nodup_cow_sheep_diag_info_claimed_sig2 = cbind.data.frame(nn_nc_nodup_cow_sheep_diag_info_claimed_sig, claimed_sig_nn_nc_nodup_cow_sheep_lm_diag_info)
nn_nc_dup_no_dup_cow_diag_info_claimed_sig2 = cbind.data.frame(nn_nc_dup_no_dup_cow_diag_info_claimed_sig, claimed_sig_nn_nc_dup_no_dup_cow_lm_diag_info)
nn_nc_dup_no_dup_sheep_diag_info_claimed_sig2 = cbind.data.frame(nn_nc_dup_no_dup_sheep_diag_info_claimed_sig, claimed_sig_nn_nc_dup_no_dup_sheep_lm_diag_info)

nn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2 = cbind.data.frame(nn_nc_dup_cow_sheep_stat_diag_info_claimed_sig, claimed_sig_nn_nc_dup_cow_sheep_stat_lm_diag_info)
nn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2 = cbind.data.frame(nn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig, claimed_sig_nn_nc_nodup_cow_sheep_stat_lm_diag_info)
nn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2 = cbind.data.frame(nn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig, claimed_sig_nn_nc_dup_no_dup_cow_stat_lm_diag_info)
nn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2 = cbind.data.frame(nn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig, claimed_sig_nn_nc_dup_no_dup_sheep_stat_lm_diag_info)

nn_nc_dup_cow_sheep_diag_info_claimed_sig2$gene_label = ifelse(nn_nc_dup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig %in% nn_nc_dup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig[order(nn_nc_dup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], nn_nc_dup_cow_sheep_diag_info_claimed_sig2$V1, "")
nn_nc_nodup_cow_sheep_diag_info_claimed_sig2$gene_label = ifelse(nn_nc_nodup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig %in% nn_nc_nodup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig[order(nn_nc_nodup_cow_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], nn_nc_nodup_cow_sheep_diag_info_claimed_sig2$V1, "")
nn_nc_dup_no_dup_cow_diag_info_claimed_sig2$gene_label = ifelse(nn_nc_dup_no_dup_cow_diag_info_claimed_sig2$.cooksd.claimed_sig %in% nn_nc_dup_no_dup_cow_diag_info_claimed_sig2$.cooksd.claimed_sig[order(nn_nc_dup_no_dup_cow_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], nn_nc_dup_no_dup_cow_diag_info_claimed_sig2$V1, "")
nn_nc_dup_no_dup_sheep_diag_info_claimed_sig2$gene_label = ifelse(nn_nc_dup_no_dup_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig %in% nn_nc_dup_no_dup_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig[order(nn_nc_dup_no_dup_sheep_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], nn_nc_dup_no_dup_sheep_diag_info_claimed_sig2$V1, "")

nn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2$gene_label = ifelse(nn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig %in% nn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig[order(nn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], nn_nc_dup_cow_sheep_stat_diag_info_claimed_sig2$V1, "")
nn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2$gene_label = ifelse(nn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig %in% nn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig[order(nn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], nn_nc_nodup_cow_sheep_stat_diag_info_claimed_sig2$V1, "")
nn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2$gene_label = ifelse(nn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2$.cooksd.claimed_sig %in% nn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2$.cooksd.claimed_sig[order(nn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], nn_nc_dup_no_dup_cow_stat_diag_info_claimed_sig2$V1, "")
nn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2$gene_label = ifelse(nn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig %in% nn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig[order(nn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2$.cooksd.claimed_sig)[1:5]], nn_nc_dup_no_dup_sheep_stat_diag_info_claimed_sig2$V1, "")


# Make the new plots.
ggplot(nn_nc_dup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point(aes(color = is_claimed_fdr_sig_paper, shape = is_claimed_fdr_sig_paper)) +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, aes(linetype = is_claimed_fdr_sig_paper, group = is_claimed_fdr_sig_paper)) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  stat_ellipse() +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "NN Vs. NC Alignment log2-Fold-Change Comparison Across Species - With Duplicates") +
  theme_classic()
ggplot(nn_nc_dup_cow_sheep_diag_info_claimed_sig2, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "NN Vs. NC Alignment log2-Fold-Change Comparison Across Species - With Duplicates") +
  theme_classic()


ggplot(nn_nc_nodup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep, color = is_claimed_fdr_sig_paper, shape = is_claimed_fdr_sig_paper)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  stat_ellipse() +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "NN Vs. NC Alignment log2-Fold-Change Comparison Across Species - No Duplicates") +
  theme_classic()
ggplot(nn_nc_nodup_cow_sheep_diag_info_claimed_sig2, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  labs(x = "log2-Fold-Change in Cow Alignment", y = "log2-Fold-Change in Sheep Alignment", title = "NN Vs. NC Alignment log2-Fold-Change Comparison Across Species - No Duplicates") +
  theme_classic()

ggplot(nn_nc_dup_no_dup_cow_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup, color = is_claimed_fdr_sig_paper, shape = is_claimed_fdr_sig_paper)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  stat_ellipse() +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "NN Vs. NC Duplicates log2-Fold-Change Comparison - Cows") +
  theme_classic()
ggplot(nn_nc_dup_no_dup_cow_diag_info_claimed_sig2, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "NN Vs. NC Duplicates log2-Fold-Change Comparison - Cows") +
  theme_classic()

ggplot(nn_nc_dup_no_dup_sheep_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup, color = is_claimed_fdr_sig_paper, shape = is_claimed_fdr_sig_paper)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  stat_ellipse() +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "NN Vs. NC Duplicates log2-Fold-Change Comparison - Sheep") +
  theme_classic()
ggplot(nn_nc_dup_no_dup_sheep_diag_info_claimed_sig2, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point() +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change with Duplicates", y = "log2-Fold-Change without Duplicates", title = "NN Vs. NC Duplicates log2-Fold-Change Comparison - Sheep") +
  theme_classic()

"PDE4D" %in% nn_nc_dup_cow$V1
"PDE4D" %in% nn_nc_dup_sheep$V1

nn_nc_dup_sheep$V1[which(nn_nc_dup_sheep$sig == "sig")]
nn_nc_nodup_sheep$V1[which(nn_nc_nodup_sheep$sig == "sig")]

meta_analysis_df = fread("/mnt/chromeos/removable/HASAN LEXAR/BIOSTAT 646/Project/Other_Papers_Meta_Data.csv")
str(meta_analysis_df)

nn_nc_dup_sheep$V1[which(nn_nc_dup_sheep$sig == "sig" & nn_nc_dup_sheep$V1 %in% meta_analysis_df$Gene)]

nn_nc_nodup_sheep$V1[which(nn_nc_nodup_sheep$sig == "sig")]

nn_nc_dup_cow$V1[which(nn_nc_dup_cow$sig == "sig")]
# [1] "LOC100847413" "LOC101907941"

nn_nc_nodup_cow$V1[which(nn_nc_nodup_cow$sig == "sig")]

# Check something.
nn_nc_dup_sheep$is_claimed_fdr_sig = ifelse(nn_nc_dup_sheep$V1 %in% fdr_sig_genes_paper, 1, 0)
nn_nc_nodup_sheep$is_claimed_fdr_sig = ifelse(nn_nc_nodup_sheep$V1 %in% fdr_sig_genes_paper, 1, 0)
nn_nc_dup_cow$is_claimed_fdr_sig = ifelse(nn_nc_dup_cow$V1 %in% fdr_sig_genes_paper, 1, 0)
nn_nc_nodup_cow$is_claimed_fdr_sig = ifelse(nn_nc_nodup_cow$V1 %in% fdr_sig_genes_paper, 1, 0)

table(nn_nc_dup_sheep$is_claimed_fdr_sig, nn_nc_dup_sheep$sig)
table(nn_nc_nodup_sheep$is_claimed_fdr_sig, nn_nc_nodup_sheep$sig)
table(nn_nc_dup_cow$is_claimed_fdr_sig, nn_nc_dup_cow$sig)
table(nn_nc_nodup_cow$is_claimed_fdr_sig, nn_nc_nodup_cow$sig)

# 4-20-2023: Let's make one figure with ALL EIGHT PANELS.

l2fc_fig1 = ggplot(nn_nc_dup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point(alpha = 0.1) +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471, ylim = c(NA, 23)) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change (Cow)", y = "log2-Fold-Change\n(Sheep)") +
  theme_classic()
l2fc_fig2 = ggplot(nn_nc_nodup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
  geom_point(alpha = 0.1) +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471, ylim = c(-4, 7), xlim = c(NA, 6)) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change (Cow)", y = "log2-Fold-Change\n(Sheep)") +
  theme_classic()
l2fc_fig3 = ggplot(nn_nc_dup_no_dup_cow_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point(alpha = 0.1) +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change (Duplicated)", y = "log2-Fold-Change\n(Deduplicated)") +
  theme_classic()
l2fc_fig4 = ggplot(nn_nc_dup_no_dup_sheep_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
  geom_point(alpha = 0.1) +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "log2-Fold-Change (Duplicated)", y = "log2-Fold-Change\n(Deduplicated)") +
  theme_classic()
l2fc_fig5 = ggplot(nn_nc_dup_cow_sheep_stat_diag_info, aes(x = stat.cow, y = stat.sheep)) +
  geom_point(alpha = 0.1) +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471, ylim = c(-2.5, 8), xlim = c(NA, 3)) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Wald Statistic (Cow)", y = "Wald Statistic\n(Sheep)") +
  theme_classic()
l2fc_fig6 = ggplot(nn_nc_nodup_cow_sheep_stat_diag_info, aes(x = stat.cow, y = stat.sheep)) +
  geom_point(alpha = 0.1) +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471, ylim = c(NA, 7)) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Wald Statistic (Cow)", y = "Wald Statistic\n(Sheep)") +
  theme_classic()
l2fc_fig7 = ggplot(nn_nc_dup_no_dup_cow_stat_diag_info, aes(x = stat.dup, y = stat.nodup)) +
  geom_point(alpha = 0.1) +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Wald Statistic (Duplicated)", y = "Wald Statistic\n(Deduplicated)") +
  theme_classic()
l2fc_fig8 = ggplot(nn_nc_dup_no_dup_sheep_stat_diag_info, aes(x = stat.dup, y = stat.nodup)) +
  geom_point(alpha = 0.1) +
  geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471) +
  geom_smooth(method = lm, se = FALSE) +
  labs(x = "Wald Statistic (Duplicated)", y = "Wald Statistic\n(Deduplicated)") +
  theme_classic()

x11()
ggarrange(l2fc_fig1,
          l2fc_fig2,
          l2fc_fig3,
          l2fc_fig4,
          l2fc_fig5,
          l2fc_fig6,
          l2fc_fig7,
          l2fc_fig8,
          labels = c("A - Sheep Vs. Cow - Duplicated",
                     "B - Sheep Vs. Cow - Deduplicated",
                     "C - Deduplicated Vs. Duplicated - Cow",
                     "D - Deduplicated Vs. Duplicated - Sheep",
                     "E - Sheep Vs. Cow - Duplicated",
                     "F - Sheep Vs. Cow - Deduplicated",
                     "G - Deduplicated Vs. Duplicated - Cow",
                     "H - Deduplicated Vs. Duplicated - Sheep"),
          ncol = 2,
          nrow = 4,
          hjust = 0,
          vjust = 0.8,
          font.label = list(size = 12, color = "black", face = "bold", family = NULL))

ggsave("Figure ALPHA NN vs NC Scatterplots 2 Test 3.jpeg",
       device = "jpeg",
       path = "H:/BIOSTAT 646/Project/",
       width = 16.7,
       height = 8.51,
       dpi = 1200)

# l2fc_fig1 = ggplot(nn_nc_dup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
#   geom_point() +
#   geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471, ylim = c(NA, 23)) +
#   geom_smooth(method = lm, se = FALSE) +
#   labs(x = "log2-Fold-Change (Cow)", y = "log2-Fold-Change\n(Sheep)") +
#   theme_classic()
# l2fc_fig1_2 = ggscatter(nn_nc_dup_cow_sheep_diag_info,
#                       x = "log2FoldChange.cow",
#                       y = "log2FoldChange.sheep",
#                       label = "gene_label",
#                       repel = TRUE,
#                       add = "reg.line",
#                       xlab = "log2-Fold-Change (Cow)",
#                       ylab = "log2-Fold-Change\n(Sheep)",
#                       title = "Sheep Vs. Cow - Duplicated")
# l2fc_fig2 = ggplot(nn_nc_nodup_cow_sheep_diag_info, aes(x = log2FoldChange.cow, y = log2FoldChange.sheep)) +
#   geom_point() +
#   geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471, ylim = c(-4, 7), xlim = c(NA, 6)) +
#   geom_smooth(method = lm, se = FALSE) +
#   labs(x = "log2-Fold-Change (Cow)", y = "log2-Fold-Change\n(Sheep)") +
#   theme_classic()
# l2fc_fig2_2 = ggscatter(nn_nc_nodup_cow_sheep_diag_info,
#                         x = "log2FoldChange.cow",
#                         y = "log2FoldChange.sheep",
#                         label = "gene_label",
#                         repel = TRUE,
#                         add = "reg.line",
#                         xlab = "log2-Fold-Change (Cow)",
#                         ylab = "log2-Fold-Change\n(Sheep)",
#                         title = "Sheep Vs. Cow - Deduplicated")
# l2fc_fig3 = ggplot(nn_nc_dup_no_dup_cow_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
#   geom_point() +
#   geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471) +
#   geom_smooth(method = lm, se = FALSE) +
#   labs(x = "log2-Fold-Change (Duplicated)", y = "log2-Fold-Change\n(Deduplicated)") +
#   theme_classic()
# l2fc_fig3_2 = ggscatter(nn_nc_dup_no_dup_cow_diag_info,
#                         x = "log2FoldChange.dup",
#                         y = "log2FoldChange.nodup",
#                         label = "gene_label",
#                         repel = TRUE,
#                         add = "reg.line",
#                         xlab = "log2-Fold-Change (Duplicated)",
#                         ylab = "log2-Fold-Change\n(Deduplicated)",
#                         title = "Deduplicated Vs. Duplicated - Cow")
# l2fc_fig4 = ggplot(nn_nc_dup_no_dup_sheep_diag_info, aes(x = log2FoldChange.dup, y = log2FoldChange.nodup)) +
#   geom_point() +
#   geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471) +
#   geom_smooth(method = lm, se = FALSE) +
#   labs(x = "log2-Fold-Change (Duplicated)", y = "log2-Fold-Change\n(Deduplicated)") +
#   theme_classic()
# l2fc_fig4_2 = ggscatter(nn_nc_dup_no_dup_sheep_diag_info,
#                         x = "log2FoldChange.dup",
#                         y = "log2FoldChange.nodup",
#                         label = "gene_label",
#                         repel = TRUE,
#                         add = "reg.line",
#                         xlab = "log2-Fold-Change (Duplicated)",
#                         ylab = "log2-Fold-Change\n(Deduplicated)",
#                         title = "Deduplicated Vs. Duplicated - Sheep")
# l2fc_fig5 = ggplot(nn_nc_dup_cow_sheep_stat_diag_info, aes(x = stat.cow, y = stat.sheep)) +
#   geom_point() +
#   geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471, ylim = c(-2.5, 8), xlim = c(NA, 3)) +
#   geom_smooth(method = lm, se = FALSE) +
#   labs(x = "Wald Statistic (Cow)", y = "Wald Statistic\n(Sheep)") +
#   theme_classic()
# l2fc_fig5_2 = ggscatter(nn_nc_dup_cow_sheep_diag_info,
#                         x = "stat.cow",
#                         y = "stat.sheep",
#                         label = "gene_label",
#                         repel = TRUE,
#                         add = "reg.line",
#                         xlab = "Wald Statistic (Cow)",
#                         ylab = "Wald Statistic\n(Sheep)",
#                         title = "Sheep Vs. Cow - Duplicated")
# l2fc_fig6 = ggplot(nn_nc_nodup_cow_sheep_stat_diag_info, aes(x = stat.cow, y = stat.sheep)) +
#   geom_point() +
#   geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471, ylim = c(NA, 7)) +
#   geom_smooth(method = lm, se = FALSE) +
#   labs(x = "Wald Statistic (Cow)", y = "Wald Statistic\n(Sheep)") +
#   theme_classic()
# l2fc_fig6_2 = ggscatter(nn_nc_nodup_cow_sheep_diag_info,
#                         x = "stat.cow",
#                         y = "stat.sheep",
#                         label = "gene_label",
#                         repel = TRUE,
#                         add = "reg.line",
#                         xlab = "Wald Statistic (Cow)",
#                         ylab = "Wald Statistic\n(Sheep)",
#                         title = "Sheep Vs. Cow - Deduplicated")
# l2fc_fig7 = ggplot(nn_nc_dup_no_dup_cow_stat_diag_info, aes(x = stat.dup, y = stat.nodup)) +
#   geom_point() +
#   geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471) +
#   geom_smooth(method = lm, se = FALSE) +
#   labs(x = "Wald Statistic (Duplicated)", y = "Wald Statistic\n(Deduplicated)") +
#   theme_classic()
# l2fc_fig7_2 = ggscatter(nn_nc_dup_no_dup_cow_diag_info,
#                         x = "stat.dup",
#                         y = "stat.nodup",
#                         label = "gene_label",
#                         repel = TRUE,
#                         add = "reg.line",
#                         xlab = "Wald Statistic (Duplicated)",
#                         ylab = "Wald Statistic\n(Deduplicated)",
#                         title = "Deduplicated Vs. Duplicated - Cow")
# l2fc_fig8 = ggplot(nn_nc_dup_no_dup_sheep_stat_diag_info, aes(x = stat.dup, y = stat.nodup)) +
#   geom_point() +
#   geom_text_repel(aes(label = gene_label), box.padding = 2, max.overlaps = Inf, force = 1, max.time = 1, seed = 29471) +
#   geom_smooth(method = lm, se = FALSE) +
#   labs(x = "Wald Statistic (Duplicated)", y = "Wald Statistic\n(Deduplicated)") +
#   theme_classic()
# l2fc_fig8_2 = ggscatter(nn_nc_dup_no_dup_sheep_diag_info,
#                         x = "stat.dup",
#                         y = "stat.nodup",
#                         label = "gene_label",
#                         repel = TRUE,
#                         add = "reg.line",
#                         xlab = "Wald Statistic (Duplicated)",
#                         ylab = "Wald Statistic\n(Deduplicated)",
#                         title = "Deduplicated Vs. Duplicated - Sheep")
#
# x11()
# ggarrange(l2fc_fig1_2,
#           l2fc_fig2_2,
#           l2fc_fig3_2,
#           l2fc_fig4_2,
#           l2fc_fig5_2,
#           l2fc_fig6_2,
#           l2fc_fig7_2,
#           l2fc_fig8_2,
#           labels = LETTERS[1:8],
#           ncol = 2,
#           nrow = 4)
#           # hjust = 0,
#           # vjust = 0.8,
#           # font.label = list(size = 12, color = "black", face = "bold", family = NULL))
