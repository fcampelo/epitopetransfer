library(dplyr)
library(ggplot2)
library(multcomp)
library(tidyr)
library(yardstick)
library(tidyr)
library(pROC)
library(wrappedtools)
library(stringr)
library(see)
library(ggrepel)
source("utils.R")

### Read and consolidate results data
# Get classification thresholds
intthres <- read.csv("../results/best_mcc_thresholds.csv", header = TRUE) %>%
  rename(Pathogen  = Dataset, 
         Predictor = Method)

extpreds <- c("bepipred3", "epidope", "epitopevec", "epitope1d")
extthres <- data.frame(Predictor = rep(extpreds, 
                                       each = length(unique(intthres$Pathogen))),
                       Pathogen  = rep(unique(intthres$Pathogen), times = length(extpreds)),
                       Threshold = rep(c(0.1512, 0.818, 0.5, 0.5), each = length(unique(intthres$Pathogen))))

pthres <- bind_rows(intthres, extthres)

# Read information of train/test splits
testFolds <- read.csv("../data/TrainTestFolds.csv", header = TRUE)

# Read results and add info of thresholds and data splits
X <- readRDS("../results/consolidated_results_v8.rds") %>%
  filter(!is.na(Class)) %>%
  left_join(pthres, by = c("Predictor", "Pathogen")) %>%
  left_join(testFolds, by = c("Pathogen", "Fold")) %>%
  mutate(Pred = factor(Prob >= Threshold, 
                       levels = c("TRUE", "FALSE"), 
                       labels = c("TRUE", "FALSE"),
                       ordered = TRUE),
         Class = factor(Class == 1, 
                        levels = c("TRUE", "FALSE"), 
                        labels = c("TRUE", "FALSE"),
                        ordered = TRUE)) %>%
  dplyr::select(-Threshold) %>%
  as_tibble()

# Calculate metrics using EpitopeTransfer holdout set
# Ignore warnings due to absence of positive/negative predictions in some cases.
Xmetrics.HO <- X %>%
  filter(Holdout == 1) %>%
  dplyr::select(-Holdout, -Fold) %>%
  group_by(Predictor, Pathogen) %>%
  arrange(Class) %>%
  summarise(AUC = as.numeric(pROC::roc(response = Class, 
                                       predictor = Prob, 
                                       levels = c("FALSE", "TRUE"))$auc),
            F1  = f_meas_vec(truth = Class, estimate = Pred, event_level = "first"),
            MCC = mcc_vec(truth = Class, estimate = Pred, event_level = "first"), 
            PPV = ppv_vec(truth = Class, estimate = Pred, event_level = "first"),
            NPV = npv_vec(truth = Class, estimate = Pred, event_level = "first"),
            Sensitivity = sens_vec(truth = Class, estimate = Pred, event_level = "first"),
            Specificity = spec_vec(truth = Class, estimate = Pred, event_level = "first"),
            `Balanced accuracy` = (Sensitivity + Specificity) / 2,
            Accuracy = accuracy_vec(truth = Class, estimate = Pred, event_level = "first"),
            .groups = "drop") %>%
  pivot_longer(AUC:Accuracy, names_to = "Metric", values_to = "Value") %>%
  mutate(Set = "Holdout")


# Calculate metrics using a no-leakage test set (for external baselines - results are identical in the cases of EpitopeTransfer and internal baselines)
# Ignore warnings due to absence of positive/negative predictions in some cases.
Xmetrics.noLeak <- X %>%
  filter(InTraining == 0) %>%
  dplyr::select(-Holdout, -Fold) %>%
  group_by(Predictor, Pathogen) %>%
  arrange(Class) %>%
  summarise(AUC = as.numeric(pROC::roc(response = Class, 
                                       predictor = Prob, 
                                       levels = c("FALSE", "TRUE"))$auc),
            F1  = f_meas_vec(truth = Class, estimate = Pred, event_level = "first"),
            MCC = mcc_vec(truth = Class, estimate = Pred, event_level = "first"), 
            PPV = ppv_vec(truth = Class, estimate = Pred, event_level = "first"),
            NPV = npv_vec(truth = Class, estimate = Pred, event_level = "first"),
            Sensitivity = sens_vec(truth = Class, estimate = Pred, event_level = "first"),
            Specificity = spec_vec(truth = Class, estimate = Pred, event_level = "first"),
            `Balanced accuracy` = (Sensitivity + Specificity) / 2,
            Accuracy = accuracy_vec(truth = Class, estimate = Pred, event_level = "first"),
            .groups = "drop") %>%
  pivot_longer(AUC:Accuracy, names_to = "Metric", values_to = "Value") %>%
  mutate(Set = "NoLeakage")

# Consolidate all metrics and fix dataset and method names
Xmetrics <- bind_rows(Xmetrics.HO, Xmetrics.noLeak) %>%
  rename(Dataset = Pathogen, Method = Predictor) %>%
  mutate(Dataset = fix_dataset_names(Dataset, "../data/pathogen_aliases.csv"),
         Method  = fix_method_names(Method, "../data/method_aliases.csv"),
         Value = ifelse(is.na(Value), 0, Value))

write.table(Xmetrics, 
            "../results/Performance_metrics_all.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

### Perform statistical tests:
# All-vs.one Wilcoxon signed rank (paired samples) tests with FDR correction on p-values
# Ignore warnings due to ties.
myres <- calc_stats(Xmetrics %>% filter(Set == "NoLeakage")) %>%
  mutate(Set = "NoLeakage") %>%
  bind_rows(calc_stats(Xmetrics %>% filter(Set == "Holdout")) %>%
              mutate(Set = "Holdout")) %>%
  mutate(Challenger = fix_method_names(Challenger, "../data/method_aliases.csv"))

# Save table with statistical results
stats.table <- myres %>%
  mutate(Comparison = paste0(Reference, " vs. ", Challenger),
         `Set Used` = Set) %>%
  dplyr::select(`Set Used`, Comparison, Metric, pValue, CI = ci, pVal.FDR = p.adj)

write.table(stats.table, 
            "../results/comparisons.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)


# Estimate 95% confidence intervals on median performances
# CIs estimated using the method implemented in wilcox.test
# Ignore warnings due to ties.
myest <- est_CIs(Xmetrics %>% filter(Set == "NoLeakage")) %>%
  mutate(Set = "NoLeakage") %>%
  bind_rows(est_CIs(Xmetrics %>% filter(Set == "Holdout")) %>%
              mutate(Set = "Holdout")) %>%
  mutate(Method = fix_method_names(Method, "../data/method_aliases.csv"))

# Save table with confidence intervals
est.table <- myest %>%
  rename(`Set Used` = Set) %>%
  dplyr::select(`Set Used`, Metric, Method, CI = ci)

write.table(est.table, 
            "../results/stats.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)


### ======  Generate plots
stats.gg <- myres %>%
  dplyr::select(Set, Metric, Method = Challenger, label = p.adj.sci, p.adj, pValue, diff) %>%
  mutate(colour = ifelse(p.adj < 0.05, ifelse(diff > 0, "#ff8888","#880088"), "#666666"),
         shape  = ifelse(p.adj < 0.05, ifelse(diff > 0, 24, 25), 15))

Xmetrics.summary <- Xmetrics %>%
  group_by(Metric, Method, Set) %>%
  summarise(Est = median(Value),
            se  = medianse(Value),
            ymin = min(Value),
            .groups = "drop") %>%
  left_join(stats.gg, by = c("Metric", "Method", "Set")) %>%
  mutate(label  = ifelse(Set == "Holdout", NA, label),
         colour = ifelse(Set == "Holdout", NA, colour),
         shape  = ifelse(Set == "Holdout", NA, shape),
         Est = ifelse(Set == "Holdout" & (Method %in% c("EpitopeTransfer", "ESM-1b bl.", "NPTransfer bl.")),
                      NA, Est),
         se = ifelse(Set == "Holdout" & (Method %in% c("EpitopeTransfer", "ESM-1b bl.", "NPTransfer bl.")),
                     NA, se))

for (metr in unique(Xmetrics$Metric)){
  
  df.gg <- Xmetrics %>%
    left_join(stats.gg, by = c("Metric", "Method", "Set")) %>%
    filter(Metric == metr)
  
  Xsumm.HO <- filter(Xmetrics.summary, Metric == metr, Set == "Holdout")
  Xsumm.NL <- filter(Xmetrics.summary, Metric ==metr, Set == "NoLeakage") %>%
    mutate(colour = ifelse(Method == "EpitopeTransfer", "#66bb66", colour),
           shape  = ifelse(Method == "EpitopeTransfer", 16, shape))
  
  
  df.gg %>%
    filter(Set == "NoLeakage") %>%
    ggplot(aes(x = Method, y = Value)) +
    geom_violin(aes(fill = Method), alpha = .15, linewidth = NA, 
                show.legend = FALSE) +
    geom_jitter(show.legend = FALSE,
                size = 2, pch = 20, width = .1, height = 0) +
    geom_pointrange(data = Xsumm.HO,
                    mapping = aes(y = Est, 
                                  ymin = Est - se,
                                  ymax = Est + se),
                    position = position_nudge(x = 0.1),
                    shape    = 16, 
                    linewidth = 1, 
                    size = 1,
                    colour = "#66666666") +
    geom_pointrange(data = Xsumm.NL,
                    mapping = aes(y = Est, 
                                  ymin = Est - se,
                                  ymax = Est + se),
                    show.legend = FALSE,
                    colour = Xsumm.NL$colour,
                    fill   = Xsumm.NL$colour,
                    shape  = Xsumm.NL$shape,
                    size = 1, linewidth = 1, stroke = 0) +
    #facet_wrap(Metric ~ .) +
    theme_light() +
    geom_vline(xintercept = 3.5, lty = 2, alpha = .5) +
    theme_light() + 
    theme(strip.text = element_text(face = "bold", colour = "black", size = 14),
          axis.text.x = element_text(hjust = .5, face = "bold", size = 14),
          axis.text.y = element_text(size = 12),
          axis.title  = element_text(face = "bold", size = 14),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "bottom") +
    xlab("") + ylab(metr) + 
    geom_text(data = Xsumm.NL, 
              mapping = aes(y = ymin, x = Method, label = label),
              size = 4, position = position_nudge(y = -0.04)) + 
    scale_fill_brewer(palette = "Dark2") 
  
  # ggsave(paste0("../figures/violin plots/results_violin_", metr, ".pdf"), 
  #        width = 4500, height = 2000, units = "px")
  
}



# plot AUC gains/losses of epitopetransfer vs. baselines for each dataset
gainthres <- 0.1
lossthres <- -0.05
neutralband <- c(-0.05, 0.05)
aliases <- read.csv("../data/method_aliases.csv")
tmpX <- read.table("../results/Performance_metrics_all.tsv",
                   sep = "\t", header = TRUE) %>%
  as_tibble() %>%
  filter(Set == "NoLeakage") %>%
  group_by(Metric, Dataset) %>%
  mutate(Diff = nth(Value, 4) - Value,
         Method = factor(Method,
                         levels = aliases$replacement[aliases$order],
                         labels = aliases$replacement[aliases$order],
                         ordered = TRUE)) %>%
  ungroup() %>%
  filter(!grepl("Transfer|bl", Method), Metric == "AUC") 

tmpXsumm <- tmpX  %>%
  group_by(Metric, Method) %>%
  summarise(Ngain = sum(Diff > 0),
            Nloss = sum(Diff < 0),
            GainRatio = 100 * Ngain / n(),
            meanGain = mean(Diff[Diff > 0]),
            meanLoss = mean(Diff[Diff < 0]),
            seGain   = sd(Diff[Diff > 0]) / sqrt(Ngain),
            seLoss   = sd(Diff[Diff < 0]) / sqrt(Nloss),
            .groups = "drop") %>%
  mutate(across(starts_with("se"), ~ifelse(is.na(.x), 0, .x))) %>%
  ungroup() %>% 
  pivot_longer(meanGain:meanLoss, names_to = "type", values_to = "mean") %>%
  pivot_longer(seGain:seLoss, names_to = "setype", values_to = "se") %>%
  pivot_longer(Ngain:Nloss, names_to = "gain", values_to = "N") %>%
  mutate(type = gsub("mean", "", type),
         setype = gsub("se", "", setype),
         gain   = ifelse(gain == "Ngain", "Gain", "Loss")) %>%
  filter(type == setype,
         type == gain) %>%
  dplyr::select(Metric, Method, type, mean, se, N)

tmpXsumm %>%
  ggplot(aes(x = Method, y = mean)) + 
  geom_rect(xmin = 0, xmax = 5, ymin = neutralband[1], ymax = neutralband[2],
            colour = NA, fill = "#cccccc22") +
  geom_violinhalf(data = tmpX %>% filter(Diff > 0), 
                  aes(x = Method, y = Diff, group = Method), 
                  flip = TRUE, scale = "count", linewidth = NA,
                  fill = "#bbffbbaa") +
  geom_violinhalf(data = tmpX %>% filter(Diff < 0), 
                  aes(x = Method, y = Diff, group = Method), 
                  flip = FALSE, scale = "count", linewidth = NA,
                  fill = "#ffbbbbaa") +
  geom_pointrange(aes(ymin = mean - se, ymax = mean + se, group = type), 
                  position = position_dodge(width = .15), 
                  size = 1, linewidth = 1, 
                  colour = "#44444455", stroke = 0,
                  pch = 15) +
  geom_text_repel(data = tmpX %>% group_by(Dataset) %>% filter(all(Diff > gainthres)),
                  aes(x = Method, y = abs(Diff), label = Dataset),
                  nudge_x = .1,
                  direction = "y",
                  hjust = 0,
                  min.segment.length = 0, force = 5,
                  segment.colour = "#66666699",
                  size = 4.5, box.padding = 0.5,
                  show.legend = FALSE) +
  geom_text_repel(data = tmpX %>% filter(Diff < min(neutralband)),
                  aes(x = Method, y = Diff, label = Dataset),
                  size = 4.5, , box.padding = 0.5,
                  min.segment.length = 0,  force = 5,
                  nudge_x = -0.1,
                  direction = "y",
                  hjust = 1,
                  show.legend = FALSE,
                  segment.colour = "#66666699",
                  colour = "#888888") +
  geom_point(data = tmpX %>% filter(Diff > 0), 
             aes(x = Method, y = abs(Diff), colour = (Diff > gainthres)),
             size = 2, pch = 20,
             show.legend = FALSE) + 
  geom_point(data = tmpX %>% filter(Diff < 0),
             aes(x = Method, y = Diff, colour = (Diff < -0.2)),
             size = 2, pch = 20,
             show.legend = FALSE) +
  theme_light() + 
  theme(strip.text = element_text(face = "bold", colour = "black", size = 12),
        axis.text.x = element_text(hjust = .5, face = "bold", size = 14),
        axis.text.y = element_text(size = 14),
        axis.title  = element_text(face = "bold", size = 14),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom") + 
  xlab("") + ylab("EpitopeTransfer gains in AUC") +
  scale_colour_manual(values = c("#777777", "#ff0000"))

# ggsave(paste0("../figures/violin plots/winloss_violin_AUC.pdf"), 
#        width = 4500, height = 2000, units = "px")