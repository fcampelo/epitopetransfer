fix_dataset_names <- function(Dataset, alias_csv){
  aliases <- read.csv(alias_csv, header = TRUE)
  for (i in seq_along(aliases$original)){
    Dataset <- ifelse(Dataset == aliases$original[i], aliases$replacement[i], Dataset)
  }
  return(Dataset)
}

fix_method_names <- function(Method, alias_csv){
  aliases <- read.csv(alias_csv, header = TRUE)
  for (i in seq_along(aliases$original)){
    Method <- ifelse(Method == aliases$original[i], aliases$replacement[i], Method)
  }
  return(factor(Method, 
                levels = aliases$replacement[aliases$order],
                labels = aliases$replacement[aliases$order],
                ordered = TRUE))
}

calc_stats <- function(mydf){
  
  mydf <- mydf %>% 
    group_by(Method) %>%
    arrange(Metric, Method, Dataset) %>%
    ungroup()
  
  methods <- as.character(unique(mydf$Method))
  methods <- methods[-which(methods=="EpitopeTransfer")]
  lapply(seq_along(unique(mydf$Metric)),
         function(i){
           tmp <- mydf %>% filter(Metric == unique(mydf$Metric)[i])
           x   <- tmp$Value[tmp$Method == "EpitopeTransfer"]
           res <- data.frame(Reference  = rep("EpitopeTransfer", length(methods)),
                             Challenger = NA_character_,
                             Metric     = unique(mydf$Metric)[i],
                             pValue = NA_real_,
                             diff     = NA_real_,
                             ci.lower = NA_real_,
                             ci.upper = NA_real_)
           for (j in seq_along(methods)){
             y <- tmp$Value[tmp$Method == methods[j]]
             cmp <- wilcox.test(x, y, paired = TRUE, conf.int = TRUE)
             res$Challenger[j] <- methods[j]
             res$pValue[j]     <- cmp$p.value
             res$diff[j]       <- cmp$estimate
             res$ci.lower[j]   <- cmp$conf.int[1]
             res$ci.upper[j]   <- cmp$conf.int[2]
           }
           res$p.adj <- p.adjust(res$pValue, method = "fdr")
           return(res)
         }) %>%
    bind_rows() %>%
    mutate(p.adj.sci = paste0("p = ", format(signif(p.adj, 3),
                                             scientific = TRUE)),
           ci = sprintf("Median of paired diffs: %2.3f (%2.3f, %2.3f)", 
                        diff, ci.lower, ci.upper))
}

est_CIs <- function(mydf){
  
  mydf <- mydf %>% 
    group_by(Method) %>%
    arrange(Metric, Method, Dataset) %>%
    ungroup()
  
  methods <- as.character(unique(mydf$Method))
  lapply(seq_along(unique(mydf$Metric)),
         function(i){
           tmp <- mydf %>% filter(Metric == unique(mydf$Metric)[i])
           res <- data.frame(Method   = rep(NA_real_, length(methods)),
                             Metric   = unique(mydf$Metric)[i],
                             median   = NA_real_,
                             ci.lower = NA_real_,
                             ci.upper = NA_real_)
           for (j in seq_along(methods)){
             x <- tmp$Value[tmp$Method == methods[j]]
             cmp <- wilcox.test(x, conf.int = TRUE)
             res$Method[j]   <- methods[j]
             res$median[j]   <- cmp$estimate
             res$ci.lower[j] <- cmp$conf.int[1]
             res$ci.upper[j] <- cmp$conf.int[2]
           }
           return(res)
         }) %>%
    bind_rows() %>%
    mutate(ci = sprintf("Median: %2.3f (%2.3f, %2.3f)", 
                        median, ci.lower, ci.upper))
}

process_epitope1d <- function(resdir, datadir, cl = 1){
  ### Process Epitope1D predictions
  ### CSV files with predictions was retrieved from Epitope1D in September 2024.
  
  dirs         <- dir(resdir, full.names = TRUE)
  dirnames     <- dir(resdir, full.names = FALSE)
  datadirs     <- dir(datadir, full.names = TRUE)
  datadirnames <- dir(datadir, full.names = FALSE)
  
  mymean <- function(z){ifelse(all(is.na(z)), NA, mean(z, na.rm = TRUE))}
  mymax <- function(z){ifelse(all(is.na(z)), NA, max(z, na.rm = TRUE))}
  
  for (i in seq_along(dirs)){
    preds <- read.csv(dir(dirs[i], pattern = "output\\_prediction.csv", full.names = TRUE), header = TRUE) %>%
      dplyr::select(Info_protein_id = Fasta_header, 
                    Info_Peptide = Peptide, Prob = Score_Epitope)
    
    df    <- readRDS(dir(datadirs[i], pattern = "^peptides.+\\.rds", full.names = TRUE)) %>%
      dplyr::select(Info_PepID, Info_protein_id, Info_peptide, Info_start_pos, Info_end_pos, Class)
    
    prots <- readRDS(dir(datadirs[i], pattern = "^proteins.+\\.rds", full.names = TRUE)) %>%
      dplyr::select(Info_protein_id, Info_sequence) %>%
      dplyr::filter(Info_protein_id %in% unique(preds$Info_protein_id))
    
    cat(sprintf("\nProcessing dir %03d/%03d (%s | %s): %03d proteins\n",
                i, length(dirs),
                dirnames[i], datadirnames[i],
                nrow(prots)))
    
    # Aggregate Epitope1D data by residue.
    x <- pbapply::pblapply(seq_along(prots$Info_protein_id), 
                           function(j, prots, df, preds){
                             require(dplyr)
                             myseq <- prots$Info_sequence[j]
                             tmp   <- preds %>%
                               filter(Info_protein_id == prots$Info_protein_id[j])
                             tmp <- tmp %>%
                               bind_cols(stringr::str_locate(myseq, tmp$Info_Peptide)) %>%
                               filter(!is.na(start))
                             
                             labels <- df %>% 
                               filter(Info_protein_id == prots$Info_protein_id[j]) %>%
                               mutate(N = nchar(Info_peptide)) %>%
                               tidyr::uncount(N) %>%
                               group_by(Info_PepID) %>%
                               mutate(Position = Info_start_pos + (1:n()) - 1) %>%
                               ungroup() %>%
                               group_by(Position) %>%
                               summarise(True_class = sign(sum(Class)), .groups = "drop") %>%
                               mutate(True_class = ifelse(True_class == 0, NA, True_class))
                             
                             x <- data.frame(Pathogen  = "", 
                                             Protein_id = prots$Info_protein_id[j],
                                             Position = 1:nchar(myseq),
                                             Predictor = "", 
                                             Predicted_prob = NA) %>%
                               left_join(labels, by = "Position") %>%
                               rowwise() %>%
                               mutate(Predicted_prob = mymax(tmp$Prob[tmp$start <= Position & tmp$end >= Position]))
                             
                             return(x)
                           }, 
                           prots = prots, df = df, preds = preds, cl = cl) %>%
      bind_rows() %>%
      dplyr::select(Pathogen, Protein_id, Position,
                    True_class, Predictor, Predicted_prob)
    
    write.csv(x, paste0(dirs[i], "/prob_by_position.csv"), 
              row.names = FALSE, quote = FALSE)
    
    
  }
  
  invisible(TRUE)
}