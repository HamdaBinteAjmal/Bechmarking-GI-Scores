## Automate
#Input Datasets:
## Dataset


## To simply, first we need to organise the gene1 and gene2 location. 
## if gene1 > gene2, keep order else swap

library(readxl)
library(tidyverse)
library(tidylog)
library(RColorBrewer) # for heatmap colors
library(kableExtra) # for formatting kables
library(corrr)
library(tibble)
library(tidyr)
library(dplyr)
nc_genes = read.csv(file = "InputData/Dede/pan-species-control-nonessentials-50genes.txt", sep = "\t")
nc_genes = unname(unlist(nc_genes,recursive = TRUE))

pc_genes = read.csv(file = "InputData/Dede/pan-species-control-essentials-50genes.txt",sep = "\t")
pc_genes = unname(unlist(pc_genes,recursive = TRUE))
dede = read_delim("InputData/Dede/counts.txt", delim = "\t")


dede = dede %>%
  mutate(pgRNA_id = row_number()) %>%
  separate(GENE_CLONE, into = c("gene1", "gRNA1_seq",  "gene2", "gRNA2_seq"), sep = "_", remove = TRUE, extra = "merge") %>%
  select(-c( GENE)) %>%
   
  mutate(gene1 = ifelse(gene1 %in% nc_genes, "NT", gene1),
         gene2 = if_else(gene2 %in% nc_genes, "NT", gene2))  %>%
  mutate(g1_ = ifelse(gene1 < gene2, gene1, gene2),
         g2_ = ifelse(gene1 < gene2, gene2, gene1)) %>%
  mutate(gRNA1_seq_ = ifelse(gene1 < gene2, gRNA1_seq, gRNA2_seq),
         gRNA2_seq_ = ifelse(gene1 < gene2, gRNA2_seq, gRNA1_seq)) %>%
  mutate(gene1 = g1_, gene2 = g2_, gRNA1_seq = gRNA1_seq_, gRNA2_seq = gRNA2_seq_) %>%
  select(-c(g1_, g2_, gRNA1_seq_, gRNA2_seq_)) %>%
  
  mutate(
    gene1_is_essential = gene1 %in% pc_genes,
    gene2_is_essential = gene2 %in% pc_genes
  )  %>%
  mutate(target_type = NA) %>% 
  mutate(target_type = case_when(
    gene1 == "NT" & gene2 == "NT" ~ "NC",
    ((gene1 == "NT") != (gene2 == "NT")) & (gene2_is_essential | gene1_is_essential) ~ "PC",
    (gene1 != "NT" & gene2 !="NT") & (gene1 != gene2) ~ "DT",
    TRUE ~ "ST"))  %>% 
  #discarded. pgRNAs with < 2 reads per million (RPM) in the plasmid pool or with a read count of zero at any time point were also removed.
  
  mutate(plasmid_RPM = plasmid.T0.Ex/sum(plasmid.T0.Ex)*10^6) %>%
  filter(plasmid_RPM >= 2) %>% 
  select(-plasmid_RPM) %>%
  mutate(across(where(is.numeric), ~na_if(.x, 0))) %>%

  mutate(A549.T2A.Ex = log2(A549.T2A.Ex / plasmid.T0.Ex),
         A549.T2B.Ex = log2(A549.T2B.Ex / plasmid.T0.Ex),
         A549.T2C.Ex = log2(A549.T2C.Ex / plasmid.T0.Ex),
         HT29.T2A.Ex = log2(HT29.T2A.Ex / plasmid.T0.Ex),
         HT29.T2B.Ex = log2(HT29.T2B.Ex / plasmid.T0.Ex),
         HT29.T2C.Ex = log2(HT29.T2C.Ex / plasmid.T0.Ex),
         
         OVCAR8.T2A.Ex = log2(OVCAR8.T2A.Ex / plasmid.T0.Ex),
         OVCAR8.T2B.Ex = log2(OVCAR8.T2B.Ex / plasmid.T0.Ex),
         OVCAR8.T2C.Ex = log2(OVCAR8.T2C.Ex / plasmid.T0.Ex)
         
         
         
         )%>% 
  mutate(A549= rowMeans(select(., A549.T2A.Ex, A549.T2B.Ex, A549.T2C.Ex), na.rm = TRUE)) %>%
  mutate(HT29= rowMeans(select(., HT29.T2A.Ex, HT29.T2B.Ex, HT29.T2C.Ex), na.rm = TRUE)) %>%
  mutate(OVCAR8= rowMeans(select(., OVCAR8.T2A.Ex, OVCAR8.T2B.Ex, OVCAR8.T2C.Ex), na.rm = TRUE)) %>%
  
  mutate(gene1_is_essential = gene1 %in% pc_genes,
         gene2_is_essential = gene2 %in% pc_genes) %>%
  mutate(target_type = NA) %>% 
  mutate(target_type = case_when(
    gene1 == "NT" & gene2 == "NT" ~ "NC",
    ((gene1 == "NT") != (gene2 == "NT")) & (gene2_is_essential | gene1_is_essential) ~ "PC",
    (gene1 != "NT" & gene2 !="NT") & (gene1 != gene2) ~ "DT",
    TRUE ~ "ST")) %>%
  mutate(paralog_pair = paste0(gene1, "_", gene2)) %>%
  mutate(pgRNA_id = row_number()) %>%
  select(c(pgRNA_id,gRNA1_seq,gRNA2_seq,gene1, gene2,target_type,paralog_pair,
           "A549", "OVCAR8", "HT29", gene2_is_essential, gene1_is_essential  ))




  
## No NC in this dataset
dede_long = dede %>%  pivot_longer(
      cols = starts_with("A549") | starts_with("HT29") | starts_with("OVCAR8"), # Select columns that start with PC9 or HeLa
      names_to = "rep",  # Column where the old column names will be stored
      values_to = "LFC"   # Column where the values will be stored
    ) #%>%
  #mutate(sample = str_extract(rep, "^[^.]*")) 

    
d.reps <- dede_long %>%
  ungroup() %>%
  dplyr::select(rep) %>%
  distinct()


d.lfc_rep_cor <- dede_long %>%
  ungroup() %>%
  dplyr::select(pgRNA_id, rep, LFC) %>%
  pivot_wider(names_from = rep,
              values_from = LFC) 
# d.lfc_rep_cor

d.cor <- d.lfc_rep_cor %>%
  dplyr::select(-pgRNA_id) %>%
  corrr::correlate() %>%
  shave() %>%
  stretch() %>%
  filter(!is.na(r)) %>%
  unite(c(x, y), col = "comparison", sep = "_vs_", remove = FALSE) %>%
  rename("sample1" = x, "sample2" = y)

comparisons <- d.cor %>% pull(comparison)

results <- lapply(comparisons, function(i){
  print(i)
  
  sample1 <- d.cor %>% filter(comparison == i) %>% pull(sample1)
  print(sample1)
  
  sample2 <- d.cor %>% filter(comparison == i) %>% pull(sample2)
  print(sample2)
  
  d.comparison_lfc <- d.lfc_rep_cor %>%
    dplyr::select(sample1, sample2) %>%
    rename("sample1" = sample1, "sample2" = sample2) %>%
    mutate("comparison" = i)
})

d.lfc_rep_cor_plot <- bind_rows(results)
contour_palette <- colorRampPalette(brewer.pal(n = 9, name ="Spectral"))(50)

d.lfc_rep_cor_plot %>%
  ggplot(aes(x = sample1, y = sample2)) + 
  geom_point(size = 1) +
  geom_density2d(aes(color = ..level..), alpha = 0.7) +
  geom_text(data = d.cor, 
            mapping = aes(x = -Inf, y = Inf, label = paste0("r=", round(r, 3))), 
            hjust = -0.25, vjust = 1.75, size = 3.5) +
  scale_colour_gradientn(colors = rev(contour_palette)) +
  labs(x = "sample1_LFC",
       y = "sample2_LFC", 
       color = "density") +
  facet_wrap(~comparison)




## There are no double non targetting in Dede 
## So this is not goning work

d.control_group_medians <- dede_long %>%
  group_by(rep, target_type ) %>%
  filter(target_type == "NC" | target_type == "PC") %>%
  summarize(median_lfc = median(LFC, na.rm = TRUE))
d.control_group_medians

df_LFC <- dede_long

df_LFC %>%
    filter(rep == "A549") %>%
  
  ggplot( aes(x = target_type, y = LFC , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "A549") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


## Now load the data if the genes are expressed or not. 

TPMs = read.csv("InputData/TPM.csv")

TPMs = TPMs %>% filter(X == "ACH-000681" | X == "ACH-000552" | X == "ACH-000696") %>%

  column_to_rownames(var = "X") %>%  # 'X' is the name of the new column
  t() %>%
  as.data.frame() %>%
  rownames_to_column("GENE") %>%
  mutate(GENE= sub("\\.{2}.*$", "", GENE)) %>% 
  rename(OVCAR8  = `ACH-000696` ) %>%
  rename(A549 = `ACH-000681`) %>%
  rename(HT29 = `ACH-000552`) %>%
  mutate(across(where(is.numeric), ~ .x >= log2(2+1))) %>% 
  pivot_longer(
  cols = starts_with("OVCAR8") | starts_with("A549") | starts_with("HT29"), # Select columns that start with PC9 or HeLa
  names_to = "rep",  # Column where the old column names will be stored
  values_to = "is_expressed"   # Column where the values will be stored
) #%>%



df_LFC = df_LFC %>%
  left_join(TPMs, by = c("gene1" = "GENE", "rep" = "rep")) %>%
  rename(gene1_expressed_flag = is_expressed ) %>%
  left_join(TPMs, by = c("gene2" = "GENE", "rep" = "rep")) %>%
  rename(gene2_expressed_flag = is_expressed ) 


df_LFC = df_LFC %>%
  mutate(gene1_expressed_flag = ifelse(is.na(gene1_expressed_flag) & gene1 != "NT", TRUE, gene1_expressed_flag))%>%
  mutate(gene2_expressed_flag = ifelse(is.na(gene2_expressed_flag) & gene2 != "NT", TRUE, gene2_expressed_flag)) 
# we assuming that if NA, means its expressed.

### Expression

#Since the pgPEN library uses non-targeting controls, 
#we adjusted for the fact that single-targeting pgRNAs generate only
#two double-strand breaks (1 per allele), whereas the double-targeting pgRNAs 
#generate four DSBs. To do this, we set the median (adjusted) LFC for 
#unexpressed genes of each group to zero. 

d.lfc_annot_adj_single <- df_LFC %>% # Dont need PC as no such pairs , all 
  filter(target_type == "PC" | target_type == "ST") %>%
  ## make a flag variable to indicate which pgRNAs are targeting unexpressed
  ## single targets
  mutate(unexpressed_ctrl_flag = case_when(
    gene2 == "NT" & gene1_expressed_flag == FALSE ~ TRUE, # for NAs, we assuming its expressed
    gene1 == "NT" & gene2_expressed_flag == FALSE ~ TRUE,
    TRUE ~ FALSE 
  )) %>%
  group_by(rep) %>%
  mutate(lfc_adj2 = LFC) %>% # as this adjustments was not applied before
  mutate(lfc_adj3 = lfc_adj2 - median(lfc_adj2[unexpressed_ctrl_flag == TRUE], na.rm = TRUE))


d.lfc_annot_adj_single_summary <- d.lfc_annot_adj_single %>%
  group_by(rep, unexpressed_ctrl_flag) %>%
  summarize(median = median(lfc_adj3, na.rm = TRUE))
d.lfc_annot_adj_single_summary

d.lfc_annot_adj_double <- df_LFC %>%
  filter(target_type == "DT") %>%
  ## make a flag variable to indicate which pgRNAs are targeting double
  ## unexpressed targets
  mutate(unexpressed_ctrl_flag = case_when(
    gene1_expressed_flag == FALSE & gene2_expressed_flag == FALSE ~ TRUE,
    TRUE ~ FALSE)) %>%
  mutate(lfc_adj2 = LFC) %>% # as this adjustments was not applied before
  
  group_by(rep) %>%
  mutate(lfc_adj3 = lfc_adj2 - median(lfc_adj2[unexpressed_ctrl_flag == TRUE],na.rm = TRUE))

d.lfc_annot_adj_double_summary <- d.lfc_annot_adj_double %>%
  group_by(rep, unexpressed_ctrl_flag) %>%
  summarize(median = median(lfc_adj3, na.rm = TRUE))
d.lfc_annot_adj_double_summary


### ntc_ntc
d.lfc_annot_adj_control <- df_LFC %>%
  filter(target_type == "NC") %>%
  mutate(lfc_adj2 = LFC) %>% # as this adjustments was not applied before
  
  mutate(lfc_adj3 = lfc_adj2)


d.lfc_annot_adj_single <- d.lfc_annot_adj_single %>%
  dplyr::select(-unexpressed_ctrl_flag)

### double targeting
d.lfc_annot_adj_double <- d.lfc_annot_adj_double %>%
  dplyr::select(-unexpressed_ctrl_flag)

## bind rows
d.lfc_annot_adj_pgRNA <- bind_rows(d.lfc_annot_adj_double, d.lfc_annot_adj_single, d.lfc_annot_adj_control)

d.lfc_annot_adj_pgRNA %>%
  filter(rep == "A549") %>%
  
  ggplot( aes(x = target_type, y = LFC , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "A549") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

d.lfc_annot_adj_pgRNA %>%
  filter(rep == "A549") %>%
  
  ggplot( aes(x = target_type, y = lfc_adj3 , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "A549") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


d.lfc_annot_adj_pgRNA <- d.lfc_annot_adj_pgRNA %>%
  ungroup() %>%
  #dplyr::select(pgRNA_id, rep, paralog_pair, lfc_adj3, target_type:gene2_essential_flag) %>%
  ## rename final adjusted column to CRISPR_score
  rename(CRISPR_score = lfc_adj3) %>%
  select(-c("lfc_adj2"  , "gene1_expressed_flag", "gene2_expressed_flag" ,"LFC", "gene1_is_essential", "gene2_is_essential"))



d.lfc_annot_adj_pgRNA %>% ungroup() %>% group_by(rep, target_type) %>%
  summarise(m = median( CRISPR_score,na.rm = TRUE))


d.lfc_pgRNA <- d.lfc_annot_adj_pgRNA
d.mean_single_target_CS <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "PC" | target_type == "ST") %>%
  mutate(targeting_gRNA_seq = case_when(
    gene2 == "NT" ~ gRNA1_seq,
    gene1 == "NT" ~ gRNA2_seq
  )) %>%
  mutate(pgRNA_target = ifelse(gene1 == "NT", gene2, gene1)) %>%
  group_by(rep, pgRNA_target, targeting_gRNA_seq) %>%
  mutate(mean_single_target_CS = mean(CRISPR_score, na.rm = TRUE)) %>%
  dplyr::select(rep, targeting_gRNA_seq, mean_single_target_CS, pgRNA_target) %>%
  distinct(rep, targeting_gRNA_seq,pgRNA_target, .keep_all = TRUE) %>% 
  ungroup()



d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "DT") %>%
  dplyr::select(pgRNA_id, rep, paralog_pair, CRISPR_score, target_type,
                gRNA1_seq, gRNA2_seq, gene1, gene2)


d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_double_targeting %>%
  rename(double_target_CS = CRISPR_score) %>%
  left_join(d.mean_single_target_CS, by = c("rep", "gRNA1_seq" = "targeting_gRNA_seq")) %>%
  rename(mean_gRNA1_single_target_CS = mean_single_target_CS) %>%
  left_join(d.mean_single_target_CS, by = c("rep",  "gRNA2_seq" = "targeting_gRNA_seq")) %>%
  rename(mean_gRNA2_single_target_CS = mean_single_target_CS)

d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_double_targeting %>%
  ## filter out pgRNAs where either single-targeting mean CRISPR score is NA
  filter(!is.na(mean_gRNA1_single_target_CS) & !is.na(mean_gRNA2_single_target_CS)) %>%
  mutate(expected_CS = mean_gRNA1_single_target_CS + mean_gRNA2_single_target_CS)
### Single-targeting
## get just single-targeting pgRNAs
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "ST" | target_type == "PC") %>%
  dplyr::select(pgRNA_id, rep, paralog_pair, CRISPR_score, target_type,
                gRNA1_seq, gRNA2_seq, gene1, gene2) %>%
  mutate(targeting_gRNA_seq = case_when(
    gene2 == "NT" ~ gRNA1_seq,
    gene1 == "NT" ~ gRNA2_seq
  )) %>%
  mutate(control_gRNA_seq = case_when(
    gene2 == "NT" ~ gRNA2_seq,
    gene1 == "NT" ~ gRNA1_seq
  ))


d.mean_double_control_CS <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "NC") %>%
  pivot_longer(cols = c(gRNA1_seq, gRNA2_seq),
               names_to = "position", 
               values_to = "control_gRNA_seq") %>%
  group_by(rep, control_gRNA_seq) %>%
  mutate(mean_double_control_CS = mean(CRISPR_score, na.rm = TRUE)) %>%
  dplyr::select(rep, control_gRNA_seq, mean_double_control_CS) %>%
  distinct(rep, control_gRNA_seq, .keep_all = TRUE)

## get targeting gRNAs to add back to DF to calculate expected GI scores
d.other_single_targeting_CS <- d.lfc_pgRNA_single_targeting %>%
  dplyr::select(rep, paralog_pair, targeting_gRNA_seq, control_gRNA_seq, CRISPR_score) %>%
  rename(other_single_target_CS = CRISPR_score, other_control_seq = control_gRNA_seq)

## add back "other" single-targeting pgRNA CRISPR scores into the main DF, 
## get rid of duplicates (same target seq and same control seq)

## most of the rows here are removed because theres hardly any different targettng seq
# for the same control. So we are losing a lot of information here. 


#d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
#  left_join(d.other_single_targeting_CS, by = c("rep", "paralog_pair", "targeting_gRNA_seq")) %>%
# mutate(same_control_seq = ifelse(control_gRNA_seq == other_control_seq, TRUE, FALSE)) %>%
#mutate(other_single_target_CS = ifelse(same_control_seq, NA, other_single_target_CS)) %>%
#dplyr::select(-same_control_seq)

d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  mutate(other_single_target_CS = NA)


## also there are NO matching double control 
 # so here we will do the same
 # Expected = observed for single targetting linear model.
 # its useless to use same number of rows with other control seq in the linear model.
 # just dont
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  left_join(y = d.mean_double_control_CS, by = c("rep", "control_gRNA_seq"))




# We will simply make a linear model with expeced  = observed
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  mutate(expected_CS = CRISPR_score) %>%
  rename(single_target_CS = CRISPR_score)


d.lfc_pgRNA_single_targeting_mean <- d.lfc_pgRNA_single_targeting %>%
  ungroup() %>%
  mutate(pgRNA_target = ifelse(gene1 == "NT", gene2, gene1)) %>%
  group_by(rep, pgRNA_target) %>%
  summarize(mean_expected_CS = mean(expected_CS),
            mean_observed_CS = mean(single_target_CS)) %>%
  ungroup()
d.lfc_pgRNA_single_targeting_mean


# Perfect linear model
d.lfc_pgRNA_single_targeting_mean_lm_summary <- d.lfc_pgRNA_single_targeting_mean %>%
  group_by(rep) %>%
  group_modify(~ broom::tidy(lm(mean_observed_CS ~ mean_expected_CS, data = .x))) %>%
  dplyr::ungroup() %>%
  dplyr::select(rep, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept = "(Intercept)", slope = mean_expected_CS)
d.lfc_pgRNA_single_targeting_mean_lm_summary


d.lfc_pgRNA_single_targeting_mean %>%
  ggplot(aes(x = mean_expected_CS, y = mean_observed_CS)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "gray55") +
  labs(x = "expected_CRISPR_score", y = "observed_CRISPR_score") +
  
  facet_wrap(~rep)


ggplot(data = d.lfc_pgRNA_single_targeting, 
       aes(x = expected_CS, y = single_target_CS)) +
  geom_point() + 
  geom_density2d(aes(color = ..level..), alpha = 0.9) +
  geom_abline(data = d.lfc_pgRNA_single_targeting_mean_lm_summary,
              aes(slope = slope, intercept = intercept), 
              color = "gray55", size = 1) +
  scale_colour_gradientn(colors = rev(contour_palette)) +
  labs(x = "expected_CRISPR_score", 
       y = "observed_CRISPR_score",
       color = "density") +
  facet_wrap(~rep)



d.lfc_pgRNA_single_targeting_GI <- d.lfc_pgRNA_single_targeting %>%
  left_join(d.lfc_pgRNA_single_targeting_mean_lm_summary, by = "rep") %>%
  mutate(GI_score = single_target_CS - (intercept + slope * expected_CS)) 

### Calculate GI score for double-targeting pgRNAs

d.lfc_pgRNA_double_targeting_GI <- d.lfc_pgRNA_double_targeting %>%
  left_join(d.lfc_pgRNA_single_targeting_mean_lm_summary, by = "rep") %>%
  mutate(GI_score = double_target_CS - (intercept + slope * expected_CS)) # almost same as observed - expected


## reformat to match
d.GI_scores_pgRNA_double <- d.lfc_pgRNA_double_targeting_GI %>%
  rename(observed_CS = double_target_CS) %>%
  dplyr::select(-c(mean_gRNA1_single_target_CS, mean_gRNA2_single_target_CS)) %>%
  mutate(broad_target_type = "double_targeting") %>%
  ungroup()

d.GI_scores_pgRNA_single <- d.lfc_pgRNA_single_targeting_GI %>%
  rename(observed_CS = single_target_CS) %>%
  dplyr::select(-c(other_single_target_CS, mean_double_control_CS,
                   targeting_gRNA_seq, control_gRNA_seq)) %>%
  mutate(broad_target_type = "single_targeting") %>%
  ungroup()

reps <- d.GI_scores_pgRNA_double %>%
  ungroup() %>%
  distinct(rep) %>% 
  pull()



results <- lapply(reps, function(i){
  ## get a vector of GI scores for all single-targeting ("control") pgRNAs for each rep
  single_GI_scores <- d.GI_scores_pgRNA_single %>%
    filter(rep == i) %>%
    pull(GI_score)
  
  ## get double-targeting pgRNAs for this rep, do a t-test to compare the double-
  ## targeting GI scores for each paralog pair to the control vector
  d.double_GI_scores <- d.GI_scores_pgRNA_double %>%
    filter(rep == i) %>%
    group_by(paralog_pair) %>%
    mutate(p_val = t.test(x = single_GI_scores,
                          y = GI_score,
                          paired = FALSE)$p.value) 
  
  ## adjust for multiple testing using the Benjamini-Hochberg method
  d.p_val <- d.double_GI_scores %>%
    dplyr::select(paralog_pair, p_val) %>%
    arrange(p_val) %>%
    distinct(p_val, .keep_all = TRUE) 
  
  p_vals <- d.p_val %>% 
    pull(p_val)
  
  fdr_vals <- p.adjust(p_vals, method = "BH")
  
  d.fdr <- tibble("fdr" = fdr_vals) %>%
    bind_cols(d.p_val) %>%
    dplyr::select(-p_val)
  
  ## add FDR values back into the double-targeting DF
  d.double_GI_scores <- left_join(d.double_GI_scores, d.fdr, by = "paralog_pair")
  
  return(d.double_GI_scores)
  
})

d.GI_scores_pgRNA_double <- bind_rows(results)
d.stats <- d.GI_scores_pgRNA_double %>%
  dplyr::select(paralog_pair, p_val, fdr) %>%
  distinct(paralog_pair, .keep_all = TRUE) 

d.GI_scores_pgRNA_single <- d.GI_scores_pgRNA_single %>%
  left_join(d.stats, by = "paralog_pair")

d.GI_scores_pgRNA <- bind_rows(d.GI_scores_pgRNA_double, d.GI_scores_pgRNA_single)

d.GI_scores_target <- d.GI_scores_pgRNA %>%
  ungroup() %>%
  # group_by(rep, pgRNA_target) %>% WARNING, major change here
  group_by(rep, paralog_pair) %>%
  
  mutate(mean_observed_CS = mean(observed_CS),
         mean_expected_CS = mean(expected_CS),
         mean_GI_score = mean(GI_score)) %>%
  #distinct(pgRNA_target, .keep_all = TRUE) %>% WARNING, major change here
  distinct(paralog_pair, .keep_all = TRUE) %>%
  
  dplyr::select(-c(pgRNA_id, contains("seq"), observed_CS, expected_CS, GI_score)) %>%
  ungroup()


d.GI_scores_target %>%
  group_by(rep) %>%
  summarize(n = n())

ggplot(d.GI_scores_target, aes(x = mean_expected_CS, y = mean_observed_CS, color = broad_target_type)) +
  geom_point() + 
  geom_abline(data = d.lfc_pgRNA_single_targeting_mean_lm_summary,
              aes(slope = slope, intercept = intercept),
              color = "gray55", size = 1) +
  labs(x = "expected_CRISPR_score", y = "observed_CRISPR_score") +
  
  facet_wrap(~rep)


final_df = d.GI_scores_target %>% filter(target_type == "DT")  %>% select(c(rep, paralog_pair, mean_GI_score,fdr, p_val))
final_df = final_df %>% 
  rename(GI_score = mean_GI_score )%>%
  pivot_wider(
  names_from = rep,   # Column that turns into new column names
  values_from = c(GI_score, fdr, p_val) ,   # Column that populates values in the new columns
  names_sep = "_"
)  
write.csv(final_df, file = "Parrish Score Scripts/ParrishOutput/Dede_Parrish.csv", row.names = FALSE)
