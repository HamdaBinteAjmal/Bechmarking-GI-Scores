library(tidyverse)
library(data.table)
library(corrr)
library(ggplot2)
library(RColorBrewer)
## From the paper, I took this data directly from Ito Paper.
## This is the list of essential genes and non essential genes.

Hart_essential <- read.delim("C:/Users/Hamda/OneDrive/Documents/GitHub/PostDoc/Conway/Research/GIScoring/Ito et al/Data/Data/CEGv2.txt", stringsAsFactors = F) %>% 
  pull(GENE)

Hart_nonessential <- read.delim("C:/Users/Hamda/OneDrive/Documents/GitHub/PostDoc/Conway/Research/GIScoring/Ito et al/Data/Data/NEGv1.txt", stringsAsFactors = F) %>% 
  pull(Gene)



ito <- readxl::read_excel("InputData/Ito/Count_data_ParalogV1.xlsx") 
ito <- ito %>%  separate(`Left-sgRNA_Right-sgRNA`, into = c("gene1_guide", "gene2_guide"), sep = "_", remove = FALSE)

rownames(ito) =ito$`Left-sgRNA_Right-sgRNA`


ito_genes = c(ito$Aureus_gene, ito$Pyogenes_gene)
ito_genes = ito_genes  %>% unique() 
ito_genes[ito_genes %in% Hart_essential]




ito = ito %>% 
  mutate(pgRNA_id = row_number()) %>%
  select(-c(`Left-sgRNA_Right-sgRNA`)) %>%
  rename(gRNA1_seq = gene1_guide ,gRNA2_seq = gene2_guide ) %>%
  rename(gene1 = Aureus_gene, gene2 =  Pyogenes_gene  ) %>%
  mutate(g1_ = ifelse(gene1 < gene2, gene1, gene2),
         g2_ = ifelse(gene1 < gene2, gene2, gene1)) %>%
  mutate(gRNA1_seq_ = ifelse(gene1 < gene2, gRNA1_seq, gRNA2_seq),
         gRNA2_seq_ = ifelse(gene1 < gene2, gRNA2_seq, gRNA1_seq)) %>%
  mutate(gene1 = g1_, gene2 = g2_, gRNA1_seq = gRNA1_seq_, gRNA2_seq = gRNA2_seq_) %>%
  select(-c(g1_, g2_, gRNA1_seq_, gRNA2_seq_)) %>%
  
  mutate(
    gene1_is_essential = gene1 %in% Hart_essential,
    gene2_is_essential = gene2 %in% Hart_essential
  )  %>%
  mutate(target_type = NA) %>% 
  mutate(target_type = case_when(
    gene1 == "AAVS1" & gene2 == "AAVS1" ~ "NC",
    ((gene1 == "AAVS1") != (gene2 == "AAVS1")) & (gene2_is_essential | gene1_is_essential) ~ "PC",
    (gene1 != "AAVS1" & gene2 !="AAVS1") & (gene1 != gene2) ~ "DT",
    TRUE ~ "ST"))  %>% 
  #discarded. pgRNAs with < 2 reads per million (RPM) in the plasmid pool or with a read count of zero at any time point were also removed.
  
  mutate(pDNA_RPM = pDNA /sum(pDNA )*10^6) %>%
  
  mutate(across(
    .cols = contains("_REP"),
    .fns = ~ifelse(pDNA_RPM < 2, NA_real_,.)  # Example transformation: add 1
  )) %>% 
  select(-c(contains("RPM"))) %>%
  
  mutate(across(
    .cols = contains("_REP"),
    .fns = ~ifelse(. == 0, NA_real_,.)  # Example transformation: add 1
  )) %>%
  mutate(across(
    .cols =contains("_REP"),
    .fns = ~log2(./pDNA)  # Example transformation: add 1
  )) %>%
  mutate(Meljuso= rowMeans(select(.,contains("Meljuso")), na.rm = TRUE)) %>%
  mutate(GI1= rowMeans(select(.,contains("GI1")), na.rm = TRUE)) %>%
  
  mutate(PK1= rowMeans(select(.,contains("PK1")), na.rm = TRUE)) %>%
  mutate(MEWO= rowMeans(select(.,contains("MEWO")), na.rm = TRUE)) %>%
  mutate(HS944T= rowMeans(select(.,contains("HS944T")), na.rm = TRUE)) %>%
  mutate(IPC298= rowMeans(select(.,contains("IPC298")), na.rm = TRUE)) %>%
  mutate(A549= rowMeans(select(.,contains("A549")), na.rm = TRUE)) %>%
  mutate(HSC5= rowMeans(select(.,contains("HSC5")), na.rm = TRUE)) %>%
  mutate(HS936T= rowMeans(select(.,contains("HS936T")), na.rm = TRUE)) %>%
  mutate(PATU8988S= rowMeans(select(.,contains("PATU8988S")), na.rm = TRUE)) %>%
  mutate(MEL202= rowMeans(select(.,contains("MEL202_003")), na.rm = TRUE)) %>%
  
  select(-c(contains("REP"))) 

## ITO does not have double negative controls
## same fate as dede

ito_long = ito %>%
  pivot_longer(
    cols = c("Meljuso","GI1","PK1","MEWO","HS944T","IPC298","A549" ,             
             "HSC5","HS936T","PATU8988S", "MEL202"), # Select columns that start with PC9 or HeLa
    names_to = "rep",  # Column where the old column names will be stored
    values_to = "LFC"   # Column where the values will be stored
  ) %>%
  filter(!is.na(LFC))
d.reps <- ito_long %>%
  ungroup() %>%
  dplyr::select(rep) %>%
  distinct()

d.lfc_rep_cor <- ito_long %>%
  ungroup() %>%
  dplyr::select(pgRNA_id, rep, LFC) %>%
  pivot_wider(names_from = rep,
              values_from = LFC) 


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
#top). Nonessential (n = 119) or pan-essential (n = 96) genes paired with sgAAVS1


d.control_group_medians <- ito_long %>%
  group_by(rep, target_type ) %>%
  filter(target_type == "NC" | target_type == "PC") %>%
  summarize(median_lfc = median(LFC, na.rm = TRUE))
d.control_group_medians

df_LFC <- ito_long

df_LFC %>%
  filter(rep == "GI1") %>%
  
  ggplot( aes(x = target_type, y = LFC , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "Meljuso") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


## There are no double non targeting in Ito 
TPMs = read.csv("InputData/TPM.csv")


TPMs = TPMs %>% filter(X == "ACH-000881" | X == "ACH-000756" | X == "ACH-000307" |
                         X == "ACH-000987" | X == "ACH-000632" | X == "ACH-000915"|
                         X=="ACH-000681" | X == "ACH-001524" | X == "ACH-000801" | X == "ACH-001554" |
                         X == "ACH-000022") %>%
  
  column_to_rownames(var = "X") %>%  # 'X' is the name of the new column
  t() %>%
  as.data.frame() %>%
  rownames_to_column("GENE") %>%
  mutate(GENE= sub("\\.{2}.*$", "", GENE)) %>% 
  rename(Meljuso  = `ACH-000881` ) %>%
  rename(GI1 = `ACH-000756`) %>%
  rename(PK1 = `ACH-000307`) %>%
  rename(MEWO      = `ACH-000987`) %>%
  rename( HS944T  = `ACH-000632`) %>%
  rename( IPC298 = `ACH-000915`) %>%
  rename(A549     = `ACH-000681`) %>%
  rename(HSC5= `ACH-001524`) %>%
  rename(HS936T   = `ACH-000801`) %>%
  rename(PATU8988S= `ACH-000022`, MEL202 = `ACH-001554`) %>%
  mutate(across(where(is.numeric), ~ .x >= log2(2+1))) %>% 
  pivot_longer(
    cols = c("A549","HS944T","PATU8988S","GI1","HSC5","IPC298","Meljuso","HS936T", "MEWO","PK1", "MEL202") , # Select columns that start with PC9 or HeLa
    names_to = "rep",  # Column where the old column names will be stored
    values_to = "is_expressed"   # Column where the values will be stored
  ) #%>%


df_LFC = df_LFC %>%
  left_join(TPMs, by = c("gene1" = "GENE", "rep" = "rep")) %>%
  rename(gene1_expressed_flag = is_expressed ) %>%
  left_join(TPMs, by = c("gene2" = "GENE", "rep" = "rep")) %>%
  rename(gene2_expressed_flag = is_expressed ) 

# Maybe not right but well..
df_LFC = df_LFC %>%
  mutate(gene1_expressed_flag = ifelse(is.na(gene1_expressed_flag) & gene1 != "AAVS1", TRUE, gene1_expressed_flag))%>%
  mutate(gene2_expressed_flag = ifelse(is.na(gene2_expressed_flag) & gene2 != "AAVS1", TRUE, gene2_expressed_flag)) # we assuming that if NA, means its expressed.


d.lfc_annot_adj_single <- df_LFC %>%  
  filter(target_type == "PC" | target_type == "ST") %>%
  ## make a flag variable to indicate which pgRNAs are targeting unexpressed
  ## single targets
  mutate(unexpressed_ctrl_flag = case_when(
    gene2 == "AAVS1" & gene1_expressed_flag == FALSE ~ TRUE, # for NAs, we assuming its expressed
    gene1 == "AAVS1" & gene2_expressed_flag == FALSE ~ TRUE,
    TRUE ~ FALSE 
  )) %>%
  group_by(rep) %>%
  mutate(lfc_adj2 = LFC) %>% # as this adjustments was not applied before
  mutate(mm = median(lfc_adj2[unexpressed_ctrl_flag == TRUE], na.rm = TRUE)) %>%
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


### ntc_ntc # No rows here
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
  rename(CRISPR_score = lfc_adj3) %>%
  select(-c("lfc_adj2"  , "gene1_expressed_flag", "gene2_expressed_flag" ,"LFC", "gene1_is_essential", "gene2_is_essential"))


d.lfc_annot_adj_pgRNA %>% ungroup() %>% group_by(rep, target_type) %>%
  summarise(m = median( CRISPR_score,na.rm = TRUE))



d.lfc_pgRNA <- d.lfc_annot_adj_pgRNA
d.mean_single_target_CS <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "PC" | target_type == "ST") %>%
  mutate(targeting_gRNA_seq = case_when(
    gene2 == "AAVS1" ~ gRNA1_seq,
    gene1 == "AAVS1" ~ gRNA2_seq
  )) %>%
  mutate(pgRNA_target = ifelse(gene1 == "AAVS1", gene2, gene1)) %>%
  group_by(rep, pgRNA_target, targeting_gRNA_seq) %>%
  mutate(mean_single_target_CS = mean(CRISPR_score, na.rm = TRUE)) %>%
  dplyr::select(rep, targeting_gRNA_seq, mean_single_target_CS, pgRNA_target) %>%
  distinct(rep, targeting_gRNA_seq,pgRNA_target, .keep_all = TRUE) %>% 
  ungroup()



d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "DT") %>% 
  mutate(paralog_pair = paste0(gene1, "_", gene2)) %>%
  dplyr::select(pgRNA_id, rep, paralog_pair, CRISPR_score, target_type,
                gRNA1_seq, gRNA2_seq, gene1, gene2)



d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_double_targeting %>%
  rename(double_target_CS = CRISPR_score) %>%
  left_join(d.mean_single_target_CS, by = c("rep", "gRNA1_seq" = "targeting_gRNA_seq")) %>%
  rename(mean_gRNA1_single_target_CS = mean_single_target_CS) %>%
  left_join(d.mean_single_target_CS, by = c("rep",  "gRNA2_seq" = "targeting_gRNA_seq")) %>%
  rename(mean_gRNA2_single_target_CS = mean_single_target_CS)

## filter out pgRNAs where either single-targeting mean CRISPR score is NA
## calculate expected double-targeting GI score by summing the two mean single-targeting
## CRISPR scores for that paralog pair


d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_double_targeting %>%
  filter(!is.na(mean_gRNA1_single_target_CS) & !is.na(mean_gRNA2_single_target_CS)) %>%
  mutate(expected_CS = mean_gRNA1_single_target_CS + mean_gRNA2_single_target_CS)
### Single-targeting
## get just single-targeting pgRNAs



d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA %>%
  ungroup() %>%
  mutate(paralog_pair = paste0(gene1, "_", gene2)) %>%
  
  filter(target_type == "ST" | target_type == "PC") %>%
  dplyr::select(pgRNA_id, rep, paralog_pair, CRISPR_score, target_type,
                gRNA1_seq, gRNA2_seq, gene1, gene2) %>%
  mutate(targeting_gRNA_seq = case_when(
    gene2 == "AAVS1" ~ gRNA1_seq,
    gene1 == "AAVS1" ~ gRNA2_seq
  )) %>%
  mutate(control_gRNA_seq = case_when(
    gene2 == "AAVS1" ~ gRNA2_seq,
    gene1 == "AAVS1" ~ gRNA1_seq
  ))


## should be NOTHING here
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

## OK so there are NO other control gRNAs (perhaps same as thompson)

## SKIP this step, so we have no ntc_ntc and also no 'other'
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  mutate(other_control_seq = NA, mean_double_control_CS  = NA) 
 # left_join(d.other_single_targeting_CS, by = c("rep", "paralog_pair", "targeting_gRNA_seq")) %>%
  #mutate(same_control_seq = ifelse(control_gRNA_seq == other_control_seq, TRUE, FALSE)) %>%
  #filter(same_control_seq == FALSE) %>%
  #dplyr::select(-same_control_seq) 
## NOW simple set expected same as observed 

d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  mutate(expected_CS = CRISPR_score) %>% # + mean_double_control_CS) %>%
  rename(single_target_CS = CRISPR_score)



d.lfc_pgRNA_single_targeting_mean <- d.lfc_pgRNA_single_targeting %>%
  ungroup() %>%
  mutate(pgRNA_target = ifelse(gene1 == "AAVS1", gene2, gene1)) %>%
  group_by(rep, pgRNA_target) %>%
  summarize(mean_expected_CS = mean(expected_CS),
            mean_observed_CS = mean(single_target_CS)) %>%
  ungroup()
d.lfc_pgRNA_single_targeting_mean


## Model with 100% corr
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


## Singele target GI score is almost always close to zero as we have no values.


d.lfc_pgRNA_single_targeting_GI <- d.lfc_pgRNA_single_targeting %>%
  left_join(d.lfc_pgRNA_single_targeting_mean_lm_summary, by = "rep") %>%
  mutate(GI_score = single_target_CS - (intercept + slope * expected_CS)) 

### Calculate GI score for double-targeting pgRNAs

d.lfc_pgRNA_double_targeting_GI <- d.lfc_pgRNA_double_targeting %>%
  left_join(d.lfc_pgRNA_single_targeting_mean_lm_summary, by = "rep") %>%
  mutate(GI_score = double_target_CS - (intercept + slope * expected_CS)) # almost same as observed - expected

# cor.test(d.lfc_pgRNA_double_targeting_GI$GI_score,d.lfc_pgRNA_double_targeting_GI$double_target_CS - d.lfc_pgRNA_double_targeting_GI$expected_CS)
#Result of this test is always 1, so use this or just do GI_score  = double_target_CS - expected_CS is just the same 
# as our intercept and slope are 0 and 1 respectively
d.GI_scores_pgRNA_single <- d.lfc_pgRNA_single_targeting_GI %>%
  rename(observed_CS = single_target_CS) %>%
  dplyr::select(-c( mean_double_control_CS,
                   targeting_gRNA_seq, control_gRNA_seq, other_control_seq)) %>%
  mutate(broad_target_type = "single_targeting") %>%
  ungroup()
d.GI_scores_pgRNA_double <- d.lfc_pgRNA_double_targeting_GI %>%
  rename(observed_CS = double_target_CS) %>%
  dplyr::select(-c(mean_gRNA1_single_target_CS, mean_gRNA2_single_target_CS)) %>%
  mutate(broad_target_type = "double_targeting") %>%
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


df_wide <- final_df %>%
  rename(GI_score = mean_GI_score) %>%
  pivot_wider(
    names_from = rep,
    values_from = c(GI_score, fdr, p_val),
    names_sep = "_"  # This adds an underscore between the rep value and the measure type for clarity
  )
write.csv(df_wide, row.names = FALSE, "Parrish Score Scripts/ParrishOutput/ito_parrish.csv")
