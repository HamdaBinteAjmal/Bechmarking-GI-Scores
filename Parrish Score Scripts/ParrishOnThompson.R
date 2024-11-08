library(readxl)
library(tidyverse)
library(tidylog)
library(RColorBrewer) # for heatmap colors
library(kableExtra) # for formatting kables
library(corrr)
library(tibble)
library(tidyr)
library(dplyr)
# Define the file path
file_path <- "InputData/Thompson/Guide_Sequences.xlsx"

# Read the first sheet of the Excel file
guide_seq <- read_excel(file_path, sheet = 1)

thompson = read.csv("InputData/Thompson/raw_counts.txt", sep = "\t")
thompson <- thompson[seq(2, nrow(thompson), by = 2), ]

thompson = thompson %>% 
  separate(gRNA, into = c("gRNA1_seq", "gRNA2_seq", "type"), sep = "_N_|_SL_oligo:", remove = FALSE, extra = "merge") %>%
  separate(GENE, into = c("gene1", "gene2"), sep = "_", remove = FALSE, extra = "merge") %>%
  mutate(gene1 = ifelse(type == "Non_targt_RND_control", "Fluc", gene1)) %>%
  mutate(gene2 = ifelse(type == "Non_targt_RND_control", "Fluc", gene2)) %>%
  merge(guide_seq, by.x = c("gRNA1_seq", "gRNA2_seq"), by.y = c("gRNA_A_IDs", "gRNA_B_IDs" )) %>% 
  rename(paralog_pair = GENE) %>%
  mutate(target_type = "DT") %>%
  mutate(target_type = ifelse(gene1 == "Fluc" & gene2 == "Fluc", "NC", target_type)) %>%
  mutate(target_type = ifelse(gene2 == "Fluc" & gene1 != "Fluc", "ST", target_type)) %>%
  mutate(target_type = ifelse(gene_pair_origin == "bagel_train_essential" & gene2 == "Fluc", "PC", target_type)) %>%
  rename(pgRNA_id =  gRNA_pair_lib_oligo_id) %>%
   select(-c(gRNA, type, gRNA_pair_lib_oligo_seq)) %>%
  mutate(
    gene1_is_essential = ifelse(gene_pair_origin == "bagel_train_essential", TRUE, FALSE),
    gene2_is_essential = FALSE #gene2 is always Fluc
  )  %>%
  mutate(RPM_R1 = Control_R1/sum(Control_R1)*10^6, RPM_R2 = Control_R2/sum(Control_R2)*10^6, RPM_R3 = Control_R3/sum(Control_R3)*10^6) %>%
  mutate(across(ends_with("R1"), ~ifelse(RPM_R1 < 2, NA_real_,.))) %>% 
  mutate(across(ends_with("R2"), ~ifelse(RPM_R2 < 2, NA_real_,.))) %>%
  mutate(across(ends_with("R3"), ~ifelse(RPM_R3 < 2, NA_real_,.))) %>%
  select(-c(contains("RPM"))) %>%
  mutate(A375_0_R1 = A375_D28_R1 ==0 | A375_D14_R1 == 0,
         A375_0_R2 = A375_D28_R2 ==0 | A375_D14_R2 == 0,
         A375_0_R3 = A375_D28_R3 ==0 | A375_D14_R3 == 0,
         
         MEWO_0_R1 = MEWO_D28_R1 ==0 | MEWO_D14_R1 == 0,
         MEWO_0_R2 = MEWO_D28_R2 ==0 | MEWO_D14_R2 == 0,
         MEWO_0_R3 = MEWO_D28_R3 ==0 | MEWO_D14_R3 == 0,
         
         RPE_0_R1 = RPE_D28_R1 ==0 | RPE_D14_R1 == 0,
         RPE_0_R2 = RPE_D28_R2 ==0 | RPE_D14_R2 == 0,
         RPE_0_R3 = RPE_D28_R3 ==0 | RPE_D14_R3 == 0) %>%
  
  mutate(
      A375_D28_R1 = if_else(A375_0_R1 | is.na(A375_0_R1),  NA_real_, A375_D28_R1),
      A375_D28_R2 = if_else(A375_0_R2 | is.na(A375_0_R2),  NA_real_, A375_D28_R2),
      A375_D28_R3 = if_else(A375_0_R3 | is.na(A375_0_R3),  NA_real_, A375_D28_R3),
    
      MEWO_D28_R1 = if_else(MEWO_0_R1| is.na(MEWO_0_R1),  NA_real_, MEWO_D28_R1),
      MEWO_D28_R2 = if_else(MEWO_0_R2| is.na(MEWO_0_R2),  NA_real_, MEWO_D28_R2),
      MEWO_D28_R3 = if_else(MEWO_0_R3| is.na(MEWO_0_R3),  NA_real_, MEWO_D28_R3),
    
      RPE_D28_R1 = if_else(RPE_0_R1| is.na(RPE_0_R1),  NA_real_, RPE_D28_R1),
      RPE_D28_R2 = if_else(RPE_0_R2 | is.na(RPE_0_R2),  NA_real_, RPE_D28_R2),
      RPE_D28_R3 = if_else(RPE_0_R3 | is.na(RPE_0_R3),  NA_real_, RPE_D28_R3)
    
    
    )%>%
  select(-c(contains("D14"), contains("_0_R"), gene_pair_origin))%>%
  #  mutate(across(ends_with("R1"), ~ifelse(RPM_R1 < 2, NA_real_,.))) %>% 

  mutate(across(ends_with("_D28_R1"), ~ log(.x / .data$Control_R1))) %>%
  mutate(across(ends_with("_D28_R2"), ~ log(.x / .data$Control_R2)))%>%
  mutate(across(ends_with("_D28_R3"), ~ log(.x / .data$Control_R3))) %>%

           

  mutate(A375= rowMeans(select(., A375_D28_R1, A375_D28_R2, A375_D28_R3), na.rm = TRUE)) %>%
  mutate(MEWO= rowMeans(select(., MEWO_D28_R1, MEWO_D28_R2, MEWO_D28_R3), na.rm = TRUE)) %>%
  mutate(RPE= rowMeans(select(., RPE_D28_R1, RPE_D28_R2, RPE_D28_R3), na.rm = TRUE)) %>%
  select(-c(Control_R1, Control_R2, Control_R3, contains("D28"))) %>%
  select(-c(gRNA1_seq, gRNA2_seq)) %>%
  rename(gRNA1_seq = gRNA_A_seq, gRNA2_seq = gRNA_B_seq)
  
  
  
thompson_long = thompson %>%
  pivot_longer(
    cols = c("A375", "MEWO", "RPE"), # Select columns that start with PC9 or HeLa
    names_to = "rep",  # Column where the old column names will be stored
    values_to = "LFC"   # Column where the values will be stored
  ) %>%
  filter(!is.na(LFC))

         
##
# The final library contained
# 41,838 different combinations of gRNAs including gRNAs targeting a total of 1191
# gene pairs, 12,803 gRNAs combined with the Fluc_gRNA control and gRNAs
# against essential/non-essential genes to aid in statistical analysis. T

# 12805 ST (with either Bagel essential or non essential)
# 498 NC (both gene1 and gene2 Fluc)
#

d.reps <- thompson_long %>%
  ungroup() %>%
  dplyr::select(rep) %>%
  distinct()

d.lfc_rep_cor <- thompson_long %>%
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
  #plot_options +
  #plot_theme +
  #theme(aspect.ratio = square_ar) +
  facet_wrap(~comparison)




d.control_group_medians <- thompson_long %>%
  group_by(rep, target_type ) %>%
  filter(target_type == "NC" | target_type == "PC") %>%
  summarize(median_lfc = median(LFC, na.rm = TRUE))
d.control_group_medians

df_LFC <- thompson_long %>%
  group_by(rep) %>%
  mutate(lfc_adj1 = LFC - median(LFC[target_type == "NC"], na.rm = TRUE),
         lfc_adj2 = lfc_adj1 / (median(lfc_adj1[target_type == "NC"], na.rm = TRUE) -
                                  median(lfc_adj1[target_type == "PC"], na.rm = TRUE)))%>%
  ungroup()


df_LFC %>%
  group_by(rep, target_type ) %>%
  filter(target_type == "NC" | target_type == "PC") %>%
  summarize(median_lfc = median(lfc_adj2, na.rm = TRUE))



df_LFC %>%
  filter(str_detect(rep, "A375")) %>%
  
  ggplot( aes(x = target_type, y = LFC , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "PC9") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

df_LFC %>%
  filter(str_detect(rep, "A375")) %>%
  
  ggplot( aes(x = target_type, y = lfc_adj2 , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "PC9") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")



TPMs = read.csv("InputData/TPM.csv")

TPMs = TPMs %>%
  
  column_to_rownames(var = "X") %>%  # 'X' is the name of the new column
  t() %>%
  as.data.frame() %>%
  rownames_to_column("GENE") %>%
  mutate(GENE= sub("\\.{2}.*$", "", GENE)) %>% 
  rename(A375  = `ACH-000219` ) %>%
  rename(MEWO = `ACH-000987`) %>%
  rename(RPE = `ACH-002464`) %>%
  select(c(GENE, A375, MEWO, RPE)) %>%
  mutate(across(where(is.numeric), ~ .x >= log2(2+1))) %>% 
  pivot_longer(
    cols = starts_with("A375") | starts_with("MEWO") | starts_with("RPE"), # Select columns that start with PC9 or HeLa
    names_to = "rep",  # Column where the old column names will be stored
    values_to = "is_expressed"   # Column where the values will be stored
  ) #%>%



df_LFC  = df_LFC %>%
  left_join(TPMs, by = c("gene1" = "GENE", "rep" = "rep")) %>%
  rename(gene1_expressed_flag = is_expressed) %>%
  left_join(TPMs, by = c("gene2" = "GENE", "rep" = "rep")) %>%
  rename(gene2_expressed_flag = is_expressed)
  
df_LFC = df_LFC %>%
  mutate(gene1_expressed_flag = ifelse(is.na(gene1_expressed_flag) & gene1 != "Fluc", TRUE, gene1_expressed_flag))%>%
  mutate(gene2_expressed_flag = ifelse(is.na(gene2_expressed_flag) & gene2 != "Fluc", TRUE, gene2_expressed_flag)) # we assuming that if NA, means its expressed.

# From Parrish et al. paper
# "Since the pgPEN library uses non-targeting controls, 
# we adjusted for the fact that single-targeting pgRNAs generate only
# two double-strand breaks (1 per allele), whereas the double-targeting pgRNAs 
# generate four DSBs. To do this, we set the median (adjusted) LFC for 
# unexpressed genes of each group to zero. "

d.lfc_annot_adj_single <- df_LFC %>% 
  filter(target_type == "PC" | target_type == "ST") %>%
  mutate(unexpressed_ctrl_flag = case_when(
    gene2 == "Fluc" & gene1_expressed_flag == FALSE ~ TRUE, # for NAs, we assuming its expressed
    gene1 == "Fluc" & gene2_expressed_flag == FALSE ~ TRUE,
    TRUE ~ FALSE 
  )) %>%
  group_by(rep) %>%
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
  group_by(rep) %>%
  mutate(lfc_adj3 = lfc_adj2 - median(lfc_adj2[unexpressed_ctrl_flag == TRUE],na.rm = TRUE))
d.lfc_annot_adj_double_summary <- d.lfc_annot_adj_double %>%
  group_by(rep, unexpressed_ctrl_flag) %>%
  summarize(median = median(lfc_adj3, na.rm = TRUE))
d.lfc_annot_adj_double_summary


### ntc_ntc
d.lfc_annot_adj_control <- df_LFC %>%
  filter(target_type == "NC") %>%
  mutate(lfc_adj3 = lfc_adj2)

### single targeting
# colnames(d.lfc_annot_adj_single)

d.lfc_annot_adj_single <- d.lfc_annot_adj_single %>%
  dplyr::select(-unexpressed_ctrl_flag)

### double targeting
d.lfc_annot_adj_double <- d.lfc_annot_adj_double %>%
  dplyr::select(-unexpressed_ctrl_flag)

## bind rows
d.lfc_annot_adj_pgRNA <- bind_rows(d.lfc_annot_adj_double, d.lfc_annot_adj_single, d.lfc_annot_adj_control)

d.lfc_annot_adj_pgRNA <- d.lfc_annot_adj_pgRNA %>%
  ungroup() %>%
  #dplyr::select(pgRNA_id, rep, paralog_pair, lfc_adj3, target_type:gene2_essential_flag) %>%
  ## rename final adjusted column to CRISPR_score
  rename(CRISPR_score = lfc_adj3) %>%
  select(-c("lfc_adj2", "LFC", "lfc_adj1" ,"gene1_expressed_flag", "gene2_expressed_flag" ,"gene1_is_essential", "gene2_is_essential"))


d.lfc_annot_adj_pgRNA %>%
  filter(str_detect(rep, "A375")) %>%
  
  ggplot( aes(x = target_type, y = CRISPR_score , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "PC9") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

d.lfc_annot_adj_pgRNA %>% ungroup() %>% group_by(rep, target_type) %>%
  summarise(m = median( CRISPR_score,na.rm = TRUE))

## Begin the scoring part
d.lfc_pgRNA = d.lfc_annot_adj_pgRNA

#### Expected GI score


#### Expected GI score

##Insert explanation of how to calculate expected GI scores here** 



## confirm that getting mean single-targeting CRISPR scores is acting as expected:
d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "PC" | target_type == "ST") %>%
  mutate(targeting_gRNA_seq = case_when(
    gene2 == "Fluc" ~ gRNA1_seq,
    gene1 == "Fluc" ~ gRNA2_seq
  )) %>%
  group_by(rep, paralog_pair, targeting_gRNA_seq) %>%
  mutate(mean_single_target_CS = mean(CRISPR_score)) %>%
  dplyr::select(pgRNA_id:gRNA2_seq, targeting_gRNA_seq, mean_single_target_CS) %>%
  arrange(rep, paralog_pair, targeting_gRNA_seq)



## calculate mean CRISPR score of single-targeting pgRNAs containing the same targeting
## sgRNA sequence but different control sgRNA sequences
d.mean_single_target_CS <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "PC" | target_type == "ST") %>%
  mutate(targeting_gRNA_seq = case_when(
    gene2 == "Fluc" ~ gRNA1_seq,
    gene1 == "Fluc" ~ gRNA2_seq
  )) %>%
  group_by(rep, paralog_pair, targeting_gRNA_seq) %>%
  mutate(mean_single_target_CS = mean(CRISPR_score, na.rm = TRUE)) %>%
  dplyr::select(rep, targeting_gRNA_seq, mean_single_target_CS) %>%
  distinct(rep, targeting_gRNA_seq,paralog_pair, .keep_all = TRUE) %>% 
  ungroup() %>% 
  select(-c(paralog_pair))



## add mean single-targeting CRISPR scores into double-targeting DF so I can calculate
## expected GI scores by summing 

## join single-target CRISPR scores with double-targeting pgRNA df based on targeting
## sgRNA sequences

## get just double-targeting pgRNAs
d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "DT") %>%
  dplyr::select(pgRNA_id, rep, paralog_pair, CRISPR_score, target_type,
                gRNA1_seq, gRNA2_seq, gene1, gene2)


# removing paralog_pair in join as 1 targetting sey is always for on single pair 
d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_double_targeting %>%
  rename(double_target_CS = CRISPR_score) %>%
  left_join(d.mean_single_target_CS, by = c("rep", "gRNA1_seq" = "targeting_gRNA_seq")) %>%
  rename(mean_gRNA1_single_target_CS = mean_single_target_CS) %>%
  left_join(d.mean_single_target_CS, by = c("rep",  "gRNA2_seq" = "targeting_gRNA_seq")) %>%
  rename(mean_gRNA2_single_target_CS = mean_single_target_CS)

## save a list of pgRNAs that are filtered out based on single targeting CS being NA
d.rm_double_targeting_pgRNAs <- d.lfc_pgRNA_double_targeting %>%
  dplyr::select(pgRNA_id, rep, paralog_pair, double_target_CS, gRNA1_seq, gRNA2_seq, 
                mean_gRNA1_single_target_CS, mean_gRNA2_single_target_CS) %>%
  filter(is.na(mean_gRNA1_single_target_CS) | is.na(mean_gRNA2_single_target_CS))
d.rm_double_targeting_pgRNAs


d.lfc_pgRNA_double_targeting <- d.lfc_pgRNA_double_targeting %>%
  ## filter out pgRNAs where either single-targeting mean CRISPR score is NA
  filter(!is.na(mean_gRNA1_single_target_CS) & !is.na(mean_gRNA2_single_target_CS)) %>%
  ## calculate expected double-targeting GI score by summing the two mean single-targeting
  ## CRISPR scores for that paralog pair
  mutate(expected_CS = mean_gRNA1_single_target_CS + mean_gRNA2_single_target_CS)

### Single-targeting
## get just single-targeting pgRNAs
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "ST" | target_type == "PC") %>%
  dplyr::select(pgRNA_id, rep, paralog_pair, CRISPR_score, target_type,
                gRNA1_seq, gRNA2_seq, gene1, gene2) %>%
  mutate(targeting_gRNA_seq = case_when(
    gene2 == "Fluc" ~ gRNA1_seq,
    gene1 == "Fluc" ~ gRNA2_seq
  )) %>%
  mutate(control_gRNA_seq = case_when(
    gene2 == "Fluc" ~ gRNA2_seq,
    gene1 == "Fluc" ~ gRNA1_seq
  ))

## gRNA at any position. 
#gRNA seq at postion A or B <- take mean
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

## In thompson dataset, there is no "other_control_seq". One pair, A_Ctrl Ctrl is always a single seq . 
# For example, for AARS_ctrl, there are 5 separate seqs for AARS but only one in place of CTRL. 
# so other_single_target_CS  must always be NA
# Therefore this step here is different from others.



d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  mutate(other_single_target_CS = NA)




d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  left_join(y = d.mean_double_control_CS, by = c("rep", "control_gRNA_seq"))


## So we have no mean double control and also no other single target. 



## save a list of pgRNAs that will be removed
d.rm_single_targeting_pgRNAs <- d.lfc_pgRNA_single_targeting %>%
  filter(is.na(other_single_target_CS) | is.na(mean_double_control_CS))
d.rm_single_targeting_pgRNAs

## filter out pgRNAs whose CS = NA ### THis is all 0s here. Deal with this. 
#d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
#  filter(!is.na(other_single_target_CS) & !is.na(mean_double_control_CS)) # 0 rows here

## My way to deal this is to set expected_CS = CRISPR_score !

d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  mutate(expected_CS = CRISPR_score ) %>% #other_single_target_CS + mean_double_control_CS) %>%
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

d.lfc_pgRNA_single_targeting_GI <- d.lfc_pgRNA_single_targeting %>%
  left_join(d.lfc_pgRNA_single_targeting_mean_lm_summary, by = "rep") %>%
  mutate(GI_score = single_target_CS - (intercept + slope * expected_CS)) 

### Calculate GI score for double-targeting pgRNAs

d.lfc_pgRNA_double_targeting_GI <- d.lfc_pgRNA_double_targeting %>%
  left_join(d.lfc_pgRNA_single_targeting_mean_lm_summary, by = "rep") %>%
  mutate(GI_score = double_target_CS - (intercept + slope * expected_CS)) # almost same as observed - expected

#Result of this test is always 1, so use this or just do GI_score  = double_target_CS - expected_CS is just the same 
# as our intercept and slope are 0 and 1 respectively
d.GI_scores_pgRNA_single <- d.lfc_pgRNA_single_targeting_GI %>%
  rename(observed_CS = single_target_CS) %>%
  dplyr::select(-c( mean_double_control_CS,
                    targeting_gRNA_seq, control_gRNA_seq)) %>%
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
  ) %>% 
  separate(col = paralog_pair, into = c("gene1", "gene2"), sep = "_") %>%
  mutate(paralog_pair = ifelse(gene1 < gene2, paste0(gene1,"_", gene2), paste0(gene2,"_", gene1))) %>%
  select(-gene1, -gene2)


write.csv(df_wide, "Parrish Score Scripts/ParrishOutput/Thompson_parrish.csv", row.names = FALSE)
