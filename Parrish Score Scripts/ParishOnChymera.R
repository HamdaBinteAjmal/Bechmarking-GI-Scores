##libraries target more than 450 core essential genes33 with over
#6,000 Cas9 and Cas12a exon-targeting guides and over 35,000 exonflanking
#guides. They also contain 1,000 control constructs targeting
#intergenic regions, which control for toxicity induced by doublestranded
#(ds)DNA breaks

## 9986 unique seqs for NegControl_Ne

##For “dual-targeting” guides which target the same gene twice,
# the gene1 column must contain the name of the targeted gene and the gene2 column must contain the string None.
# Perhaps remove anything that has "none" in gene2 as its targetting same gene twice

## For single-targeting guides paired with a standardized negative control, such as an intergenic region for the dataset analyzed in Procedure 2 or a non-essential gene for the dataset analyzed in Procedure 3 (dede),
# the control must be named NegControl in the corresponding gene label column.# Dede non-essential is NegControl 

# If gene1 is NegControl, Cas5.Guide.Type is always intergenic
# if gene2 is NegControl, Cpf1.Guide.Type is always intergenic


# From Chymera doc: As for dual-targeted scoring, 
# additionally pass in “NT” to the filter_genes parameter to ignore non-targeting control genes during scoring..
# NT is non targeting
library(orthrus)
library(readxl)
library(tibble)

#Take Chymera data directly from the R library

chymera <- chymera_paralog 
rownames(chymera) = paste0(chymera$Cas9.Guide, ";", chymera$Cpf1.Guide )
chymera = chymera %>%
  filter(gene1 != "NT" & gene2 != "NT")

# Anything with "None" means its targetting the same gene twice , dont need it for now.
chymera = chymera %>%
  filter (gene1 != "None" & gene2 != "None")

# positive control, this screen includes 575 genes from the pool of Core Essential Genes (CEG2) (Hart et al., 2017). 
pc_genes = read_xlsx("InputData/Chymera/chymera_essential.xlsx")
pc_genes = pc_genes %>% pull(1)


chymera = chymera %>% 
  rename(pgRNA_id = ID) %>%
  select(-c(Library, Cas9.Target.Site, Cas9.Guide.Source,Cas9.Guide.Type,
            Cpf1.Target.Site, Cpf1.Guide.Type, CNN.Score,contains("Torin"), )) %>%
  rename(gRNA1_seq = Cas9.Guide,gRNA2_seq = Cpf1.Guide) %>%
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
    gene1 == "NegControl" & gene2 == "NegControl" ~ "NC",
    ((gene1 == "NegControl") != (gene2 == "NegControl")) & (gene2_is_essential | gene1_is_essential) ~ "PC",
    (gene1 != "NegControl" & gene2 !="NegControl") & (gene1 != gene2) ~ "DT",
    TRUE ~ "ST"))  %>% 
  #discarded. pgRNAs with < 2 reads per million (RPM) in the plasmid pool or with a read count of zero at any time point were also removed.
  
  mutate(HAP1_RPM = HAP1.T0/sum(HAP1.T0)*10^6) %>%
  mutate(RPE1_RPM = RPE1.T0/sum(RPE1.T0)*10^6) %>%
  mutate(across(starts_with("RPE1"), ~ifelse(RPE1_RPM < 2, NA_real_,.))) %>% 
  mutate(across(starts_with("HAP1"), ~ifelse(HAP1_RPM < 2, NA_real_,.))) %>% 
  select(-c(contains("RPM"))) %>%
  
  mutate(RPE1_0_A = RPE1.T18A ==0 | RPE1.T24A == 0,
         RPE1_0_B = RPE1.T18B ==0 | RPE1.T24B == 0,
         RPE1_0_C = RPE1.T18C ==0 | RPE1.T24C == 0,
         
         HAP1_0_A = HAP1.T12A ==0 | HAP1.T18A == 0,
         HAP1_0_B = HAP1.T12B ==0 | HAP1.T18B == 0,
         HAP1_0_C = HAP1.T12C ==0 | HAP1.T18C == 0) %>%
  
  mutate(
    RPE1.T24A = if_else(RPE1_0_A | is.na(RPE1_0_A),  NA_real_, RPE1.T24A),
    RPE1.T24B = if_else(RPE1_0_B | is.na(RPE1_0_B),  NA_real_, RPE1.T24B),
    RPE1.T24C = if_else(RPE1_0_C | is.na(RPE1_0_C),  NA_real_, RPE1.T24C),
    
    HAP1.T18A = if_else(HAP1_0_A| is.na(HAP1_0_A),  NA_real_, HAP1.T18A),
    HAP1.T18B = if_else(HAP1_0_B| is.na(HAP1_0_B),  NA_real_, HAP1.T18B),
    HAP1.T18C = if_else(HAP1_0_C| is.na(HAP1_0_C),  NA_real_, HAP1.T18C)) %>%
  
  
  select(-c(contains("HAP1.T12"), contains("_0_"), contains("RPE1.T18")))%>%
  

  mutate(RPE1.T24A = log2(RPE1.T24A / RPE1.T0),
         RPE1.T24B = log2(RPE1.T24B / RPE1.T0),
         RPE1.T24C = log2(RPE1.T24C / RPE1.T0),
         HAP1.T18A = log2(HAP1.T18A / HAP1.T0),
         HAP1.T18B = log2(HAP1.T18B / HAP1.T0),
         HAP1.T18C = log2(HAP1.T18C / HAP1.T0)  )%>% 
  mutate(RPE1= rowMeans(select(., RPE1.T24A, RPE1.T24B, RPE1.T24C), na.rm = TRUE)) %>%
  mutate(HAP1= rowMeans(select(., HAP1.T18A, HAP1.T18B, HAP1.T18C), na.rm = TRUE)) %>%
  select(-c(contains("HAP1.T18"), contains("RPE1.T24"))) 


chymera_long = chymera %>%
  pivot_longer(
    cols = c("RPE1", "HAP1"), # Select columns that start with PC9 or HeLa
    names_to = "rep",  # Column where the old column names will be stored
    values_to = "LFC"   # Column where the values will be stored
  ) %>%
  filter(!is.na(LFC))
d.reps <- chymera_long %>%
  ungroup() %>%
  dplyr::select(rep) %>%
  distinct()

d.lfc_rep_cor <- chymera_long %>%
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


d.control_group_medians <- chymera_long %>%
  group_by(rep, target_type ) %>%
  filter(target_type == "NC" | target_type == "PC") %>%
  summarize(median_lfc = median(LFC, na.rm = TRUE))
d.control_group_medians


df_LFC <- chymera_long %>%
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
  filter(str_detect(rep, "HAP1")) %>%
  
  ggplot( aes(x = target_type, y = LFC , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "PC9") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

df_LFC %>%
  filter(str_detect(rep, "HAP1")) %>%
  
  ggplot( aes(x = target_type, y = lfc_adj2 , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "PC9") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# This file is downloaded from DepMap "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
TPMs = read.csv("InputData/TPM.csv")

TPMs = TPMs %>%
  
  column_to_rownames(var = "X") %>%  # 'X' is the name of the new column
  t() %>%
  as.data.frame() %>%
  rownames_to_column("GENE") %>%
  mutate(GENE= sub("\\.{2}.*$", "", GENE)) %>% 
  rename(HAP1  = `ACH-002475` ) %>%
  rename(RPE1 = `ACH-002464`) %>%
  select(c(GENE, HAP1, RPE1)) %>%
  mutate(across(where(is.numeric), ~ .x >= log2(2+1))) %>% 
  pivot_longer(
    cols = starts_with("HAP1") | starts_with("RPE1"), # Select columns that start with PC9 or HeLa
    names_to = "rep",  # Column where the old column names will be stored
    values_to = "is_expressed"   # Column where the values will be stored
  ) #%>%

df_LFC  = df_LFC %>%
  left_join(TPMs, by = c("gene1" = "GENE", "rep" = "rep")) %>%
  rename(gene1_expressed_flag = is_expressed) %>%
  left_join(TPMs, by = c("gene2" = "GENE", "rep" = "rep")) %>%
  rename(gene2_expressed_flag = is_expressed)

df_LFC = df_LFC %>%
  mutate(gene1_expressed_flag = ifelse(is.na(gene1_expressed_flag) & gene1 != "NegControl", TRUE, gene1_expressed_flag))%>%
  mutate(gene2_expressed_flag = ifelse(is.na(gene2_expressed_flag) & gene2 != "NegControl", TRUE, gene2_expressed_flag)) # we assuming that if NA, means its expressed.

d.lfc_annot_adj_single <- df_LFC %>% 
  filter(target_type == "PC" | target_type == "ST") %>%
  ## make a flag variable to indicate which pgRNAs are targeting unexpressed
  ## single targets
  mutate(unexpressed_ctrl_flag = case_when(
    gene2 == "NegControl" & gene1_expressed_flag == FALSE ~ TRUE, # for NAs, we assuming its expressed
    gene1 == "NegControl" & gene2_expressed_flag == FALSE ~ TRUE,
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



d.lfc_annot_adj_single <- d.lfc_annot_adj_single %>%
  dplyr::select(-unexpressed_ctrl_flag)

### double targeting
d.lfc_annot_adj_double <- d.lfc_annot_adj_double %>%
  dplyr::select(-unexpressed_ctrl_flag)

## bind rows
d.lfc_annot_adj_pgRNA <- bind_rows(d.lfc_annot_adj_double, d.lfc_annot_adj_single, d.lfc_annot_adj_control)


d.lfc_annot_adj_pgRNA <- d.lfc_annot_adj_pgRNA %>%
  ungroup() %>%
  rename(CRISPR_score = lfc_adj3) %>%
  select(-c("lfc_adj2", "LFC", "lfc_adj1" ,"gene1_expressed_flag", "gene2_expressed_flag" ,"gene1_is_essential", "gene2_is_essential"))



d.lfc_annot_adj_pgRNA %>%
  filter(str_detect(rep, "HAP1")) %>%
  
  ggplot( aes(x = target_type, y = CRISPR_score , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "PC9") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


d.lfc_annot_adj_pgRNA %>% ungroup() %>% group_by(rep, target_type) %>%
  summarise(m = median( CRISPR_score,na.rm = TRUE))


d.lfc_pgRNA = d.lfc_annot_adj_pgRNA


d.lfc_pgRNA = d.lfc_pgRNA %>%
  mutate(paralog_pair = paste0(gene1,"_", gene2))


d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "PC" | target_type == "ST") %>%
  mutate(targeting_gRNA_seq = case_when(
    gene2 == "NegControl" ~ gRNA1_seq,
    gene1 == "NegControl" ~ gRNA2_seq
  )) %>%
  group_by(rep, paralog_pair, targeting_gRNA_seq) %>%
  mutate(mean_single_target_CS = mean(CRISPR_score)) %>%
  dplyr::select(pgRNA_id:gRNA2_seq, targeting_gRNA_seq, mean_single_target_CS) %>%
  arrange(rep, paralog_pair, targeting_gRNA_seq)




d.mean_single_target_CS <- d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "PC" | target_type == "ST") %>%
  mutate(targeting_gRNA_seq = case_when(
    gene2 == "NegControl" ~ gRNA1_seq,
    gene1 == "NegControl" ~ gRNA2_seq
  )) %>%
  group_by(rep, paralog_pair, targeting_gRNA_seq) %>%
  mutate(mean_single_target_CS = mean(CRISPR_score, na.rm = TRUE)) %>%
  dplyr::select(rep, targeting_gRNA_seq, mean_single_target_CS) %>%
  distinct(rep, targeting_gRNA_seq,paralog_pair, .keep_all = TRUE) %>% 
  ungroup() %>% 
  select(-c(paralog_pair))

## get just double-targeting pgRNAs
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
    gene2 == "NegControl" ~ gRNA1_seq,
    gene1 == "NegControl" ~ gRNA2_seq
  )) %>%
  mutate(control_gRNA_seq = case_when(
    gene2 == "NegControl" ~ gRNA2_seq,
    gene1 == "NegControl" ~ gRNA1_seq
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


d.other_single_targeting_CS <- d.lfc_pgRNA_single_targeting %>%
  dplyr::select(rep, paralog_pair, targeting_gRNA_seq, control_gRNA_seq, CRISPR_score) %>%
  rename(other_single_target_CS = CRISPR_score, other_control_seq = control_gRNA_seq)

## At this step, i m losing a large number of data points.
## 

d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  left_join(d.other_single_targeting_CS, by = c("rep", "paralog_pair", "targeting_gRNA_seq")) %>%
  mutate(same_control_seq = ifelse(control_gRNA_seq == other_control_seq, TRUE, FALSE)) %>%
  filter(same_control_seq == FALSE) %>%
  dplyr::select(-same_control_seq)


d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  left_join(y = d.mean_double_control_CS, by = c("rep", "control_gRNA_seq"))



d.rm_single_targeting_pgRNAs <- d.lfc_pgRNA_single_targeting %>%
  filter(is.na(other_single_target_CS) | is.na(mean_double_control_CS))
d.rm_single_targeting_pgRNAs




## This has a lot of NAs for double control
# filter: removed 4,452 rows (8%), 49,842 rows remaining. 94% NAs for double control
# Again, if 94% are NAs for double contro, then I would rather set expected = observed again #
# Otherwise too much data lost.

#d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
#  filter(!is.na(other_single_target_CS) & !is.na(mean_double_control_CS))



d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  mutate(expected_CS = CRISPR_score) %>% #other_single_target_CS + mean_double_control_CS) %>%
  rename(single_target_CS = CRISPR_score)


d.lfc_pgRNA_single_targeting

d.lfc_pgRNA_single_targeting_mean <- d.lfc_pgRNA_single_targeting %>%
  mutate(pgRNA_target = ifelse(gene1 == "NegControl", gene2, gene1)) %>% ## i think so :D
  group_by(rep, pgRNA_target) %>%
  summarize(mean_expected_CS = mean(expected_CS),
            mean_observed_CS = mean(single_target_CS)) %>%
  ungroup()
d.lfc_pgRNA_single_targeting_mean ## only left with 120


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






### Single-target GI scores

d.lfc_pgRNA_single_targeting_GI <- d.lfc_pgRNA_single_targeting %>%
  left_join(d.lfc_pgRNA_single_targeting_mean_lm_summary, by = "rep") %>%
  mutate(GI_score = single_target_CS - (intercept + slope * expected_CS)) 

### Calculate GI score for double-targeting pgRNAs

d.lfc_pgRNA_double_targeting_GI <- d.lfc_pgRNA_double_targeting %>%
  left_join(d.lfc_pgRNA_single_targeting_mean_lm_summary, by = "rep") %>%
  mutate(GI_score = double_target_CS - (intercept + slope * expected_CS))

## reformat to match
d.GI_scores_pgRNA_double <- d.lfc_pgRNA_double_targeting_GI %>%
  rename(observed_CS = double_target_CS) %>%
  dplyr::select(-c(mean_gRNA1_single_target_CS, mean_gRNA2_single_target_CS)) %>%
  mutate(broad_target_type = "double_targeting") %>%
  ungroup()

d.GI_scores_pgRNA_single <- d.lfc_pgRNA_single_targeting_GI %>%
  rename(observed_CS = single_target_CS) %>%
  dplyr::select(-c(other_single_target_CS, mean_double_control_CS,
                   targeting_gRNA_seq, control_gRNA_seq, other_control_seq)) %>%
  mutate(broad_target_type = "single_targeting") %>%
  ungroup()
### Calculate p-values


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
    mutate(p_val = if(length(GI_score) >= 2) {
      t.test(x = single_GI_scores, y = GI_score, paired = FALSE)$p.value
    } else {
      NA_real_  # Return NA of type double
    })
  
  
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
# results
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



final_df = d.GI_scores_target %>% filter(broad_target_type == "double_targeting") %>%
  select(c(rep, paralog_pair, mean_GI_score,fdr, p_val)) %>%
  rename(GI_score = mean_GI_score) %>%
  pivot_wider(
    names_from = rep,
    values_from = c(GI_score, fdr, p_val),
    names_sep = "_"  # This adds an underscore between the rep value and the measure type for clarity
  )
write.csv(final_df, row.names = FALSE, file = "Parrish Score Scripts/ParrishOutput/Chymera_Parrish.csv")
