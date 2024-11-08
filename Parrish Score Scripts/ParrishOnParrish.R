library(ggplot2)
library(RColorBrewer)
library(corrr)
library(stringr)
library(tidyr)
library(dplyr)

## Note: This is my own implementation of the Parrish score. Most parts of the code is taken from the code provided by Parrish et al. 

parrish = read.csv("InputData/Parrish/GSE178179_pgPEN_counts_PC9.txt", sep = "\t")
PC9_columns <- c("PC9_ETP_RepA", "PC9_ETP_RepB", "PC9_ETP_RepC", 
                   "PC9_LTP_RepA", "PC9_LTP_RepB", "PC9_LTP_RepC")
HeLa_columns <-  c("HeLa_ETP", 
                   "HeLa_LTP_RepA", "HeLa_LTP_RepB", "HeLa_LTP_RepC")

essentials = read.csv("InputData/Parrish/AchillesCommonEssentialControls.csv")
essentials <- essentials %>%
  separate(Gene, into = c("gene_symbol", "entrez_id"), sep = "\\s*\\(", convert = TRUE) %>%
  mutate(entrez_id = str_replace(entrez_id, "\\)", ""))  %>% select(gene_symbol) %>% unlist() %>% unname()

## apply pre processing as guided by the Parrish et al. paper
parrish = parrish %>%
  separate(paralog_pair, into = c("gene1", "gene2"), sep = "\\|", remove = FALSE, extra = "merge") %>%
  mutate(gene1  = if_else(str_starts(gene1, "NTpg" ), "NT", gene1)) %>%
  mutate(gene2  = if_else(str_starts(gene2, "NTpg" ), "NT", gene2)) %>%
  mutate(gene1 = str_replace(gene1, "^nt[0-9]$", "NT")) %>%
  mutate(gene2 = str_replace(gene2, "^nt[0-9]$", "NT")) %>%
  mutate(gene1 = ifelse(gene1 == "NA", "NT", gene1)) %>%
  mutate(gene2 = ifelse(gene2 == "NA", "NT", gene2)) %>%
  ## RPM is reads/all reads in the sample * 10^6
  ## pgRNAs with < 2 reads per million (RPM) in the plasmid pool or with a read count of zero at any time point were also removed
  
  mutate(PC9_RPM = PC9_plasmid/sum(PC9_plasmid)*10^6, HeLa_RPM = HeLa_plasmid/sum(HeLa_plasmid)*10^6) %>%
  mutate(across(starts_with("PC9"), ~ifelse(PC9_RPM < 2, NA_real_,.)  )) %>%
  mutate(across(starts_with("HeLa"), ~ifelse(HeLa_RPM < 2, NA_real_,.)  )) %>%

  
  mutate(PC9_0_RepA = PC9_LTP_RepA ==0 | PC9_ETP_RepA == 0,
         PC9_0_RepB = PC9_LTP_RepB ==0 | PC9_ETP_RepB == 0,
         PC9_0_RepC = PC9_LTP_RepC ==0 | PC9_ETP_RepC == 0,
         HeLa_0_RepA = HeLa_ETP ==0 | HeLa_LTP_RepA == 0,
         HeLa_0_RepB = HeLa_ETP ==0 | HeLa_LTP_RepB == 0,
         HeLa_0_RepC = HeLa_ETP ==0 | HeLa_LTP_RepC == 0,
         
         
         
         ) %>% # Wrong, we need to do it by replicates
  
  mutate(
    PC9_LTP_RepA = if_else(PC9_0_RepA,  NA_real_, PC9_LTP_RepA),
    PC9_LTP_RepB = if_else(PC9_0_RepB,  NA_real_, PC9_LTP_RepB),
    PC9_LTP_RepC = if_else(PC9_0_RepC,  NA_real_, PC9_LTP_RepC),
    HeLa_LTP_RepA = if_else(HeLa_0_RepA,  NA_real_, HeLa_LTP_RepA),
    HeLa_LTP_RepB = if_else(HeLa_0_RepB,  NA_real_, HeLa_LTP_RepB),
    HeLa_LTP_RepC = if_else(HeLa_0_RepC,  NA_real_, HeLa_LTP_RepC) )%>%
  
  
  
  select(-c("HeLa_ETP", starts_with("PC9_ETP"), starts_with("PC9_0"), starts_with("HeLa_0"))) %>%
  mutate(PC9_LTP_RepA = log2(PC9_LTP_RepA / PC9_plasmid),
         PC9_LTP_RepB = log2(PC9_LTP_RepB / PC9_plasmid),
         PC9_LTP_RepC = log2(PC9_LTP_RepC / PC9_plasmid),
         HeLa_LTP_RepA = log2(HeLa_LTP_RepA / HeLa_plasmid),
         HeLa_LTP_RepB = log2(HeLa_LTP_RepB / HeLa_plasmid),
         HeLa_LTP_RepC = log2(HeLa_LTP_RepC / HeLa_plasmid)) %>% 
  mutate(PC9= rowMeans(select(., PC9_LTP_RepA, PC9_LTP_RepB, PC9_LTP_RepC), na.rm = TRUE)) %>%
  mutate(HeLa= rowMeans(select(., HeLa_LTP_RepA, HeLa_LTP_RepB, HeLa_LTP_RepC), na.rm = TRUE)) %>%
  
  mutate(gene1_is_essential = gene1 %in% essentials,
         gene2_is_essential = gene2 %in% essentials) %>%
  mutate(target_type = NA) %>% 
  mutate(target_type = case_when(
    gene1 == "NT" & gene2 == "NT" ~ "NC",
    ((gene1 == "NT") != (gene2 == "NT")) & (gene2_is_essential | gene1_is_essential) ~ "PC",
    (gene1 != "NT" & gene2 !="NT") & (gene1 != gene2) ~ "DT",
    TRUE ~ "ST")) %>%
  mutate(paralog_pair = paste0(gene1, "_", gene2)) %>%
  select(c(pgRNA_id,gRNA1_seq,gRNA2_seq,gene1, gene2,target_type,paralog_pair,
           PC9_LTP_RepA, PC9_LTP_RepB,PC9_LTP_RepC,HeLa_LTP_RepA, HeLa_LTP_RepB, HeLa_LTP_RepC, PC9, HeLa, gene2_is_essential, gene1_is_essential  ))

  
  




parrish_long = parrish %>%
  pivot_longer(
    cols = starts_with("PC9") | starts_with("HeLa"), # Select columns that start with PC9 or HeLa
    names_to = "rep",  # Column where the old column names will be stored
    values_to = "LFC"   # Column where the values will be stored
  ) %>%
  filter(!is.na(LFC))
  
  
  
  # Remove the total counts column after calculation
d.reps <- parrish_long %>%
  ungroup() %>%
  dplyr::select(rep) %>%
  distinct()

d.lfc_rep_cor <- parrish_long %>%
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


d.control_group_medians <- parrish_long %>%
  group_by(rep, target_type ) %>%
  filter(target_type == "NC" | target_type == "PC") %>%
  summarize(median_lfc = median(LFC, na.rm = TRUE))
d.control_group_medians

df_LFC <- parrish_long %>%
  group_by(rep) %>%
  mutate(lfc_adj1 = LFC - median(LFC[target_type == "NC"], na.rm = TRUE),
         lfc_adj2 = lfc_adj1 / (median(lfc_adj1[target_type == "NC"], na.rm = TRUE) -
                                  median(lfc_adj1[target_type == "PC"], na.rm = TRUE)))%>%
  ungroup()


df_LFC %>%
  filter(str_detect(rep, "PC9")) %>%

ggplot( aes(x = target_type, y = LFC , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "PC9") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

df_LFC %>%
  filter(str_detect(rep, "PC9")) %>%
  
  ggplot( aes(x = target_type, y = lfc_adj2 , fill = target_type)) +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, coef = 0, width = 0.1) +
  labs(x = "pgRNA_category", y = "PC9") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
##########################
mapping = read.delim("InputData/genenames.txt", header = TRUE)

mapping = mapping %>% dplyr::select(c("Approved.symbol", "Ensembl.gene.ID", "Previous.symbols")) %>%
  separate_rows(Previous.symbols, sep = ",") %>% 
  mutate(Previous.symbols = trimws(Previous.symbols)) %>% 
  dplyr::select(c(Ensembl.gene.ID, Approved.symbol,Previous.symbols )) %>%
  filter(Ensembl.gene.ID != "")





parrish_genes = unique(c(df_LFC$gene1, df_LFC$gene2))
parrish_genes <- data.frame(gene  = parrish_genes)


parrish_genes = parrish_genes %>% 
  filter(gene != "NT") %>%
  left_join(mapping, by = c("gene"= "Approved.symbol")) %>%
  dplyr::select(c("gene" , "Ensembl.gene.ID" )) %>%
  distinct() %>%
  
  group_by(gene) %>%
  mutate(is_duplicate = n() > 1) %>%
  filter(!(is_duplicate == TRUE & Ensembl.gene.ID == "")) %>%
  ungroup()

df_without_na <- parrish_genes %>%
  filter(!is.na(Ensembl.gene.ID))


df_filled <- parrish_genes %>%
  filter(is.na(Ensembl.gene.ID)) %>%
  left_join(mapping, by = c("gene" = "Previous.symbols")) %>%
  mutate(Ensembl.gene.ID = coalesce(Ensembl.gene.ID.x, Ensembl.gene.ID.y)) %>%
  dplyr::select(c(-Ensembl.gene.ID.y, -Ensembl.gene.ID.x)) # Remove the duplicate column from the mapping dataframe

parrish_genes <- bind_rows(df_without_na, df_filled) %>% 
  dplyr::select(c(gene, Ensembl.gene.ID))

## This dataset has been downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120703


TPM = as_tibble(read.table("InputData/Parrish/GSE120703_d.expr.txt",header = TRUE))
TPM = TPM %>%
  mutate(is_expressed_PC9 = TPM_PC9 >= 2, is_expressed_HeLa = TPM_PC9 >=  2)


df_LFC = df_LFC %>%
  mutate(sample = ifelse(str_detect(rep, "PC9"), "PC9", "HeLa"))


parrish_genes = parrish_genes %>% left_join(TPM, by = c("Ensembl.gene.ID" = "gene")) 
PC9_isExp = parrish_genes %>% select(c(gene,is_expressed_PC9 ))
HeLa_isExp = parrish_genes %>% select(c(gene,is_expressed_HeLa ))

df_LFC <- df_LFC %>% 
  mutate(gene1_expressed_flag = NA, gene2_expressed_flag = NA) %>%
  left_join(PC9_isExp, by= c("gene1" = "gene")) %>%
  mutate(gene1_expressed_flag = ifelse(sample == "PC9", is_expressed_PC9, gene1_expressed_flag)) %>%
  select(-c(is_expressed_PC9) ) %>%
  left_join(PC9_isExp, by= c("gene2" = "gene")) %>%
  mutate(gene2_expressed_flag = ifelse(sample == "PC9", is_expressed_PC9, gene2_expressed_flag)) %>%
  select(-c(is_expressed_PC9) ) %>%
  left_join(HeLa_isExp, by= c("gene1" = "gene")) %>%
  mutate(gene1_expressed_flag = ifelse(sample == "HeLa", is_expressed_HeLa, gene1_expressed_flag)) %>%
  select(-c(is_expressed_HeLa) ) %>%
  left_join(HeLa_isExp, by= c("gene2" = "gene")) %>%
  mutate(gene2_expressed_flag = ifelse(sample == "HeLa", is_expressed_HeLa, gene2_expressed_flag)) %>%
  select(-c(is_expressed_HeLa )) %>%
  mutate(gene1_expressed_flag = ifelse(is.na(gene1_expressed_flag) & gene1 != "NT", TRUE, gene1_expressed_flag))%>%
  mutate(gene2_expressed_flag = ifelse(is.na(gene2_expressed_flag) & gene2 != "NT", TRUE, gene2_expressed_flag)) # we assuming that if NA, means its expressed.


####################################
### Expression
# Text from the paper:
#"Since the pgPEN library uses non-targeting controls, 
#we adjusted for the fact that single-targeting pgRNAs generate only
#two double-strand breaks (1 per allele), whereas the double-targeting pgRNAs 
#generate four DSBs. To do this, we set the median (adjusted) LFC for 
#unexpressed genes of each group to zero."

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
  mutate(mm = median(lfc_adj2[unexpressed_ctrl_flag == TRUE], na.rm = TRUE)) %>%
  mutate(lfc_adj3 = lfc_adj2 - median(lfc_adj2[unexpressed_ctrl_flag == TRUE], na.rm = TRUE))

d.lfc_annot_adj_single_summary <- d.lfc_annot_adj_single %>%
  group_by(rep, unexpressed_ctrl_flag) %>%
  summarize(median = median(lfc_adj3, na.rm = TRUE))
d.lfc_annot_adj_single_summary


#### Double-targeting


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
  select(-c("lfc_adj2", "gene1_expressed_flag", "gene2_expressed_flag" ,"gene1_is_essential", "gene2_is_essential"))


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

## confirm that getting mean single-targeting CRISPR scores is acting as expected:
d.lfc_pgRNA %>%
  ungroup() %>%
  filter(target_type == "PC" | target_type == "ST") %>%
  mutate(targeting_gRNA_seq = case_when(
    gene2 == "NT" ~ gRNA1_seq,
    gene1 == "NT" ~ gRNA2_seq
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
    gene2 == "NT" ~ gRNA1_seq,
    gene1 == "NT" ~ gRNA2_seq
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

## filter out single-targeting CS == NA rows
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
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  left_join(d.other_single_targeting_CS, by = c("rep", "paralog_pair", "targeting_gRNA_seq")) %>%
  mutate(same_control_seq = ifelse(control_gRNA_seq == other_control_seq, TRUE, FALSE)) %>%
  filter(same_control_seq == FALSE) %>%
  dplyr::select(-same_control_seq)

d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  left_join(y = d.mean_double_control_CS, by = c("rep", "control_gRNA_seq"))




## save a list of pgRNAs that will be removed
d.rm_single_targeting_pgRNAs <- d.lfc_pgRNA_single_targeting %>%
  filter(is.na(other_single_target_CS) | is.na(mean_double_control_CS))
d.rm_single_targeting_pgRNAs

## filter out pgRNAs whose CS = NA
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  filter(!is.na(other_single_target_CS) & !is.na(mean_double_control_CS))
d.lfc_pgRNA_single_targeting <- d.lfc_pgRNA_single_targeting %>%
  mutate(expected_CS = other_single_target_CS + mean_double_control_CS) %>%
  rename(single_target_CS = CRISPR_score)



#### Linear model



d.lfc_pgRNA_single_targeting

d.lfc_pgRNA_single_targeting_mean <- d.lfc_pgRNA_single_targeting %>%
  mutate(pgRNA_target = ifelse(gene1 == "NT", gene2, gene1)) %>% ## i think so :D
  group_by(rep, pgRNA_target) %>%
  summarize(mean_expected_CS = mean(expected_CS),
            mean_observed_CS = mean(single_target_CS)) %>%
  ungroup()
d.lfc_pgRNA_single_targeting_mean

## fit linear model to target-level mean single-targeting pgRNA expected vs. 
## observed values and extract slope and intercept values
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


ggplot(data = d.lfc_pgRNA_single_targeting, 
       aes(x = expected_CS, y = single_target_CS)) +
  geom_point() + 
  geom_density2d(aes(color = ..level..), alpha = 0.9, contour_var = "count") +
  geom_abline(data = d.lfc_pgRNA_single_targeting_mean_lm_summary,
              aes(slope = slope, intercept = intercept), 
              color = "gray55", size = 1) +
  scale_colour_gradientn(colors = rev(contour_palette)) +
  labs(x = "expected_CRISPR_score", 
       y = "observed_CRISPR_score",
       color = "count") +
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
# results


d.stats <- d.GI_scores_pgRNA_double %>%
  dplyr::select(paralog_pair, p_val, fdr) %>%
  distinct(paralog_pair, .keep_all = TRUE) 

## add p-val and fdr to single-targeting
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
## why only 3,086 targets not 3,090? some had too many single-targeting pgRNAs filtered out

ggplot(d.GI_scores_target, aes(x = mean_expected_CS, y = mean_observed_CS, color = broad_target_type)) +
  geom_point() + 
  geom_abline(data = d.lfc_pgRNA_single_targeting_mean_lm_summary,
              aes(slope = slope, intercept = intercept),
              color = "gray55", size = 1) +
  labs(x = "expected_CRISPR_score", y = "observed_CRISPR_score") +

  facet_wrap(~rep)



score = d.GI_scores_target %>%
  filter(broad_target_type == "double_targeting") %>%
  filter(rep == "PC9" | rep == "HeLa") %>%
  group_by(paralog_pair, rep) %>%
  mutate(mean_final_score = mean(mean_GI_score),
         mean_fdr = mean(fdr),
         .groups = "drop") %>%
  select(c(paralog_pair, rep, mean_final_score , mean_fdr)) %>%
  distinct() %>%
  pivot_wider(
    names_from = rep,   # Column that turns into new column names
    values_from = c(mean_final_score, mean_fdr) ,   # Column that populates values in the new columns
    names_sep = "_"
  )  



score = score %>% 
  rename(PC9 = mean_final_score_PC9,
         HeLa = mean_final_score_HeLa,
         PC9_fdr = mean_fdr_PC9,
         HeLa_fdr = mean_fdr_HeLa)
write.csv(x = score, file =  "Parrish Score Scripts/ParrishOutput/Parrish_Parrish.csv", row.names = FALSE)

test = read.csv(file= "Parrish Score Scripts/ParrishOutput/Parrish.csv")
