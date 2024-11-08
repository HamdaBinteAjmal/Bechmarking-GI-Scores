

library(orthrus)
library(stringr) 
library(readxl)

output_folder <- file.path("Orthrus Scripts/OrthrusOutput/")

if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }


df <- data.frame(read_excel(file.path("InputData/Ito/Count_data_ParalogV1.xlsx")))



nonessentials <- c("AAVS1") 
df = df %>% rename(gene1 = Aureus_gene, gene2 =  Pyogenes_gene) %>%
  rename(Left.sgRNA_Right.sgRNA = last_col()) %>%
  mutate(Guide1 =  str_trim(str_split_fixed(Left.sgRNA_Right.sgRNA, "[_]", 2)[,1]), 
         Guide2 =  str_trim(str_split_fixed(Left.sgRNA_Right.sgRNA, "[_]", 2)[,2])) %>%
           select(-c("Left.sgRNA_Right.sgRNA")) %>%
  mutate(gene1 = ifelse(gene1 %in% nonessentials, "NegControl", gene1)) %>%
  mutate(gene2 = ifelse(gene2 %in% nonessentials, "NegControl", gene2))
  
sample_table <- data.frame(
  Screen = c("T0", "Meljuso", "GT1", "MEL202_203", "PK1", "MEWO", "HS944T", "IPC298", "A549", "HSC5", "HS936T", "PATU8988S"),
  Replicates = c("pDNA", "Meljuso_RepA;Meljuso_RepB;Meljuso_RepC", "GI1_004_RepA;GI1_004_RepB;GI1_004_RepC", 
                 "MEL202_003_RepA;MEL202_003_RepB;MEL202_003_RepC", "PK1_REPA;PK1_REPB;PK1_REPC", "MEWO_REPA;MEWO_REPB;MEWO_REPC",
                 "HS944T_REPA;HS944T_REPB;HS944T_REPC", "IPC298_REPA;IPC298_REPB;IPC298_REPC", "A549_REPA;A549_REPB;A549_REPC",
                 "HSC5_REPA;HSC5_REPB;HSC5_REPC", "HS936T_REPA;HS936T_REPB;HS936T_REPC", "PATU8988S_REPA;PATU8988S_REPB"),
  NormalizeTo = c(NA, "T0", "T0", "T0", "T0", "T0", "T0", "T0", "T0", "T0", "T0", "T0")
)

batch_table <- data.frame(
  Screen = c("Meljuso", "GT1", "MEL202_203", "PK1", "MEWO", "HS944T", "IPC298", 
             "A549", "HSC5", "HS936T", "PATU8988S"),
  Control = rep("combn", 11)
)


if (!dir.exists(output_folder)) { dir.create(output_folder) }

# Processes data 
screens <- add_screens_from_table(sample_table)

df <- normalize_screens(df, screens, filter_names = "T0") 

guides <- split_guides(df, screens, "Guide1", "Guide2") 

dual <- guides[["dual"]] 
single <- guides[["single"]] 
combn <- guides[["combn"]] 

temp <- score_combn_vs_single(combn, single, screens,
                              screens_to_score, test = "moderated-t", 
                              return_residuals = FALSE, filter_genes = c("NegControl"),ignore_orientation = TRUE, ) 
paralog_scores <- temp[["scored_data"]] 
write.table(paralog_scores, file.path(output_folder,
                                      "ito_orthrus.tsv"), sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE) 


score_combn_batch(combn, single, screens, batch_table, output_folder, test = "moderated-t", loess = FALSE, ignore_orientation = TRUE, filter_genes = c("NegControl"), neg_type = "Sensitizer", pos_type = "Suppressor", fdr_threshold = 0.2, differential_threshold = 0.5) 
