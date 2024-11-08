## This dataset is different from others in a way that is gene1 = A2M  and gene2 = NegControl, then there are no rows with the reverse:
## that is, gene1 = NegControl and gene2 = A2M. Therefore a hack was appplied. 
## This is also different from Thompson et al data where a negative control is always on position A.

library(dplyr)
library(orthrus)
library(stringr)
library(tibble)
library(tidyr)


output_folder <- file.path("Orthrus Scripts/OrthrusOutput/")

if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }



df = read.csv("InputData/Parrish/GSE178179_pgPEN_counts_PC9.txt", sep = "\t")


df = df %>%
  separate(paralog_pair, into = c("gene1", "gene2"), sep = "\\|", remove = FALSE, extra = "merge") %>%
  mutate( gene1 = ifelse(grepl("^nt[1-8]$",gene1), "NegControl", gene1)) %>%
  mutate( gene2 = ifelse(grepl("^nt[1-8]$",gene2), "NegControl", gene2)) %>%
  mutate(gene1 = ifelse(grepl("^NTpg\\d+$", gene1), "NegControl", gene1)) %>%
  mutate(gene2 = ifelse(gene2 == "NA", "NegControl", gene2 )) %>%
  mutate(rownames = paste0(gRNA1_seq,":",gRNA2_seq)) %>%
  column_to_rownames("rownames")

df = df %>%   select( -c(PC9_plasmid, HeLa_plasmid)) 


df = df %>%
  select(-c(pgRNA_id)) %>%
  rename(GENE = paralog_pair) %>%
  rename(Guide1 = gRNA1_seq) %>%
  rename(Guide2 = gRNA2_seq ) 


df = df %>% 
  select(GENE, PC9_LTP_RepA, PC9_LTP_RepB, PC9_LTP_RepC ,PC9_ETP_RepA, PC9_ETP_RepB, PC9_ETP_RepC,
         HeLa_LTP_RepA ,HeLa_LTP_RepB ,HeLa_LTP_RepC, HeLa_ETP, 
         Guide1, Guide2, gene1, gene2 )

df2 = df %>% mutate(gene3 = gene1, gene4 = gene2) %>% select(-c(gene1,gene2)) %>% rename(gene1 = gene4, gene2 = gene3) %>% mutate(Guide3 = Guide1, Guide4 = Guide2) %>% select(-c(Guide1,Guide2)) %>% rename(Guide1 = Guide4, Guide2 = Guide3)
df2 = df2 %>% rbind(df)



sample_table <- data.frame(
  Screen = c("T0_PC9", 
             "T0_HeLa", 
             "PC9", "HeLa"),
  Replicates = c("PC9_ETP_RepA;PC9_ETP_RepB;PC9_ETP_RepC", 
                 "HeLa_ETP", 
                 "PC9_LTP_RepA;PC9_LTP_RepB;PC9_LTP_RepC", 
                 "HeLa_LTP_RepA;HeLa_LTP_RepB;HeLa_LTP_RepC"),
  NormalizeTo = c(NA, NA, "T0_PC9", "T0_HeLa")
)

batch_table = data.frame(
  Screen = c("PC9", "HeLa"),
  Control = c("combn", "combn")
  )


output_folder <- file.path("Orthrus Scripts/OrthrusOutput/") 

if (!dir.exists(output_folder)) { dir.create(output_folder) }


screens <- add_screens_from_table(sample_table)

df2 <- normalize_screens(df2, screens, min_reads = -Inf ,  max_reads = Inf) 

screens_to_score <- c("PC9", "HeLa") 

guides <- split_guides(df2, screens, "Guide1", "Guide2") 

dual <- guides[["dual"]] 
single <- guides[["single"]] 
combn <- guides[["combn"]] 

score_combn_batch(combn, single, screens, batch_file, output_folder, test = "moderated-t", loess = FALSE, neg_type = "Sensitizer", pos_type = "Suppressor", fdr_threshold = 0.2, differential_threshold = 0.5)


temp <- score_combn_vs_single(combn, single, screens,
                              screens_to_score, test = "moderated-t", 
                              return_residuals = TRUE) 
paralog_scores <- temp[["scored_data"]] 
write.table(paralog_scores, file.path(output_folder,
                                      "parrish_orthrus.tsv"), sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE) 

