library(orthrus)
library(stringr) 
library(dplyr)
library(tidyr)

## In a separate file, I ran the same code by not removing the Fluc_FLuc pairs. and also ran the same code by changing
# filter_genes from "FLuc" to NULL. These changes do not make any difference to the results. All are highly correlated.

output_folder <- file.path("Orthrus Scripts/OrthrusOutput/")

if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }



df <- read.csv("InputData/Thompson/raw_counts.txt", sep = "\t", 
               header = TRUE, stringsAsFactors = FALSE) 
df <- df %>%  drop_na() %>%
  select(-c("gRNA")) %>% 
  filter(GENE != "Non_Tar_RND_control_Non_Tar_RND_control") %>%
  separate(GENE, into = c("gene1", "gene2"), sep = "_") %>%
  mutate(gene1 = ifelse(gene1 == "Fluc", "NegControl", gene1)) %>%
  mutate(gene2 = ifelse(gene2 == "Fluc", "NegControl", gene2)) %>%
  mutate(Guide1 = gene1, Guide2 = gene2)

head(df)

## Applying the same hack as we did for Parrish study
df2 = df %>% mutate(gene3 = gene1, gene4 = gene2) %>% select(-c(gene1,gene2)) %>% rename(gene1 = gene4, gene2 = gene3) %>% mutate(Guide3 = Guide1, Guide4 = Guide2) %>% select(-c(Guide1,Guide2)) %>% rename(Guide1 = Guide4, Guide2 = Guide3)
df2 = df2 %>% rbind(df)
df = df2



sample_table <- data.frame(
  Screen = c("T0", "A375_D14", "A375_D28", "MEWO_D14", "MEWO_D28", "RPE_D14", "RPE_D28"),
  Replicates = c(
    "Control_R1;Control_R2;Control_R3", 
    "A375_D14_R1;A375_D14_R2;A375_D14_R3", 
    "A375_D28_R1;A375_D28_R2;A375_D28_R3", 
    "MEWO_D14_R1;MEWO_D14_R2;MEWO_D14_R3", 
    "MEWO_D28_R1;MEWO_D28_R2;MEWO_D28_R3", 
    "RPE_D14_R1;RPE_D14_R2;RPE_D14_R3", 
    "RPE_D28_R1;RPE_D28_R2;RPE_D28_R3"
  ),
  NormalizeTo = c(NA, "T0", "T0", "T0", "T0", "T0", "T0")
)

screens <- add_screens_from_table(sample_table)

df <- normalize_screens(df, screens, filter_names = "T0", min_reads = 0, max_reads =  Inf) 

guides <- split_guides(df, screens, "Guide1", "Guide2") 

dual <- guides[["dual"]] # When they trget same gene twice
single <- guides[["single"]] 
combn <- guides[["combn"]] 


score_combn_batch(combn, single, screens, batch_file, output_folder, test = "moderated-t", loess = FALSE, filter_genes = c("NegControl"), neg_type = "Sensitizer", pos_type = "Suppressor", fdr_threshold = 0.2, differential_threshold = 0.5) 

screens_to_score <- c("A375_D14", "MEWO_D14", "RPE_D14", "A375_D28", "MEWO_D28", "RPE_D28") 

temp <- score_combn_vs_single(combn, single, screens,
                              screens_to_score, test = "moderated-t", 
                              return_residuals = TRUE) 
paralog_scores <- temp[["scored_data"]]
write.table(paralog_scores, file.path(output_folder,
                                      "thompson_orthrus.tsv"), sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE) 


