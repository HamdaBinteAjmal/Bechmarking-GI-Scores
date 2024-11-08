library(orthrus)
library(stringr) 
library(dplyr)

output_folder <- file.path("Orthrus Scripts/OrthrusOutput/")


if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }



df <- read.csv(file.path("InputData/Dede/counts.txt"), sep = "\t", 
               header = TRUE, stringsAsFactors = FALSE) 
essentials <- read.csv(file.path("InputData/Dede/pan-species-control-essentials-50genes.txt"), 
                       header = TRUE, stringsAsFactors = FALSE,sep = "\t") 
nonessentials <- read.csv(file.path("InputData/Dede/pan-species-control-nonessentials-50genes.txt"), 
                          header = TRUE, stringsAsFactors = FALSE,sep = "\t") 


essentials <- unlist(essentials) 
nonessentials <- unlist(nonessentials) 

split <- str_split_fixed(df$GENE, ":", 2) 
df$gene1 <- gsub("\\..*", "", split[,1]) 
df$gene2 <- gsub("\\..*", "", split[,2])
df$gene1[df$gene1 %in% nonessentials] <- "NegControl"
df$gene2[df$gene2 %in% nonessentials] <- "NegControl" 

# Adds guide columns 
split <- str_split_fixed(df$GENE_CLONE, "_", 4)
df$Guide1 <- split[,2] 
df$Guide2 <- split[,4] 



sample_table <- data.frame(
  Screen = c("T0", "A549", "HT29", "OVCAR8"),
  Replicates = c("plasmid.T0.Ex", 
                 "A549.T2A.Ex;A549.T2B.Ex;A549.T2C.Ex", 
                 "HT29.T2A.Ex;HT29.T2B.Ex;HT29.T2C.Ex", 
                 "OVCAR8.T2A.Ex;OVCAR8.T2B.Ex;OVCAR8.T2C.Ex"),
  NormalizeTo = c(NA, "T0", "T0", "T0")
)


# Processes data 
screens <- add_screens_from_table(sample_table)
# Do not filter based on min counts or max counts

df <- normalize_screens(df, screens, filter_names = "T0", min_reads = 0, max_reads = Inf) 

df = df %>% select(-c("GENE_CLONE", "GENE"))

guides <- split_guides(df, screens, "Guide1", "Guide2") 

dual <- guides[["dual"]] 
single <- guides[["single"]] 
combn <- guides[["combn"]] 


score_combn_batch(combn, single, screens, batch_file, output_folder, test = "moderated-t", loess = FALSE, filter_genes = c("NegControl"), neg_type = "Sensitizer", pos_type = "Suppressor", fdr_threshold = 0.2, differential_threshold = 0.5) 
## Alternative way

screens_to_score <- c("A549", "HT29", "OVCAR8") 

temp <- score_combn_vs_single(combn, single, screens,
                              screens_to_score, test = "moderated-t", 
                              return_residuals = TRUE, filter_genes = c("NegControl")) 
paralog_scores <- temp[["scored_data"]] 
write.table(paralog_scores, file.path(output_folder,
                                      "dede_orthrus.tsv"), sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE) 

