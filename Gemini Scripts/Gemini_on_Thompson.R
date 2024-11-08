library(orthrus)
library("gemini")
library(dplyr)
library(tidyr)
library(stringr) 

thompson = read.csv("InputData/Thompson/raw_counts.txt", sep = "\t")
thompson <- thompson[seq(2, nrow(thompson), by = 2), ]
thompson_guides = readxl::read_excel("InputData/Thompson//Guide_Sequences.xlsx")



thompson = thompson %>%
  separate(gRNA, into = c("Guide1", "Guide2", "type"), sep = "_N_|_SL_oligo:", remove = FALSE, extra = "merge") %>%
  separate(GENE, into = c("gene1", "gene2"), sep = "_", remove = FALSE, extra = "merge") %>%
  mutate(gene1 = ifelse(type == "Non_targt_RND_control", "Fluc", gene1)) %>%
  mutate(gene2 = ifelse(type == "Non_targt_RND_control", "Fluc", gene2)) %>%
  merge(thompson_guides, by.x = c("Guide1", "Guide2"), by.y = c("gRNA_A_IDs", "gRNA_B_IDs" )) 

# Did not have to use it
#duplicates <- thompson %>%
#  group_by(gRNA_A_seq, gRNA_B_seq) %>%
#  filter(n() > 1) %>%
#  arrange(gRNA_A_seq, gRNA_B_seq) %>%
#  ungroup() 
#duplicates

## IMPORTANT STEP, here you will see that some combinations of guide sequences are duplicate, that hinders the dataset as rownames that we are 
#trying to set are not unique. Upoin further invetigation, I saw that the duplicated rows have all 0's in them. So better we remove them, 
# 168 of such rows are duplicares for guide seqs. around 84 will have 0's#



#thompson <- thompson %>%
  # Select only numeric columns
  # Filter rows where all selected numeric columns are 0
#  filter(if_all(where(is.numeric), ~. != 0))
# Now its not causing issue
rownames(thompson) = paste0(thompson$gRNA_A_seq,"_",thompson$Guide1, ":", thompson$gRNA_B_seq, "_",thompson$Guide2)


guide.annotation_thompson = thompson %>% select(c("gene1", "gene2", "gRNA_A_seq", "gRNA_B_seq")) %>%
  rename(gene1.guide = gRNA_A_seq, gene2.guide = gRNA_B_seq ) %>%
  tibble::rownames_to_column() 

thompson = thompson %>%
  select(where(is.numeric)) %>%
  select(c("A375_D14_R1", "A375_D14_R2", "A375_D14_R3", 
           "MEWO_D14_R1", "MEWO_D14_R2", "MEWO_D14_R3",
           "RPE_D14_R1", "RPE_D14_R2", "RPE_D14_R3" , 
           "A375_D28_R1", "A375_D28_R2", "A375_D28_R3",
           "MEWO_D28_R1" ,"MEWO_D28_R2", "MEWO_D28_R3",
           "RPE_D28_R1", "RPE_D28_R2", "RPE_D28_R3",
           "Control_R1", "Control_R2", "Control_R3"))

thompson = as.matrix(thompson)


sample.replicate.annotation_thompson <- data.frame(
  colname = c("A375_D14_R1", "A375_D14_R2", "A375_D14_R3", 
              "MEWO_D14_R1", "MEWO_D14_R2", "MEWO_D14_R3",
              "RPE_D14_R1", "RPE_D14_R2", "RPE_D14_R3" , 
              "A375_D28_R1", "A375_D28_R2", "A375_D28_R3",
              "MEWO_D28_R1" ,"MEWO_D28_R2", "MEWO_D28_R3",
              "RPE_D28_R1", "RPE_D28_R2", "RPE_D28_R3",
              "Control_R1", "Control_R2", "Control_R3"),
  samplename = c("A375_D14", "A375_D14", "A375_D14",
                 "MEWO_D14", "MEWO_D14", "MEWO_D14",
                 "RPE_D14", "RPE_D14", "RPE_D14",
                 "A375_D28", "A375_D28", "A375_D28",
                 "MEWO_D28", "MEWO_D28", "MEWO_D28",
                 "RPE_D28", "RPE_D28", "RPE_D28",
                 "Control",  "Control",  "Control"),
  replicate = c("R1", "R2", "R3" ,"R1" ,"R2", "R3", "R1",  "R2",  "R3",
                "R1", "R2", "R3" ,"R1" ,"R2", "R3", "R1",  "R2",  "R3",
                "R1", "R2", "R3")
)


nc_genes = c("Fluc")
LTP = c("A375_D14_R1", "A375_D14_R2", "A375_D14_R3",
        "A375_D28_R1", "A375_D28_R2", "A375_D28_R3",
        "MEWO_D14_R1", "MEWO_D14_R2", "MEWO_D14_R3",
        "MEWO_D28_R1", "MEWO_D28_R2","MEWO_D28_R3",
        "RPE_D14_R1", "RPE_D14_R2", "RPE_D14_R3",
        "RPE_D28_R1","RPE_D28_R2","RPE_D28_R3")
ETP = grep("Control_", colnames(thompson))
LTP = 1:18
Input <- gemini_create_input(counts.matrix = thompson ,
                             sample.replicate.annotation = sample.replicate.annotation_thompson,
                             guide.annotation = guide.annotation_thompson,
                             ETP.column = ETP, 
                             LTP.column = LTP,
                             gene.column.names = c("gene1", "gene2"),
                             sample.column.name = "samplename",
                             verbose = TRUE)  ##
Input %<>% gemini_calculate_lfc(normalize = TRUE, 
                                CONSTANT = 32)
 
Model <- gemini_initialize(Input = Input, 
                           nc_gene = c("Fluc"), 
                           pattern_join = ':',
                           pattern_split = ':', 
                           cores = 1,
                           verbose = TRUE)

Model %<>% gemini_inference(cores = 1,
                            verbose = TRUE, force_results = TRUE)

gemini_plot_mae(Model)
Score <- gemini_score(Model = Model,pc_gene = essentials)

write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Thompson_Strong.csv", Score$strong)
write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Thompson_Sensitive.csv", Score$sensitive_lethality)
write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Thompson_sensitive_recovery.csv", Score$sensitive_recovery)

