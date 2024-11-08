library(orthrus)
library("gemini")
library(dplyr)
data("counts", "guide.annotation", "sample.replicate.annotation", package = "gemini")
library(readxl)
library(tidyr)
library(dplyr)
library(tibble) 
library(stringr) 



chymera_genes = c(chymera_paralog$gene1, chymera_paralog$gene2) %>% unique()
pc_genes = pc_genes[which((pc_genes %in% chymera_genes))]

chymera <- chymera_paralog # get directly from Orthrus package
rownames(chymera) = paste0(chymera$Cas9.Guide, ";", chymera$Cpf1.Guide )
chymera = chymera %>%
  filter(gene1 != "NT" & gene2 != "NT")
chymera = chymera %>%
  filter (! (gene1 == "NegControl" & gene2 == "NegControl"))

chymera = chymera %>%
  filter (gene1 != "None" & gene2 != "None")

guide.annotation_chymera = data.frame(rowname = rownames(chymera) , gene1 = chymera$gene1, gene2 = chymera$gene2, gene1.guide = chymera$Cas9.Guide, gene2.guide = chymera$Cpf1.Guide )
chymera_HAP = chymera[,-c(1:12, 19:24,32)]
chymera_HAP = data.matrix(chymera_HAP)
head(chymera_HAP)

guide.annotation_chymera = data.frame(rowname = rownames(chymera) , gene1 = chymera$gene1, gene2 = chymera$gene2, gene1.guide = chymera$Cas9.Guide, gene2.guide = chymera$Cpf1.Guide )
chymera_RPE = chymera[,-c(1:18, 25:31)]
chymera_RPE = data.matrix(chymera_RPE)
head(chymera_RPE)




## Run Separately separately for HAP and RPE
sample.replicate.annotation_chymera_HAP <- data.frame(
  colname = c("HAP1.Torin.T12A", "HAP1.Torin.T12B" ,"HAP1.Torin.T12C", "HAP1.Torin.T18A", "HAP1.Torin.T18B",
              "HAP1.Torin.T18C",      
              "HAP1.T12A","HAP1.T12B","HAP1.T12C","HAP1.T18A","HAP1.T18B","HAP1.T18C","HAP1.T0"),
  samplename = c("HAP1.Torin.T12", "HAP1.Torin.T12" ,"HAP1.Torin.T12", "HAP1.Torin.T18", "HAP1.Torin.T18",
                 "HAP1.Torin.T18",      
                 "HAP1.T12","HAP1.T12","HAP1.T12","HAP1.T18","HAP1.T18","HAP1.T18","HAP1.T0"),
  replicate = c("A" ,"B" ,"C", "A", "B",
                "C",       
                "A","B","C","A","B","C", NA)
)

sample.replicate.annotation_chymera_RPE <- data.frame(
  colname = c( "RPE1.T18A", "RPE1.T18B", "RPE1.T18C", "RPE1.T24A", "RPE1.T24B", "RPE1.T24C",
               "RPE1.T0"),
  samplename = c("RPE1.T18", "RPE1.T18", "RPE1.T18", "RPE1.T24", "RPE1.T24", "RPE1.T24","RPE1.T0"),
  replicate = c("A", "B", "C", "A", "B", "C", NA)
)

nc_genes = "NegControl"


## HAP

Input_HAP <- gemini_create_input(counts.matrix = chymera_HAP,
                                 sample.replicate.annotation = sample.replicate.annotation_chymera_HAP,
                                 guide.annotation = guide.annotation_chymera,
                                 ETP.column = 'HAP1.T0', 
                                 gene.column.names = c("gene1", "gene2"),
                                 sample.column.name = "samplename",
                                 verbose = TRUE)
Input_HAP %<>% gemini_calculate_lfc(normalize = TRUE, 
                                    CONSTANT = 32)


Model_HAP <- gemini_initialize(Input = Input_HAP, 
                               nc_gene = "NegControl", 
                               pattern_join = ';',
                               pattern_split = ';', 
                               cores = 1,
                               verbose = TRUE)
Model_HAP %<>% gemini_inference(cores = 1,
                                verbose = FALSE)


gemini_plot_mae(Model_HAP)

Score_HAP1 <- gemini_score(Model = Model_HAP, pc_threshold = -Inf) # pc_threshold removes the essential_gene


write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Chymera_HAP1_Strong.csv", Score_HAP1$strong)
write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Chymera_HAP1_sensitive_lethality.csv", Score_HAP1$sensitive_lethality)
write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Chymera_HAP1_sensitive_recovery.csv", Score_HAP1$sensitive_recovery)

###


## RPE

Input_RPE <- gemini_create_input(counts.matrix = chymera_RPE,
                                 sample.replicate.annotation = sample.replicate.annotation_chymera_RPE,
                                 guide.annotation = guide.annotation_chymera,
                                 ETP.column = 'RPE1.T0', 
                                 gene.column.names = c("gene1", "gene2"),
                                 sample.column.name = "samplename",
                                 verbose = TRUE)
Input_RPE %<>% gemini_calculate_lfc(normalize = TRUE, 
                                    CONSTANT = 32)


Model_RPE <- gemini_initialize(Input = Input_RPE, 
                               nc_gene = "NegControl", 
                               pattern_join = ';',
                               pattern_split = ';', 
                               cores = 1,
                               verbose = TRUE)



Model_RPE %<>% gemini_inference(cores = 1,
                                verbose = FALSE)

############################################
gemini_plot_mae(Model_RPE)

Score_RPE <- gemini_score(Model = Model_RPE,  pc_threshold = -Inf)


write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Chymera_RPE1_Strong.csv", Score_RPE$strong)
write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Chymera_RPE1_sensitive_lethality.csv", Score_RPE$sensitive_lethality)
write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Chymera_RPE1_sensitive_recovery.csv", Score_RPE$sensitive_recovery)

