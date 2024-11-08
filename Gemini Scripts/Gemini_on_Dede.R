library(orthrus)
library(stringr)
library(gemini)
library(tibble)

dede = read.csv("InputData\\Dede\\Dede.txt", sep = "\t")

rownames_ = strsplit(dede[,1], "_")
temp_a <- data.frame(do.call(rbind, rownames_))
rownames_ = paste0(temp_a[,2], ";", temp_a[,4])
rownames(dede) = rownames_

genes = strsplit(dede[,2], "[:.]")
temp_b <- data.frame(do.call(rbind, genes))


guide.annotation_dede = data.frame(rowname = rownames_, gene1 = temp_b[,1], gene2 = temp_b[,3], gene1.guide = temp_a[,2], gene2.guide = temp_a[,4] )
dede = dede[,-c(1:2)]
dede = data.matrix(dede)
head(dede)





sample.replicate.annotation_dede <- data.frame(
  colname = c("A549.T2A.Ex", "A549.T2B.Ex", "A549.T2C.Ex", "HT29.T2A.Ex",
              "HT29.T2B.Ex", "HT29.T2C.Ex", "OVCAR8.T2A.Ex", "OVCAR8.T2B.Ex", "OVCAR8.T2C.Ex", "plasmid.T0.Ex"),
  samplename = c("A549", "A549", "A549", "HT29", "HT29", "HT29", "OVCAR8", "OVCAR8", "OVCAR8", "plasmid.T0.Ex"),
  replicate = c("T2A.Ex", "T2B.Ex", "T2C.Ex", "T2A.Ex",
                "T2B.Ex", "T2C.Ex", "T2A.Ex", "T2B.Ex", "T2C.Ex", NA)
)

nc_genes = read.csv(file = "InputData\\Dede\\control_nonessentials.csv")
nc_genes = unname(unlist(nc_genes,recursive = TRUE))

pc_genes = read.csv(file = "InputData\\Dede\\control_essentials.csv")
pc_genes = unname(unlist(pc_genes,recursive = TRUE))

Input <- gemini_create_input(counts.matrix = dede,
                             sample.replicate.annotation = sample.replicate.annotation_dede,
                             guide.annotation = guide.annotation_dede,
                             ETP.column = 'plasmid.T0.Ex', 
                             gene.column.names = c("gene1", "gene2"),
                             sample.column.name = "samplename",
                             verbose = TRUE)
Input %<>% gemini_calculate_lfc(normalize = TRUE, 
                                CONSTANT = 32)


Model <- gemini_initialize(Input = Input, 
                           nc_gene = nc_genes, 
                           pattern_join = ';',
                           pattern_split = ';', 
                           cores = 1,
                           verbose = TRUE)
save(Model,file =  "Dede_Model.RData")
Model %<>% gemini_inference(cores = 1,
                            verbose = FALSE)


gemini_plot_mae(Model)

Score<- gemini_score(Model = Model,
                      pc_threshold = -Inf) # To skip the essential gene filter, pc_threshold has been set to -Inf


write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Dede_Strong.csv", Score$strong)
write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Dede_Sensitive_Lethality.csv", Score$sensitive_lethality)
write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_Dede_Sensitive_Recovery.csv", Score$sensitive_recovery)
