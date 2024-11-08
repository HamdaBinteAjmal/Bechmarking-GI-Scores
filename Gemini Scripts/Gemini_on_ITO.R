library("gemini")
library(tidyr)
library(dplyr)

ito <- readxl::read_excel("InputData/Ito/Count_data_ParalogV1.xlsx") 
ito <- ito %>%  separate(`Left-sgRNA_Right-sgRNA`, into = c("gene1_guide", "gene2_guide"), sep = "_", remove = FALSE)
ito_genes = c(ito$Aureus_gene, ito$Pyogenes_gene) %>% unique()
pc_genes = pc_genes[pc_genes %in% ito_genes]
rownames(ito) =ito$`Left-sgRNA_Right-sgRNA`

head(ito)
sample.replicate.annotation_ito <- data.frame(
  colname = c("Meljuso_RepA", "Meljuso_RepB", "Meljuso_RepC",          
              "GI1_004_RepA","GI1_004_RepB","GI1_004_RepC",          
              "MEL202_003_RepA","MEL202_003_RepB","MEL202_003_RepC",       
              "PK1_REPA","PK1_REPB","PK1_REPC",              
              "MEWO_REPA","MEWO_REPB","MEWO_REPC",             
              "HS944T_REPA","HS944T_REPB","HS944T_REPC" ,          
              "IPC298_REPA","IPC298_REPB", "IPC298_REPC",           
              "A549_REPA", "A549_REPB", "A549_REPC",             
              "HSC5_REPA", "HSC5_REPB", "HSC5_REPC",          
              "HS936T_REPA", "HS936T_REPB", "HS936T_REPC",        
              "PATU8988S_REPA", "PATU8988S_REPB", 
              "pDNA"),
  samplename = c("Meljuso", "Meljuso", "Meljuso",          
                 "GI1_004","GI1_004","GI1_004",          
                 "MEL202_003","MEL202_003","MEL202_003",       
                 "PK1","PK1","PK1",              
                 "MEWO","MEWO","MEWO",             
                 "HS944T","HS944T","HS944T" ,          
                 "IPC298","IPC298", "IPC298",           
                 "A549", "A549", "A549",             
                 "HSC5", "HSC5", "HSC5",          
                 "HS936T", "HS936T", "HS936T",        
                 "PATU8988S", "PATU8988S", 
                 "pDNA"),
  replicate = c("RepA", "RepB", "RepC",          
                "RepA","RepB","RepC",          
                "RepA","RepB","RepC",       
                "REPA","REPB","REPC",              
                "REPA","REPB","REPC",             
                "REPA","REPB","REPC" ,          
                "REPA","REPB", "REPC",           
                "REPA", "REPB", "REPC",             
                "REPA", "REPB", "REPC",          
                "REPA", "REPB", "REPC",        
                "REPA", "REPB", 
                NA)
)

guide.annotation_ito = data.frame(rowname = rownames(ito), Aureus_gene = ito$Aureus_gene,
                                  Pyogenes_gene = ito$Pyogenes_gene, gene1.guide = ito$gene1_guide, gene2.guide = ito$gene2_guide )

counts_ito = ito %>% dplyr::select(-c("Aureus_gene","Pyogenes_gene","Left-sgRNA_Right-sgRNA","gene1_guide","gene2_guide" ))
counts_ito = as.matrix(counts_ito) 
rownames(counts_ito) = rownames(ito)

Input <- gemini_create_input(counts.matrix = counts_ito,
                             sample.replicate.annotation = sample.replicate.annotation_ito,
                             guide.annotation = guide.annotation_ito,
                             ETP.column = 'pDNA', 
                             gene.column.names = c("Aureus_gene", "Pyogenes_gene"),
                             sample.column.name = "samplename",
                             verbose = TRUE)
Input %<>% gemini_calculate_lfc(normalize = TRUE, 
                                CONSTANT = 32)


Model <- gemini_initialize(Input = Input, 
                           nc_gene = c("AAVS1"), 
                           pattern_join = ';',
                           pattern_split = '_', 
                           cores = 1,
                           verbose = TRUE)
Model %<>% gemini_inference(cores = 1,
                            verbose = FALSE)
gemini_plot_mae(Model)

Score <- gemini_score(Model = Model,pc_threshold = -Inf)

write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_ITO_Strong.csv", Score$strong)
write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_ITO_Sensitive.csv", Score$sensitive_lethality)
write.csv(file = "Gemini Scripts/GeminiOutput/Gemini_ITO_sensitive_recovery.csv", Score$sensitive_recovery)
