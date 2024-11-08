## take data directly from Orthrus manual , copy files


library(orthrus) 

df <- chymera_paralog 

output_folder <- file.path("Orthrus Scripts/OrthrusOutput/") 




#Call the add_screens_from_table function to build up a list of screens with names 
# and corresponding technical replicates, starting with T0 replicates.


sample_table <- chymera_sample_table 
screens <- add_screens_from_table(sample_table) 



# Now we need to normalize each screen in three different ways:
#   
# To their respective T0 screens by computing log fold-changes (LFCs)
# To the respective depth of each technical replicate
# The function normalize_screens automatically performs all of these normalization steps. The function infers which columns of df need to be normalized to which T0 screens based on the normalize_name parameter of each screen in screens (screens without this optional parameter will not be normalized to other screens). Log-scaling and depth-normalization is performed on each screen regardless of the normalize_name parameter. For example, after normalization T0 columns in df will contain log-scaled, depth-normalized read counts, whereas columns from later timepoints will contain depth-normalized LFCs compared to their respective T0s.
# We do not remove any guides based on low/high read counts
df <- normalize_screens(df, screens, filter_names = c("HAP1_T0", "RPE1_T0"), min_reads = 0, max_reads = 10000) 



# The last thing we need to do before scoring data is parse it into a different structure and split guides by their type, since we score dual-targeting guides separately from combinatorial-targeting guides.
# 
guides <- split_guides(df, screens, "Cas9.Guide", "Cpf1.Guide")
dual <- guides[["dual"]]
single <- guides[["single"]]
paralogs <- guides[["combn"]]

screens_to_score <- c("HAP1_T12", "HAP1_T18", "RPE1_T18", "RPE1_T24") 

temp <- score_combn_vs_single(paralogs, single, screens,
                              screens_to_score, test = "moderated-t", 
                              return_residuals = TRUE, filter_genes = c("NT")) 
paralog_scores <- temp[["scored_data"]] 
write.table(paralog_scores, file.path(output_folder,
                                      "chymera_orthrus.tsv"), sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE) 


