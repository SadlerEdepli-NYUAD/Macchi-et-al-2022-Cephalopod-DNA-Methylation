##########################################################################################
### Figure 2I - 5C - S12 ###
##########################################################################################

# import annotation files from Trinotate
##########################################################################################
trinotate_annotation <- read.delim("./Input_Data/Trinotate_annotation/trinotate_annotation_report_new_ONLY_transcript.txt")

# preparation of dataframe 
##########################################################################################
library("dplyr")
trinotate_annotation_select <- trinotate_annotation %>% select(transcript_id, custom_db_name_BLASTX) 

# Remove the "transcript:" string in front of gene name in target_id column 
trinotate_annotation_select$transcript_id <- gsub("transcript:(.*)", "\\1", trinotate_annotation_select$transcript_id)

# Separate a column in multiple column base on special characters
library("tidyr")
trinotate_annotation_select_sep <- trinotate_annotation_select %>% 
  separate(custom_db_name_BLASTX, c("Prot_ID", "Species_ID", "Prot_ID.2", "Species_ID.2", "Alignment","Perc.Simil", "pValue", "FullName"), sep = "([_^])", extra = "merge", fill = "right")

View(trinotate_annotation_select_sep)
##########################################################################################

# Filter by a string pattern and then rename as unique transcript
##########################################################################################
# Histone Lysine methyltransferases
####################################
library("dplyr")
library("stringr")
trinotate_histone_lysine <- trinotate_annotation_select_sep %>% 
  filter(str_detect(FullName, regex(c("histone-lysine|histone methyltransferase"), ignore_case = TRUE))) %>%
  distinct(transcript_id, .keep_all = TRUE) %>% 
  arrange(transcript_id)

View(trinotate_histone_lysine)

# convert identical name as unique name adding .1, .2, .3, etc.
trinotate_histone_lysine$Prot_ID <- make.names(trinotate_histone_lysine$Prot_ID, unique=T)
####################################

# Histone acetyltransferase
####################################
library("dplyr")
library("stringr")
trinotate_histone_acetyl <- trinotate_annotation_select_sep %>% 
  filter(str_detect(FullName, regex("Histone acetyltransferase", ignore_case = TRUE))) %>%
  distinct(transcript_id, .keep_all = TRUE) %>% 
  arrange(transcript_id)

# convert identical name as unique name adding .1, .2, .3, etc.
trinotate_histone_acetyl$Prot_ID <- make.names(trinotate_histone_acetyl$Prot_ID, unique=T)
####################################

# Lysine-specific demethylase 
####################################
library("dplyr")
library("stringr")
trinotate_histone_demethylase <- trinotate_annotation_select_sep %>% 
  filter(str_detect(FullName, regex("Lysine-specific demethylase|Lysine-specific histone demethylase|Bifunctional arginine demethylase", ignore_case = TRUE))) %>%
  distinct(transcript_id, .keep_all = TRUE) %>% 
  arrange(transcript_id)

View(trinotate_histone_demethylase)

# convert identical name as unique name adding .1, .2, .3, etc.
trinotate_histone_demethylase$Prot_ID <- make.names(trinotate_histone_demethylase$Prot_ID, unique=T)
####################################

# Histone deacetylase 
####################################
library("dplyr")
library("stringr")
trinotate_histone_deacetylase <- trinotate_annotation_select_sep %>% 
  filter(str_detect(FullName, regex("Histone deacetylase", ignore_case = TRUE))) %>%
  distinct(transcript_id, .keep_all = TRUE) %>% 
  arrange(transcript_id)

# convert identical name as unique name adding .1, .2, .3, etc.
trinotate_histone_deacetylase$Prot_ID <- make.names(trinotate_histone_deacetylase$Prot_ID, unique=T)
####################################

# DNA methylation machinery 
####################################
library("dplyr")
library("stringr")
trinotate_histone_DNAmethylation <- trinotate_annotation_select_sep %>% 
  filter(str_detect(FullName, regex("CpG|cytosine|mismatch-specific thymine|ligase UHRF1|ligase UHRF2", ignore_case = TRUE))) %>%
  distinct(transcript_id, .keep_all = TRUE) %>% 
  arrange(transcript_id)

# changing inaccurate annotation of UHRF1 and UHRF2
trinotate_histone_DNAmethylation$Prot_ID <- gsub("UHRF1", "YDG_domain", trinotate_histone_DNAmethylation$Prot_ID)
trinotate_histone_DNAmethylation$Prot_ID.2 <- gsub("UHRF1", "YDG_domain (NOT UHRF1)", trinotate_histone_DNAmethylation$Prot_ID.2)
trinotate_histone_DNAmethylation$FullName <- gsub("E3 ubiquitin-protein ligase UHRF1", "YDG domain-containing protein (NOT E3 ubiquitin-protein ligase UHRF1)", trinotate_histone_DNAmethylation$FullName)
# changing inaccurate annotation of UHRF1 and UHRF2
trinotate_histone_DNAmethylation$Prot_ID <- gsub("UHRF2", "UHRF1", trinotate_histone_DNAmethylation$Prot_ID)
trinotate_histone_DNAmethylation$Prot_ID.2 <- gsub("UHRF2", "UHRF1 (NOT UHRF2)", trinotate_histone_DNAmethylation$Prot_ID.2)
trinotate_histone_DNAmethylation$FullName <- gsub("E3 ubiquitin-protein ligase UHRF2", "E3 ubiquitin-protein ligase UHRF1 (NOT UHRF2)", trinotate_histone_DNAmethylation$FullName)

View(trinotate_histone_DNAmethylation)

# convert identical name as unique name adding .1, .2, .3, etc.
trinotate_histone_DNAmethylation$Prot_ID <- make.names(trinotate_histone_DNAmethylation$Prot_ID, unique=T)
####################################
##########################################################################################

# extract and export genelist
##########################################################################################
load("./Input_Data/All_TPM_tissue_our_GeneID.RData")

# Histone Lysine methyltransferases
####################################
TPM_trinotate_histone_lysine <- merge(trinotate_histone_lysine, All_TPM_tissue_our_GeneID, by.x = "transcript_id", by.y = "target_id")
TPM_trinotate_histone_lysine_export <- subset(TPM_trinotate_histone_lysine, select = c(1:2))
write.table(TPM_trinotate_histone_lysine_export, file="./Input_Data/GeneLists/TPM_trinotate_histone_lysine.txt", quote = FALSE, sep="\t", row.names = FALSE)

# Histone acetyltransferase
####################################
TPM_trinotate_histone_acetyl <- merge(trinotate_histone_acetyl, All_TPM_tissue_our_GeneID, by.x = "transcript_id", by.y = "target_id")
TPM_trinotate_histone_acetyl_export <- subset(TPM_trinotate_histone_acetyl, select = c(1:2))
write.table(TPM_trinotate_histone_acetyl_export, file="./Input_Data/GeneLists/HistoneAcetylation_Genes.txt", quote = FALSE, sep="\t", row.names = FALSE)

# Histone demethylase
####################################
TPM_trinotate_histone_demethylase <- merge(trinotate_histone_demethylase, All_TPM_tissue_our_GeneID, by.x = "transcript_id", by.y = "target_id")
TPM_trinotate_histone_demethylase_export <- subset(TPM_trinotate_histone_demethylase, select = c(1:2))
write.table(TPM_trinotate_histone_demethylase_export, file="./Input_Data/GeneLists/HistoneDemethylase_Genes.txt", quote = FALSE, sep="\t", row.names = FALSE)

# Histone deacetylase
####################################
TPM_trinotate_histone_deacetylase <- merge(trinotate_histone_deacetylase, All_TPM_tissue_our_GeneID, by.x = "transcript_id", by.y = "target_id")
TPM_trinotate_histone_deacetylase_export <- subset(TPM_trinotate_histone_deacetylase, select = c(1:2))
write.table(TPM_trinotate_histone_deacetylase_export, file="./Input_Data/GeneLists/HistoneDeAcetylation_Genes.txt", quote = FALSE, sep="\t", row.names = FALSE)

# DNA methylation machinery
####################################
TPM_trinotate_histone_DNAmethylation <- merge(trinotate_histone_DNAmethylation, All_TPM_tissue_our_GeneID, by.x = "transcript_id", by.y = "target_id")
TPM_trinotate_histone_DNAmethylation_export <- subset(TPM_trinotate_histone_DNAmethylation, select = c(1:2))
write.table(TPM_trinotate_histone_DNAmethylation_export, file="./Input_Data/GeneLists/DNAmethylation_Genes.txt", quote = FALSE, sep="\t", row.names = FALSE)
##########################################################################################


# reordering column for having the supervised clustering from heatmap in Figure 1B
##########################################################################################
library("dplyr")
All_TPM_tissue_our_GeneID <- All_TPM_tissue_our_GeneID %>% 
  select(Testes_TPM, Retina_TPM, Subesophageal_brain_TPM, 
         Supraesophageal_TPM, Axial_nerve_cord_TPM, Optic_lobe_TPM, 
         Arm_B5_L1_2_TPM, Arm_B5_L1_3_TPM, Arm_B5_L1_1_TPM, 
         Junior_OBJ1_TPM, Skin_TPM, Suckers_TPM, Viscera_TPM, 
         Ova_TPM, Posterior_salivary_gland_TPM, target_id) 

View(All_TPM_tissue_our_GeneID)

# import specific gene list (manually curated with gene functions) from list exported above
##########################################################################################


# DNA methylation (Fig 2I)
##########################################################################################
DNAmethylation_Genes <- read.table("Input_Data/GeneLists/DNAmethylation_Genes.txt", sep="\t", header = T)

# merging with sort=FALSE keep the order of X list
All_TPM_tissue_our_DNAmethylation_Genes <- merge(DNAmethylation_Genes, All_TPM_tissue_our_GeneID, by = "target_id", sort = FALSE)
write.table(All_TPM_tissue_our_DNAmethylation_Genes, file="./Output_Data/heatmap_epigenetics_TPMs/All_TPM_tissue_our_DNAmethylation_Genes.txt", quote = FALSE, sep="\t", row.names = FALSE)

# for using in the heatmap
rownames(All_TPM_tissue_our_DNAmethylation_Genes) <- All_TPM_tissue_our_DNAmethylation_Genes$ALIAS
All_TPM_tissue_our_DNAmethylation_Genes_select <- subset(All_TPM_tissue_our_DNAmethylation_Genes, select = -c(target_id, ALIAS, Function) )
All_TPM_tissue_our_DNAmethylation_Genes_select_matrix <- as.matrix(All_TPM_tissue_our_DNAmethylation_Genes_select)

# z score normalization
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
# apply function and transpose matrix
All_TPM_tissue_our_DNAmethylation_Genes_select_matrix_norm <- t(apply(All_TPM_tissue_our_DNAmethylation_Genes_select_matrix, 1, cal_z_score))

# heatmap
############################################
library(pheatmap)
postscript("./Results/Epigenetic_Machinery/DNAmethylation_Genes_heatmap.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 7.5, height = 4.35)
pheatmap(All_TPM_tissue_our_DNAmethylation_Genes_select_matrix_norm, cluster_rows=FALSE, show_rownames=T, cellwidth=10) 
dev.off()
##########################################################################################

# Histone Methyl-WriterEraser (Fig.5C)
##########################################################################################
HistoneMethylation_WriterEraser <- read.table("Input_Data/GeneLists/HistoneMethylation_WriterEraser.txt", sep="\t", header = T)

# merging with sort=FALSE keep the order of X list
All_TPM_tissue_our_HistoneMethylation_WriterEraser <- merge(HistoneMethylation_WriterEraser, All_TPM_tissue_our_GeneID, by = "target_id", sort = FALSE)
write.table(All_TPM_tissue_our_HistoneMethylation_WriterEraser, file="./Output_Data/heatmap_epigenetics_TPMs/All_TPM_tissue_our_HistoneMethylation_WriterEraser.txt", quote = FALSE, sep="\t", row.names = FALSE)

# for using in the heatmap
rownames(All_TPM_tissue_our_HistoneMethylation_WriterEraser) <- All_TPM_tissue_our_HistoneMethylation_WriterEraser$ALIAS
All_TPM_tissue_our_HistoneMethylation_WriterEraser_select <- subset(All_TPM_tissue_our_HistoneMethylation_WriterEraser, select = -c(target_id, ALIAS, Function) )
All_TPM_tissue_our_HistoneMethylation_WriterEraser_select_matrix <- as.matrix(All_TPM_tissue_our_HistoneMethylation_WriterEraser_select)

# z score normalization
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
# apply function and transpose matrix
All_TPM_tissue_our_HistoneMethylation_WriterEraser_select_matrix_norm <- t(apply(All_TPM_tissue_our_HistoneMethylation_WriterEraser_select_matrix, 1, cal_z_score))

# heatmap
############################################
library(pheatmap)
postscript("./Results/Epigenetic_Machinery/HistoneMethylation_WriterEraser_heatmap_NOclust.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 7.5, height = 6.5)
pheatmap(All_TPM_tissue_our_HistoneMethylation_WriterEraser_select_matrix_norm, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=T, cellwidth=10) 
dev.off()
##########################################################################################

# Histone Lysine methyltransferases extended heatmap (Fig.S12A)
##########################################################################################
# for using in the heatmap
rownames(TPM_trinotate_histone_lysine) <- TPM_trinotate_histone_lysine$Prot_ID
TPM_trinotate_histone_lysine_select <- subset(TPM_trinotate_histone_lysine, select = -c(1:9))
TPM_trinotate_histone_lysine_select_matrix <- as.matrix(TPM_trinotate_histone_lysine_select)

# z score normalization
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
# apply function and transpose matrix
TPM_trinotate_histone_lysine_select_matrix_norm <- t(apply(TPM_trinotate_histone_lysine_select_matrix, 1, cal_z_score))

# heatmap
############################################
library(pheatmap)
postscript("./Results/Epigenetic_Machinery/histone_lysine_methyltransferase_heatmap_RowClust.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 7.5, height = 14.25)
pheatmap(TPM_trinotate_histone_lysine_select_matrix_norm, cutree_rows = 13, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=T, cellwidth=10)
dev.off()
##########################################################################################

# Histone Acetyl-WriterEraser (Fig.S12B)
##########################################################################################
HistoneAcetylation_WriterEraser <- read.table("Input_Data/GeneLists/HistoneAcetylation_WriterEraser.txt", sep="\t", header = T)

# merging with sort=FALSE keep the order of X list
All_TPM_tissue_our_HistoneAcetylation_WriterEraser <- merge(HistoneAcetylation_WriterEraser, All_TPM_tissue_our_GeneID, by = "target_id", sort = FALSE)
write.table(All_TPM_tissue_our_HistoneAcetylation_WriterEraser, file="./Output_Data/heatmap_epigenetics_TPMs/All_TPM_tissue_our_HistoneAcetylation_WriterEraser.txt", quote = FALSE, sep="\t", row.names = FALSE)

# for using in the heatmap
rownames(All_TPM_tissue_our_HistoneAcetylation_WriterEraser) <- All_TPM_tissue_our_HistoneAcetylation_WriterEraser$ALIAS
All_TPM_tissue_our_HistoneAcetylation_WriterEraser_select <- subset(All_TPM_tissue_our_HistoneAcetylation_WriterEraser, select = -c(target_id, ALIAS, Function) )
All_TPM_tissue_our_HistoneAcetylation_WriterEraser_select_matrix <- as.matrix(All_TPM_tissue_our_HistoneAcetylation_WriterEraser_select)

# z score normalization
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
# apply function and transpose matrix
All_TPM_tissue_our_HistoneAcetylation_WriterEraser_select_matrix_norm <- t(apply(All_TPM_tissue_our_HistoneAcetylation_WriterEraser_select_matrix, 1, cal_z_score))

# heatmap
############################################
library(pheatmap)

postscript("./Results/Epigenetic_Machinery/HistoneAcetylation_WriterEraser_heatmap_NOclust.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 7.5, height = 6.25)
pheatmap(All_TPM_tissue_our_HistoneAcetylation_WriterEraser_select_matrix_norm, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=T, cellwidth=10) 
dev.off()
##########################################################################################

