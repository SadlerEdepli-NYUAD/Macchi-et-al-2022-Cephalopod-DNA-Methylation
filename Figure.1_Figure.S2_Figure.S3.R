##########################################################################################
### Figure 1A ###
##########################################################################################

# import all file *tsv with TPMs
##########################################################################################
Axial_nerve_cord_TPM <- read.table("Input_Data/Octopus_Reads_TPM/Axial_nerve_cord_abundance.tsv", sep="\t", header = T)
Optic_lobe_TPM <- read.table("Input_Data/Octopus_Reads_TPM/Optic_lobe_abundance.tsv", sep="\t", header = T)
Ova_TPM <- read.table("Input_Data/Octopus_Reads_TPM/Ova_abundance.tsv", sep="\t", header = T)
Posterior_salivary_gland_TPM <- read.table("Input_Data/Octopus_Reads_TPM/posterior_salivary_gland_abundance.tsv", sep="\t", header = T)
Retina_TPM <- read.table("Input_Data/Octopus_Reads_TPM/Retina_abundance.tsv", sep="\t", header = T)
Skin_TPM <- read.table("Input_Data/Octopus_Reads_TPM/skin_abundance.tsv", sep="\t", header = T)
Subesophageal_brain_TPM <- read.table("Input_Data/Octopus_Reads_TPM/Subesophageal_brain_abundance.tsv", sep="\t", header = T)
Suckers_TPM <- read.table("Input_Data/Octopus_Reads_TPM/Suckers_abundance.tsv", sep="\t", header = T)
Supraesophageal_TPM <- read.table("Input_Data/Octopus_Reads_TPM/Supraesophageal_abundance.tsv", sep="\t", header = T)
Testes_TPM <- read.table("Input_Data/Octopus_Reads_TPM/Testes_abundance.tsv", sep="\t", header = T)
Viscera_TPM <- read.table("Input_Data/Octopus_Reads_TPM/Viscera_abundance.tsv", sep="\t", header = T)

# L1 Arm 1=distal, 2=medial, 3=proximal
B5_L1_1_abundance_TPM <- read.table("Input_Data/Our_RNA-seq/B5_L1_1_abundance.tsv", sep="\t", header = T)
B5_L1_2_abundance_TPM <- read.table("Input_Data/Our_RNA-seq/B5_L1_2_abundance.tsv", sep="\t", header = T)
B5_L1_3_abundance_TPM <- read.table("Input_Data/Our_RNA-seq/B5_L1_3_abundance.tsv", sep="\t", header = T)

# Hatchling
OBJ1_abundance_TPM <- read.table("Input_Data/Our_RNA-seq/OBJ1-Bsby_abundance.tsv", sep="\t", header = T)
##########################################################################################

# check if tables are ordered in the same way before "cbind" the TPMs columns
##########################################################################################
identical(Axial_nerve_cord_TPM$target_id, Optic_lobe_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, Ova_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, Posterior_salivary_gland_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, Retina_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, Skin_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, Subesophageal_brain_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, Suckers_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, Supraesophageal_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, Testes_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, Viscera_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, B5_L1_1_abundance_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, B5_L1_2_abundance_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, B5_L1_3_abundance_TPM$target_id)
identical(Axial_nerve_cord_TPM$target_id, OBJ1_abundance_TPM$target_id)
##########################################################################################

# atttach all the column from the pubblic analysis
##########################################################################################
All_TPM_tissue_our <- cbind(Axial_nerve_cord_TPM$tpm, Optic_lobe_TPM$tpm, Ova_TPM$tpm, Posterior_salivary_gland_TPM$tpm, Retina_TPM$tpm, Skin_TPM$tpm,
                            Subesophageal_brain_TPM$tpm, Suckers_TPM$tpm, Supraesophageal_TPM$tpm, Testes_TPM$tpm, Viscera_TPM$tpm,
                            B5_L1_1_abundance_TPM$tpm, B5_L1_2_abundance_TPM$tpm, B5_L1_3_abundance_TPM$tpm, OBJ1_abundance_TPM$tpm)
colnames(All_TPM_tissue_our) <- c("Axial_nerve_cord_TPM", "Optic_lobe_TPM", "Ova_TPM", "Posterior_salivary_gland_TPM", "Retina_TPM", 
                                  "Skin_TPM", "Subesophageal_brain_TPM", "Suckers_TPM", "Supraesophageal_TPM", "Testes_TPM", "Viscera_TPM",
                                  "Arm_B5_L1_1_TPM", "Arm_B5_L1_2_TPM", "Arm_B5_L1_3_TPM", "Junior_OBJ1_TPM")

# for having a column with ID
target_id <- Axial_nerve_cord_TPM$target_id
All_TPM_tissue_our_GeneID <- cbind.data.frame(All_TPM_tissue_our, target_id)

save(All_TPM_tissue_our_GeneID, file="./R_Objects/All_TPM_tissue_our_GeneID.RData")
write.table(All_TPM_tissue_our_GeneID, file = "./Output_Data/All_TPM_tissue_our_GeneID.txt",
            quote=F, sep='\t', row.names = F, col.names = T)
##########################################################################################

# preparing dataset for using in the heatmap
##########################################################################################
rownames(All_TPM_tissue_our) <- target_id

# z score normalization
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# apply function and transpose matrix
All_TPM_tissue_our_norm <- t(apply(All_TPM_tissue_our, 1, cal_z_score))

# Removing row with NaN (Not a Number)
All_TPM_tissue_our_norm_filtered <- All_TPM_tissue_our_norm[-which(is.nan(All_TPM_tissue_our_norm)), ]
##########################################################################################

# creating annotation from the dendrogram and export the data of heatmap clustering 
##########################################################################################
library("pheatmap")
my_heatmap_out <-pheatmap(All_TPM_tissue_our_norm_filtered, silent = TRUE) # with silent to do NOT plot the heatmap

#Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
my_heatmap_orderedGenes <- rownames(All_TPM_tissue_our_norm_filtered[my_heatmap_out$tree_row[["order"]],])

#Re-order original data (samples) to match ordering in heatmap (left-to-right)
my_heatmap_orderedSamples <- colnames(All_TPM_tissue_our_norm_filtered[,my_heatmap_out$tree_col[["order"]]])

# Plot Dendrogram (in case many row/genes it is not useful)
plot(my_heatmap_out$tree_row)
plot(my_heatmap_out$tree_col)

# If you want something like gene-to-cluster assignment, you can 'cut' 
# your row dendrogram into a pre-selected number of groups as follows
# 13 groups sorted by cluster numbers
my_heatmap_orderedGenes_cluster <- sort(cutree(my_heatmap_out$tree_row, k=13))
my_heatmap_orderedGenes_cluster <- cutree(my_heatmap_out$tree_row, k=13)

# count number of genes or extract gene list for each cluster
length(which(my_heatmap_orderedGenes_cluster=="1")) # cluster 1: 1840 genes
my_heatmap_orderedGenes_cluster[my_heatmap_orderedGenes_cluster==1]

# convert into data.frame for using in pheatmap
my_heatmap_orderedGenes_cluster_df <- data.frame(my_heatmap_orderedGenes_cluster)
colnames(my_heatmap_orderedGenes_cluster_df) <- c("cluster")

# add a word in front of cluster number so that colors in heatmap are not a progressive gradient 
my_heatmap_orderedGenes_cluster_df$cluster <- paste0("cluster", my_heatmap_orderedGenes_cluster_df$cluster)
##########################################################################################

# plot the heatmap with the annotation from the dendrogram clustering
##########################################################################################
library("pheatmap")
my_cluster_colour <- list(cluster = c(cluster1 = "#a6cee3", cluster2 = "#1f78b4", cluster3 = "#b2df8a",
                                      cluster4 = "#33a02c", cluster5 = "#fb9a99", cluster6 = "#e31a1c",
                                      cluster7 = "#fdbf6f", cluster8 = "#ff7f00", cluster9 = "#cab2d6",
                                      cluster10 = "#6a3d9a", cluster11 = "#FFDC15", cluster12 = "#b15928",
                                      cluster13 = "#808184"))

my_heatmap <- pheatmap(All_TPM_tissue_our_norm_filtered, 
                       annotation_row = my_heatmap_orderedGenes_cluster_df, 
                       annotation_colors = my_cluster_colour, 
                       cutree_rows = 13, show_rownames=F, cellwidth=15)

# creating function to save the object of the heatmap so don't need to run again for save it
save_pheatmap_png <- function(x, filename, width= 10, height= 7.5, units="in", res=300) {
  png(filename, width = width, height = height, units = units, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_eps <- function(x, filename, horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5) {
  postscript(filename, horizontal = horizontal, onefile = onefile, paper = paper, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "./Results/Heatmap_public_plusOur/plusJunior/RNA-seq_AllpublicTPM_OcBimac_ClustAll_plusOurAll_plusJunior_RNA-seq_slim15_Annot.png")
save_pheatmap_eps(my_heatmap, "./Results/Heatmap_public_plusOur/plusJunior/RNA-seq_AllpublicTPM_OcBimac_ClustAll_plusOurAll_plusJunior_RNA-seq_slim15_Annot.eps")
##########################################################################################

# add column with Gene_ID for exporting clusters
##########################################################################################
my_heatmap_orderedGenes_cluster_df_ID <- cbind(my_heatmap_orderedGenes_cluster_df,
                                               row.names(my_heatmap_orderedGenes_cluster_df) )
colnames(my_heatmap_orderedGenes_cluster_df_ID) <- c("cluster", "transcript_ID")
write.table(my_heatmap_orderedGenes_cluster_df_ID, file = "./Output_Data/Heatmap_orderedGenes_cluster_Gene_ID.txt",
            quote=F, sep='\t', row.names = F, col.names = T)

save(my_heatmap_orderedGenes_cluster_df_ID, file = "./R_Objects/Heatmap_orderedGenes_cluster_Gene_ID.RData")
##########################################################################################



##########################################################################################
### Figure S2 ###
##########################################################################################

# import all file *tsv with TPMs
##########################################################################################
OB2_19_L1_1_1_S1_TPM <- read.table("./Input_Data/TPM/OB2-19_L1_1_1_S1_abundance.tsv", sep="\t", header = T)
OB9_19_L1_1_4_S9_TPM <- read.table("./Input_Data/TPM/0B9-19_L1_1_4_S9_abundance.tsv", sep="\t", header = T)

# check if the each *tsv file have the gene/trascript ID in exact in same order
identical(OB2_19_L1_1_1_S1_TPM$target_id, OB9_19_L1_1_4_S9_TPM$target_id)
identical(OB2_19_L1_1_1_S1_TPM$target_id, All_TPM_tissue_our_GeneID$target_id)

ArmTips_L1_comparison <- cbind(OB2_19_L1_1_1_S1_TPM$tpm,
                               OB9_19_L1_1_4_S9_TPM$tpm,
                               All_TPM_tissue_our_GeneID$Arm_B5_L1_1_TPM)

colnames(ArmTips_L1_comparison) <- c("OB2_19_L1_1_1_S1_TPM",
                                     "OB9_19_L1_1_4_S9_TPM",
                                     "Arm_B5_L1_1_TPM")

View(ArmTips_L1_comparison)

# for having a column with ID
target_id <- OB2_19_L1_1_1_S1_TPM$target_id
ArmTips_L1_comparison_GeneID <- cbind.data.frame(ArmTips_L1_comparison, target_id)
View(ArmTips_L1_comparison_GeneID)

save(ArmTips_L1_comparison_GeneID, file="./R_Objects/ArmTips_L1_comparison_GeneID.RData")
write.table(ArmTips_L1_comparison_GeneID, file = "./Output_Data/ArmTips_L1_comparison_GeneID.txt",
            quote=F, sep='\t', row.names = F, col.names = T)
##########################################################################################

# extract top 1000 transcripts from each sample
##########################################################################################
library("dplyr")
OB2_19_L1_GeneID_top1000 <- top_n(ArmTips_L1_comparison_GeneID, 1000, OB2_19_L1_1_1_S1_TPM)
OB9_19_L1_GeneID_top1000 <- top_n(ArmTips_L1_comparison_GeneID, 1000, OB9_19_L1_1_4_S9_TPM)
Arm_B5_L1_1_TPM_GeneID_top1000 <- top_n(ArmTips_L1_comparison_GeneID, 1000, Arm_B5_L1_1_TPM)

# export for analyzing in Venny app (https://bioinfogp.cnb.csic.es/tools/venny/) and scaling it in Eulerr app (https://eulerr.co)
write.table(OB2_19_L1_GeneID_top1000, file = "./Output_Data/OB2_19_L1_GeneID_top1000.txt",
            quote=F, sep='\t', row.names = F, col.names = T)
write.table(OB9_19_L1_GeneID_top1000, file = "./Output_Data/OB9_19_L1_GeneID_top1000.txt",
            quote=F, sep='\t', row.names = F, col.names = T)
write.table(Arm_B5_L1_1_TPM_GeneID_top1000, file = "./Output_Data/Arm_B5_L1_1_TPM_GeneID_top1000.txt",
            quote=F, sep='\t', row.names = F, col.names = T)
##########################################################################################

# heatmap on all the transcripts used in Venn (1405 unique IDs among 3 samples)
##########################################################################################
AllTranscriptsVenn <- read.table("./Input_Data/Venn/AllTranscriptsVenn.txt", sep="\t", header = T)

# merging with sort=FALSE keep the order of X list
ArmTips_L1_comparison_GeneID_AllTranscriptsVenn <- merge(AllTranscriptsVenn, ArmTips_L1_comparison_GeneID, by = "target_id", sort = FALSE)

# for using in the heatmap
rownames(ArmTips_L1_comparison_GeneID_AllTranscriptsVenn) <- ArmTips_L1_comparison_GeneID_AllTranscriptsVenn$target_id
ArmTips_L1_comparison_GeneID_AllTranscriptsVenn_select <- subset(ArmTips_L1_comparison_GeneID_AllTranscriptsVenn, select = -c(target_id) )
ArmTips_L1_comparison_GeneID_AllTranscriptsVenn_select_matrix <- as.matrix(ArmTips_L1_comparison_GeneID_AllTranscriptsVenn_select)

library("pheatmap")
# Modify ordering of the clusters using clustering callback option
# with - in from of sv it reverse the ordering weight as decrescent 
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = -sv)
  as.hclust(dend)
}

png("./Results/pheatmap/TPM/RNA-seq_ArmTips_L1_comparison_AllTranscriptsVenn_ClustAll_log2_reordered2.png", width= 5, height= 5, units="in", res=300)
pheatmap(log2(ArmTips_L1_comparison_GeneID_AllTranscriptsVenn_select_matrix+1), show_rownames=F, cellwidth=15, clustering_callback = callback)
dev.off()
##########################################################################################



##########################################################################################
### Figure 1B and S3 ###
##########################################################################################

# import annotation from BioMart file
# with this columns Gene stable ID, Transcript stable ID, Protein stable ID, UniProtKB/TrEMBL ID, GO term accession, GO term name, GO domain
##########################################################################################
Annotation_BioMartGO <- read.delim("Input_Data/GeneAnnotation/BioMart/mart_export_UniProt_GO.txt")

# split the table using GO domain as 3 parent ID
GO_molecular_function <- subset(Annotation_BioMartGO, GO.domain==("molecular_function"))
GO_cellular_component <- subset(Annotation_BioMartGO, GO.domain==("cellular_component"))
GO_biological_process <- subset(Annotation_BioMartGO, GO.domain==("biological_process"))
##########################################################################################

# subsetting cluster extracted from the heatmap
##########################################################################################
heatmap_cluster1 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster1"))
heatmap_cluster2 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster2"))
heatmap_cluster3 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster3"))
heatmap_cluster4 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster4"))
heatmap_cluster5 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster5"))
heatmap_cluster6 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster6"))
heatmap_cluster7 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster7"))
heatmap_cluster8 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster8"))
heatmap_cluster9 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster9"))
heatmap_cluster10 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster10"))
heatmap_cluster11 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster11"))
heatmap_cluster12 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster12"))
heatmap_cluster13 <- subset(my_heatmap_orderedGenes_cluster_df_ID, cluster==("cluster13"))

# print number of genes for cluster
lengths(heatmap_cluster1)
lengths(heatmap_cluster2)
lengths(heatmap_cluster3)
lengths(heatmap_cluster4)
lengths(heatmap_cluster5)
lengths(heatmap_cluster6)
lengths(heatmap_cluster7)
lengths(heatmap_cluster8)
lengths(heatmap_cluster9)
lengths(heatmap_cluster10)
lengths(heatmap_cluster11)
lengths(heatmap_cluster12)
lengths(heatmap_cluster13)

# saving 
dir.create("Output_Data/heatmap_clusters")
write.table(heatmap_cluster1, file="./Output_Data/heatmap_clusters/heatmap_cluster1.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster2, file="./Output_Data/heatmap_clusters/heatmap_cluster2.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster3, file="./Output_Data/heatmap_clusters/heatmap_cluster3.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster4, file="./Output_Data/heatmap_clusters/heatmap_cluster4.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster5, file="./Output_Data/heatmap_clusters/heatmap_cluster5.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster6, file="./Output_Data/heatmap_clusters/heatmap_cluster6.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster7, file="./Output_Data/heatmap_clusters/heatmap_cluster7.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster8, file="./Output_Data/heatmap_clusters/heatmap_cluster8.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster9, file="./Output_Data/heatmap_clusters/heatmap_cluster9.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster10, file="./Output_Data/heatmap_clusters/heatmap_cluster10.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster11, file="./Output_Data/heatmap_clusters/heatmap_cluster11.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster12, file="./Output_Data/heatmap_clusters/heatmap_cluster12.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(heatmap_cluster13, file="./Output_Data/heatmap_clusters/heatmap_cluster13.txt", quote = FALSE, sep="\t", row.names = FALSE)
##########################################################################################

# preparation of dataframe for Go analysis in clusterProfiler without a model organism
##########################################################################################
library("dplyr")
term2gene_MolFun <- GO_molecular_function %>% select(GO.term.accession, Transcript.stable.ID) #TERM2GENE
term2name_MolFun <- GO_molecular_function %>% select(GO.term.accession, GO.term.name) #TERM2NAME

term2gene_CellComp <- GO_cellular_component %>% select(GO.term.accession, Transcript.stable.ID) #TERM2GENE
term2name_CellComp <- GO_cellular_component %>% select(GO.term.accession, GO.term.name) #TERM2NAME

term2gene_BioProc <- GO_biological_process %>% select(GO.term.accession, Transcript.stable.ID) #TERM2GENE
term2name_BioProc <- GO_biological_process %>% select(GO.term.accession, GO.term.name) #TERM2NAME
##########################################################################################

# Universal enrichment analysis by GO over-representation test in enricher function
##########################################################################################
library("clusterProfiler")
# GO_molecular_function
enricher_MolFun_cluster1 <- enricher(heatmap_cluster1$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster2 <- enricher(heatmap_cluster2$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster3 <- enricher(heatmap_cluster3$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster4 <- enricher(heatmap_cluster4$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster5 <- enricher(heatmap_cluster5$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster6 <- enricher(heatmap_cluster6$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster7 <- enricher(heatmap_cluster7$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster8 <- enricher(heatmap_cluster8$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster9 <- enricher(heatmap_cluster9$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster10 <- enricher(heatmap_cluster10$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster11 <- enricher(heatmap_cluster11$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster12 <- enricher(heatmap_cluster12$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
enricher_MolFun_cluster13 <- enricher(heatmap_cluster13$transcript_ID, TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)

# GO_cellular_component
enricher_CellComp_cluster1 <- enricher(heatmap_cluster1$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster2 <- enricher(heatmap_cluster2$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster3 <- enricher(heatmap_cluster3$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster4 <- enricher(heatmap_cluster4$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster5 <- enricher(heatmap_cluster5$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster6 <- enricher(heatmap_cluster6$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster7 <- enricher(heatmap_cluster7$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster8 <- enricher(heatmap_cluster8$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster9 <- enricher(heatmap_cluster9$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster10 <- enricher(heatmap_cluster10$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster11 <- enricher(heatmap_cluster11$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster12 <- enricher(heatmap_cluster12$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
enricher_CellComp_cluster13 <- enricher(heatmap_cluster13$transcript_ID, TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)

# GO_biological_process
enricher_BioProc_cluster1 <- enricher(heatmap_cluster1$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster2 <- enricher(heatmap_cluster2$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster3 <- enricher(heatmap_cluster3$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster4 <- enricher(heatmap_cluster4$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster5 <- enricher(heatmap_cluster5$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster6 <- enricher(heatmap_cluster6$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster7 <- enricher(heatmap_cluster7$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster8 <- enricher(heatmap_cluster8$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster9 <- enricher(heatmap_cluster9$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster10 <- enricher(heatmap_cluster10$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster11 <- enricher(heatmap_cluster11$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster12 <- enricher(heatmap_cluster12$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
enricher_BioProc_cluster13 <- enricher(heatmap_cluster13$transcript_ID, TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
##########################################################################################


# creating a unique plot for each GO.domain with all the clusters
##########################################################################################

# create a list with all transcript IDs
geneCluster_list = list(cluster1=heatmap_cluster1$transcript_ID,
                        cluster2=heatmap_cluster2$transcript_ID,
                        cluster3=heatmap_cluster3$transcript_ID,
                        cluster4=heatmap_cluster4$transcript_ID,
                        cluster5=heatmap_cluster5$transcript_ID,
                        cluster6=heatmap_cluster6$transcript_ID,
                        cluster7=heatmap_cluster7$transcript_ID,
                        cluster8=heatmap_cluster8$transcript_ID,
                        cluster9=heatmap_cluster9$transcript_ID,
                        cluster10=heatmap_cluster10$transcript_ID,
                        cluster11=heatmap_cluster11$transcript_ID,
                        cluster12=heatmap_cluster12$transcript_ID,
                        cluster13=heatmap_cluster13$transcript_ID)

# GO_molecular_function
#######################################
library("clusterProfiler")
compareCluster_MolFun <- compareCluster(geneCluster = geneCluster_list, fun = "enricher", TERM2GENE = term2gene_MolFun, TERM2NAME = term2name_MolFun)
head(compareCluster_MolFun)
tail(compareCluster_MolFun)

library("enrichplot")
library("ggplot2")
postscript("./Results/clusterProfiler/compareCluster_MolFun_5top_scaled.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 20, height = 15)
# dotplot has default 5 top categories per cluster or use showCategory=8
dotplot(compareCluster_MolFun) +
  geom_text(aes(label=Count, size=0.002, fontface='bold'), colour='white', show.legend = F) +
  scale_size(range = c(3, 9)) +
  theme(legend.key.height = unit(30,'points'),
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=18, angle=45, hjust=1),
        text = element_text(size=18)) + 
  ggtitle("GO_molecular_function for all clusters")
dev.off() # change in postscript width, height for a good aspect ratio
#######################################

# GO_cellular_component
#######################################
library("clusterProfiler")
compareCluster_CellComp <- compareCluster(geneCluster = geneCluster_list, fun = "enricher", TERM2GENE = term2gene_CellComp, TERM2NAME = term2name_CellComp)
head(compareCluster_CellComp)
tail(compareCluster_CellComp)

library("enrichplot")
library("ggplot2")
postscript("./Results/clusterProfiler/compareCluster_CellComp_5top_scaled.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10.5, height = 10)
dotplot(compareCluster_CellComp) +
  geom_text(aes(label=Count, size=0.002, fontface='bold'), colour='white', show.legend = F) +
  scale_size(range = c(3, 9)) +
  theme(legend.key.height = unit(30,'points'),
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=18, angle=45, hjust=1),
        text = element_text(size=18)) + 
  ggtitle("GO_cellular_component for all clusters")
dev.off() # change in postscript width, height for a good aspect ratio
#######################################

# GO_biological_process
#######################################
library("clusterProfiler")
compareCluster_BioProc <- compareCluster(geneCluster = geneCluster_list, fun = "enricher", TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
head(compareCluster_BioProc)
tail(compareCluster_BioProc)

library("enrichplot")
library("ggplot2")
postscript("./Results/clusterProfiler/compareCluster_BioProc_5top_scaled.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 15.5, height = 17)
dotplot(compareCluster_BioProc) +
  geom_text(aes(label=Count, size=0.002, fontface='bold'), colour='white', show.legend = F) +
  scale_size(range = c(3, 9)) +
  theme(legend.key.height = unit(30,'points'),
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=18, angle=45, hjust=1),
        text = element_text(size=18)) + 
  ggtitle("GO_biological_process for all clusters")
dev.off() # change in postscript width, height for a good aspect ratio
##########################################################################################
