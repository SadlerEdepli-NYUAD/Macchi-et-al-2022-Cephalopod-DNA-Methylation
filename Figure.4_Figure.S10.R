##########################################################################################
### Figure 4A-B - S10 ###
##########################################################################################

# generate annotation for full-lenght transcripts with ensembl ID as rownames
##########################################################################################
library("GenomicFeatures")
# need to run once re-open project because txdb is external linked file
gff3.file = "./Input_Data/Obimaculoides_280.gene.gff3"
txdb <- makeTxDbFromGFF(file=gff3.file, organism="Octopus bimaculoides")

# extract and use the transcript names as rownames from txdb
###################################################
GR_transcripts_OB_rownames <- transcripts(txdb, use.names=TRUE)

# matrixCalculation for EnrichedHeatmap on all transcripts (using transcripts name as rownames)
########################################################################
library("EnrichedHeatmap")
# using k=1 we generate a mean value of DNA methylation for the entire gene body of each transcripts so we can use to calcuate the kmeans distribution
normMat_WGBS_Brain_SupraE_mean = normalizeToMatrix(GR_Octopus_WGBS_Brain_SupraE_select, GR_transcripts_OB_rownames, value_column = "perc.meth", mean_mode = "absolute",
                                                   extend = 0, k = 1, w = 50)

# changing the target ratio we obtain a different proportional draw size of extended window
normMat_WGBS_Brain_SupraE_Exteded = normalizeToMatrix(GR_Octopus_WGBS_Brain_SupraE_select, GR_transcripts_OB_rownames, value_column = "perc.meth", mean_mode = "absolute",
                                                      extend = 2000, k = 50, w = 50, target_ratio= 0.625, background = NA)

normMat_WGBS_Brain_SubE_Exteded = normalizeToMatrix(GR_Octopus_WGBS_Brain_SubE_select, GR_transcripts_OB_rownames, value_column = "perc.meth", mean_mode = "absolute",
                                                    extend = 2000, k = 50, w = 50, target_ratio= 0.625, background = NA)

normMat_CpGs_All_GR_transcripts_OB_Exteded = normalizeToMatrix(GR_CpGs_All, GR_transcripts_OB_rownames, value_column = "perc.meth", mean_mode = "absolute",
                                                               extend = 2000, k = 50, w = 50, target_ratio= 0.625, background = NA)

# order by rownames so that we have all the matrix the same before merging them in the ht list
# need always make sure the original row orders of the matrices are the same (before clustering)
normMat_WGBS_Brain_SupraE_mean <- as.matrix(normMat_WGBS_Brain_SupraE_mean[order(rownames(normMat_WGBS_Brain_SupraE_mean)), ])
class(normMat_WGBS_Brain_SupraE_mean)

normMat_WGBS_Brain_SupraE_Exteded <- as.matrix(normMat_WGBS_Brain_SupraE_Exteded[order(rownames(normMat_WGBS_Brain_SupraE_Exteded)), ])
class(normMat_WGBS_Brain_SupraE_Exteded)

normMat_WGBS_Brain_SubE_Exteded <- as.matrix(normMat_WGBS_Brain_SubE_Exteded[order(rownames(normMat_WGBS_Brain_SubE_Exteded)), ])
class(normMat_WGBS_Brain_SubE_Exteded)

normMat_CpGs_All_GR_transcripts_OB_Exteded <- as.matrix(normMat_CpGs_All_GR_transcripts_OB_Exteded[order(rownames(normMat_CpGs_All_GR_transcripts_OB_Exteded)), ])
class(normMat_CpGs_All_GR_transcripts_OB_Exteded)

# check for consistency if TRUE are the same name and order
identical(rownames(normMat_WGBS_Brain_SupraE_mean), rownames(normMat_WGBS_Brain_SupraE_Exteded))
identical(rownames(normMat_WGBS_Brain_SupraE_mean), rownames(normMat_WGBS_Brain_SubE_Exteded))
identical(rownames(normMat_WGBS_Brain_SupraE_mean), rownames(normMat_CpGs_All_GR_transcripts_OB_Exteded))

save(normMat_WGBS_Brain_SupraE_mean, file="./R_Objects/normMat_WGBS_Brain_SupraE_mean.RData")
save(normMat_WGBS_Brain_SupraE_Exteded, file="./R_Objects/normMat_WGBS_Brain_SupraE_Exteded.RData")
save(normMat_WGBS_Brain_SubE_Exteded, file="./R_Objects/normMat_WGBS_Brain_SubE_Exteded.RData")
save(normMat_CpGs_All_GR_transcripts_OB_Exteded, file="./R_Objects/normMat_CpGs_All_GR_transcripts_OB_Exteded.RData")
########################################################################

# use kmeans a first time to get the centers based on mean values of each gene body (Fig. 4A-B and S10A-B)
########################################################################
centers <- kmeans(normMat_WGBS_Brain_SupraE_mean, centers = 4)$centers

# order the centers, run using stats 4.0.2
# since the ordering of kmeans is randomly initialized we ordered them based on perc.meth values to have high meth as first cluster and low meth as fourth one
centers <- sort(centers, decreasing = TRUE)

# transcript ordered in clusters, run using stats 4.0.2
# extract ordered transcripts in each meth.cluster and create data.frame for later use
clusters <- kmeans(normMat_WGBS_Brain_SupraE_mean, centers = centers)$cluster
df_clusters <- as.data.frame(clusters)
colnames(df_clusters) <- "cluster"
df_clusters <- cbind(df_clusters, transcript_ID=rownames(df_clusters))
df_clusters$cluster <- paste0("cluster", df_clusters$cluster)

# call kmeans again but this time passing the centers calculated in the previous step to generate partition for heatmap
partition = paste0("cluster", kmeans(normMat_WGBS_Brain_SupraE_mean, centers = centers)$cluster)
lgd = Legend(at = c("cluster1", "cluster2", "cluster3", "cluster4"), title = "Clusters", 
             type = "lines", legend_gp = gpar(col = 2:5))
ht_list = Heatmap(partition, col = structure(2:5, names = paste0("cluster", 1:4)), name = "partition",
                  show_row_names = FALSE, width = unit(3, "mm")) +
  EnrichedHeatmap(normMat_WGBS_Brain_SupraE_Exteded, name = "methylation_SupraE",
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5), ylim = c(0, 100))), 
                  column_title = "Methylation") +
  EnrichedHeatmap(normMat_WGBS_Brain_SubE_Exteded, name = "methylation_SubE",
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5), ylim = c(0, 100))), 
                  column_title = "Methylation") +
  EnrichedHeatmap(normMat_CpGs_All_GR_transcripts_OB_Exteded, name = "methylation_Hatchlings",
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5), ylim = c(0, 100))), 
                  column_title = "Methylation") 

pdf("./Results/EnrichedHeatmap/Meth_Expres/WGBS_Brain_SupraE_Sub_Hatch_Exteded0.625_Meth.Cluster_stats_4.0.2_new.pdf", onefile = FALSE, paper = "special", width = 20, height = 7.5)
draw(ht_list, split = partition, annotation_legend_list = list(lgd), 
     ht_gap = unit(c(2, 8, 8), "mm"))
dev.off()

# count transcripts in meth.clusters using the partition generated for heatmap
partition_df <- data.frame(partition)
library("dplyr")
cluster_count <- partition_df %>%
  group_by(partition) %>% count()
########################################################################

# extracting the rownames with transcript ID to get the transcript in each cluster
########################################################################
meth.cluster1 <- subset(df_clusters, cluster=="cluster1")
meth.cluster2 <- subset(df_clusters, cluster=="cluster2")
meth.cluster3 <- subset(df_clusters, cluster=="cluster3")
meth.cluster4 <- subset(df_clusters, cluster=="cluster4")

dir.create("./Output_Data/Meth_Cluster")
write.table(meth.cluster1, file="./Output_Data/Meth_Cluster/meth.cluster1.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(meth.cluster2, file="./Output_Data/Meth_Cluster/meth.cluster2.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(meth.cluster3, file="./Output_Data/Meth_Cluster/meth.cluster3.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(meth.cluster4, file="./Output_Data/Meth_Cluster/meth.cluster4.txt", quote = FALSE, sep="\t", row.names = FALSE)

# subsetting extract tx_name by transcript_ID for each cluster to get 
# transcript position, width and expression TPM for the tissues SupraE=Supraesophageal_TPM, SubE=Subesophageal_brain_TPM, Hatchling=Junior_OBJ1_TPM
########################################################################
All_TPM_tissue_our_GeneID_commonNA <- merge(GR_transcripts_OB_rownames, All_TPM_tissue_our_GeneID, 
                                            by.x="tx_name", by.y="target_id",
                                            sort = TRUE, all.x=TRUE)
rownames(All_TPM_tissue_our_GeneID_commonNA) <- All_TPM_tissue_our_GeneID_commonNA$tx_name

df_meth.cluster1 <- subset(All_TPM_tissue_our_GeneID_commonNA, tx_name %in% meth.cluster1$transcript_ID, select=c(seqnames, start, end, width, strand, tx_id, tx_name, Supraesophageal_TPM, Subesophageal_brain_TPM, Junior_OBJ1_TPM))
GR_meth.cluster1 <- as(df_meth.cluster1, "GRanges")

df_meth.cluster2 <- subset(All_TPM_tissue_our_GeneID_commonNA, tx_name %in% meth.cluster2$transcript_ID, select=c(seqnames, start, end, width, strand, tx_id, tx_name, Supraesophageal_TPM, Subesophageal_brain_TPM, Junior_OBJ1_TPM))
GR_meth.cluster2 <- as(df_meth.cluster2, "GRanges")

df_meth.cluster3 <- subset(All_TPM_tissue_our_GeneID_commonNA, tx_name %in% meth.cluster3$transcript_ID, select=c(seqnames, start, end, width, strand, tx_id, tx_name, Supraesophageal_TPM, Subesophageal_brain_TPM, Junior_OBJ1_TPM))
GR_meth.cluster3 <- as(df_meth.cluster3, "GRanges")

df_meth.cluster4 <- subset(All_TPM_tissue_our_GeneID_commonNA, tx_name %in% meth.cluster4$transcript_ID, select=c(seqnames, start, end, width, strand, tx_id, tx_name, Supraesophageal_TPM, Subesophageal_brain_TPM, Junior_OBJ1_TPM))
GR_meth.cluster4 <- as(df_meth.cluster4, "GRanges")

write.table(df_meth.cluster1, file="./Output_Data/Meth_Cluster/df_meth.cluster1.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(df_meth.cluster2, file="./Output_Data/Meth_Cluster/df_meth.cluster2.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(df_meth.cluster3, file="./Output_Data/Meth_Cluster/df_meth.cluster3.txt", quote = FALSE, sep="\t", row.names = FALSE)
write.table(df_meth.cluster4, file="./Output_Data/Meth_Cluster/df_meth.cluster4.txt", quote = FALSE, sep="\t", row.names = FALSE)
########################################################################
##########################################################################################

# violin plot on SupraE_TPM expression of transcript divided by meth.clusters (Fig.4C)
##########################################################################################
# create the list that are different identical with NAs
sq <- seq(max(length(df_meth.cluster1$Supraesophageal_TPM),
              length(df_meth.cluster2$Supraesophageal_TPM),
              length(df_meth.cluster3$Supraesophageal_TPM),
              length(df_meth.cluster4$Supraesophageal_TPM)
))

# this is forcing to keep the same lenght for all columns and put NAs
df_meth.cluster_SupraE_TPM <- data.frame(df_meth.cluster1$Supraesophageal_TPM[sq],
                                         df_meth.cluster2$Supraesophageal_TPM[sq],
                                         df_meth.cluster3$Supraesophageal_TPM[sq],
                                         df_meth.cluster4$Supraesophageal_TPM[sq]
)

head(df_meth.cluster_SupraE_TPM)
colnames(df_meth.cluster_SupraE_TPM) <- c("cluster1_transcript_SupraE_TPM",
                                          "cluster2_transcript_SupraE_TPM",
                                          "cluster3_transcript_SupraE_TPM",
                                          "cluster4_transcript_SupraE_TPM"
)

View(df_meth.cluster_SupraE_TPM)
dir.create("Analysis")
write.table(df_meth.cluster_SupraE_TPM, file="./Analysis/df_meth.cluster_SupraE_TPM.txt", quote = FALSE, sep="\t", row.names = FALSE)

library("reshape2")
df_meth.cluster_SupraE_TPM_melt <- melt(df_meth.cluster_SupraE_TPM)
head(df_meth.cluster_SupraE_TPM_melt)
tail(df_meth.cluster_SupraE_TPM_melt)

library("ggplot2")
dir.create("./Results/ggplots")
dir.create("./Results/ggplots/Meth.Cluster_transcript_TPMs")

postscript("./Results/ggplots/Meth.Cluster_transcript_TPMs/cluster_transcript_SupraE_TPM_ViolinPlot_Area.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
ggplot(df_meth.cluster_SupraE_TPM_melt, aes(x = variable, y = log2(value+1), fill=variable)) + 
  geom_violin(trim = TRUE, scale = "area") + #scale_y_continuous(breaks=seq(0,100,25)) +
  labs(x= "Meth.Cluster_transcript_TPMs", y = "log2(transcript_SupraE_TPM)") + 
  theme_classic() + scale_fill_manual(values = c("#DF536B","#61D04F","#2297E6","#28E2E5")) + 
  theme(axis.text.x = element_text(angle=30, hjust = 1)) +
  geom_boxplot(width=0.05, fill="white", notch = TRUE)
dev.off()
########################################################################

# violin plot on transcript-width of transcript divided by meth.clusters (Fig.S10D)
########################################################################
# create the list that are different identical with NAs
sq <- seq(max(length(df_meth.cluster1$width),
              length(df_meth.cluster2$width),
              length(df_meth.cluster3$width),
              length(df_meth.cluster4$width)
))

# this is forcing to keep the same lenght for all columns and put NAs
df_meth.cluster_width <- data.frame(df_meth.cluster1$width[sq],
                                    df_meth.cluster2$width[sq],
                                    df_meth.cluster3$width[sq],
                                    df_meth.cluster4$width[sq]
)

head(df_meth.cluster_width)
colnames(df_meth.cluster_width) <- c("cluster1_transcript_width",
                                     "cluster2_transcript_width",
                                     "cluster3_transcript_width",
                                     "cluster4_transcript_width"
)

View(df_meth.cluster_width)
dir.create("Analysis")
write.table(df_meth.cluster_width, file="./Analysis/df_meth.cluster_width.txt", quote = FALSE, sep="\t", row.names = FALSE)

library("reshape2")
df_meth.cluster_width_melt <- melt(df_meth.cluster_width)
head(df_meth.cluster_width_melt)
tail(df_meth.cluster_width_melt)

library("ggplot2")
dir.create("./Results/ggplots")
dir.create("./Results/ggplots/Meth.Cluster_transcript_width")

postscript("./Results/ggplots/Meth.Cluster_transcript_width/cluster_transcript_width_ViolinPlot_Area.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
ggplot(df_meth.cluster_width_melt, aes(x = variable, y = log10(value), fill=variable)) + 
  geom_violin(trim = TRUE, scale = "area") + #scale_y_continuous(breaks=seq(0,100,25)) +
  labs(x= "Meth.Cluster_transcript_width", y = "log10(transcript_width)") + 
  theme_classic() + scale_fill_manual(values = c("#DF536B","#61D04F","#2297E6","#28E2E5")) + 
  theme(axis.text.x = element_text(angle=30, hjust = 1)) +
  geom_boxplot(width=0.05, fill="white", notch = TRUE)
dev.off()
########################################################################
##########################################################################################


# regression plot for methylation values of gene bodies and expression (Fig. 4D and S10C)
##########################################################################################

# Brain SupraEsophageal
########################################################################
library("GenomicRanges")
# merge by overlap is creating a dataset with both the input maintained in the output
WGBS_Brain_SupraE_transcripts_merged <- mergeByOverlaps(GR_Octopus_WGBS_Brain_SupraE_select, GR_transcripts_OB)
df_WGBS_Brain_SupraE_transcripts_merged <- data.frame(tx_name = WGBS_Brain_SupraE_transcripts_merged$tx_name, 
                                                      perc.meth = WGBS_Brain_SupraE_transcripts_merged$perc.meth)

# perform mean of CpGs methylation for each transcript
library("dplyr")
df_WGBS_Brain_SupraE_transcripts_perc.meth_mean <- df_WGBS_Brain_SupraE_transcripts_merged %>%
  group_by(tx_name) %>%
  summarise(perc.meth_mean = mean(perc.meth))

# extract expression data for only used dataset
load("./Output_Data/All_TPM_tissue_our_GeneID.RData")
All_TPM_tissue_our_GeneID

All_TPM_Supraesophageal <- data.frame(Supraesophageal_TPM = All_TPM_tissue_our_GeneID$Supraesophageal_TPM, 
                                      target_id = All_TPM_tissue_our_GeneID$target_id)

# merge methylation and expression
df_WGBS_Brain_SupraE_perc.meth_mean_Expres <- merge(df_WGBS_Brain_SupraE_transcripts_perc.meth_mean, 
                                                    All_TPM_Supraesophageal, by.x="tx_name", by.y="target_id")

# calculate percentile for TPM
df_WGBS_Brain_SupraE_perc.meth_mean_Expres_percentile <- df_WGBS_Brain_SupraE_perc.meth_mean_Expres %>%
  mutate(percentile = ntile(Supraesophageal_TPM,100))

# calculate mean of methylation and expression for each percentile block
df_WGBS_Brain_SupraE_perc.meth_mean_Expres_percentile_mean <- df_WGBS_Brain_SupraE_perc.meth_mean_Expres_percentile %>%
  group_by(percentile) %>%
  summarise_at(c("perc.meth_mean", "Supraesophageal_TPM"), mean, na.rm = TRUE)

# plot data
df_WGBS_Brain_SupraE_perc.meth_mean_Expres_percentile_mean_select <- data.frame(df_WGBS_Brain_SupraE_perc.meth_mean_Expres_percentile_mean$perc.meth_mean, 
                                                                                df_WGBS_Brain_SupraE_perc.meth_mean_Expres_percentile_mean$Supraesophageal_TPM)
head(df_WGBS_Brain_SupraE_perc.meth_mean_Expres_percentile_mean_select)
colnames(df_WGBS_Brain_SupraE_perc.meth_mean_Expres_percentile_mean_select) <- c("perc.meth","TPM")

library("ggplot2")
dir.create("./Results/ggplots")
dir.create("./Results/ggplots/Meth_Expres")

pdf("./Results/ggplots/Meth_Expres/WGBS_Brain_SupraE_perc.meth_mean_Expres_percentile_JitterPlot_blackDot.pdf", onefile = FALSE, paper = "special", width = 8.5, height = 7.5)
ggplot(df_WGBS_Brain_SupraE_perc.meth_mean_Expres_percentile_mean_select, aes(x = log2(TPM + 1), y = perc.meth)) + 
  geom_jitter() + geom_smooth(color="darkred", lwd=1.5) + 
  labs(x= "log2(TPM + 1)", y = "% of methylation") + 
  theme_classic()
dev.off()
########################################################################

# Brain SubEsophageal
########################################################################
WGBS_Brain_SubE_transcripts_merged <- mergeByOverlaps(GR_Octopus_WGBS_Brain_SubE_select, GR_transcripts_OB)
df_WGBS_Brain_SubE_transcripts_merged <- data.frame(tx_name = WGBS_Brain_SubE_transcripts_merged$tx_name, 
                                                    perc.meth = WGBS_Brain_SubE_transcripts_merged$perc.meth)

# perform mean of CpGs methylation for each transcript
library("dplyr")
df_WGBS_Brain_SubE_transcripts_perc.meth_mean <- df_WGBS_Brain_SubE_transcripts_merged %>%
  group_by(tx_name) %>%
  summarise(perc.meth_mean = mean(perc.meth))

# extract expression data for only used dataset
All_TPM_Subesophageal <- data.frame(Subesophageal_TPM = All_TPM_tissue_our_GeneID$Subesophageal_brain_TPM, 
                                    target_id = All_TPM_tissue_our_GeneID$target_id)

# merge methyalation and expression
df_WGBS_Brain_SubE_perc.meth_mean_Expres <- merge(df_WGBS_Brain_SubE_transcripts_perc.meth_mean, 
                                                  All_TPM_Subesophageal, by.x="tx_name", by.y="target_id")

# calculate percentile for TPM
df_WGBS_Brain_SubE_perc.meth_mean_Expres_percentile <- df_WGBS_Brain_SubE_perc.meth_mean_Expres %>%
  mutate(percentile = ntile(Subesophageal_TPM,100))

# calculate mean of methylation and expression for each percentile block
df_WGBS_Brain_SubE_perc.meth_mean_Expres_percentile_mean <- df_WGBS_Brain_SubE_perc.meth_mean_Expres_percentile %>%
  group_by(percentile) %>%
  summarise_at(c("perc.meth_mean", "Subesophageal_TPM"), mean, na.rm = TRUE)

# plot data
df_WGBS_Brain_SubE_perc.meth_mean_Expres_percentile_mean_select <- data.frame(df_WGBS_Brain_SubE_perc.meth_mean_Expres_percentile_mean$perc.meth_mean, 
                                                                              df_WGBS_Brain_SubE_perc.meth_mean_Expres_percentile_mean$Subesophageal_TPM)
head(df_WGBS_Brain_SubE_perc.meth_mean_Expres_percentile_mean_select)
colnames(df_WGBS_Brain_SubE_perc.meth_mean_Expres_percentile_mean_select) <- c("perc.meth","TPM")

library("ggplot2")
dir.create("./Results/ggplots")
dir.create("./Results/ggplots/Meth_Expres")

pdf("./Results/ggplots/Meth_Expres/WGBS_Brain_SubE_perc.meth_mean_Expres_percentile_JitterPlot_blackDot.pdf", onefile = FALSE, paper = "special", width = 8.5, height = 7.5)
ggplot(df_WGBS_Brain_SubE_perc.meth_mean_Expres_percentile_mean_select, aes(x = log2(TPM + 1), y = perc.meth)) + 
  geom_jitter() + geom_smooth(color="darkred", lwd=1.5) + 
  labs(x= "log2(TPM + 1)", y = "% of methylation") + 
  theme_classic()
dev.off()
########################################################################

# Hatchlings
########################################################################
GR_CpGs_All_transcripts_merged <- mergeByOverlaps(GR_CpGs_All, GR_transcripts_OB)
df_GR_CpGs_All_transcripts_merged <- data.frame(tx_name = GR_CpGs_All_transcripts_merged$tx_name, 
                                                perc.meth = GR_CpGs_All_transcripts_merged$perc.meth)

# perform mean of CpGs methylation for each transcript
library("dplyr")
df_GR_CpGs_All_transcripts_perc.meth_mean <- df_GR_CpGs_All_transcripts_merged %>%
  group_by(tx_name) %>%
  summarise(perc.meth_mean = mean(perc.meth))

countCpG_eachTranscript <- df_GR_CpGs_All_transcripts_merged %>%
  count(tx_name)

# extract expression data for only used dataset
load("/Users/fm1442/Sadler Edepli Lab Dropbox/Filippo Macchi/BioInformatics/06_Projects/Sadler_Lab/04_Octopus/Octopus_RNA-seq/Output_Data/All_TPM_tissue_our_GeneID.RData")
All_TPM_tissue_our_GeneID

All_TPM_Hatchlings <- data.frame(Hatchlings_TPM = All_TPM_tissue_our_GeneID$Junior_OBJ1_TPM, 
                                 target_id = All_TPM_tissue_our_GeneID$target_id)

# merge methyalation and expression
df_GR_CpGs_All_perc.meth_mean_Expres <- merge(df_GR_CpGs_All_transcripts_perc.meth_mean, 
                                              All_TPM_Hatchlings, by.x="tx_name", by.y="target_id")

# calculate percentile for TPM
df_GR_CpGs_All_perc.meth_mean_Expres_percentile <- df_GR_CpGs_All_perc.meth_mean_Expres %>%
  mutate(percentile = ntile(Hatchlings_TPM,100))

# calculate mean of methylation and expression for each percentile block
df_GR_CpGs_All_perc.meth_mean_Expres_percentile_mean <- df_GR_CpGs_All_perc.meth_mean_Expres_percentile %>%
  group_by(percentile) %>%
  summarise_at(c("perc.meth_mean", "Hatchlings_TPM"), mean, na.rm = TRUE)

# plot data
df_GR_CpGs_All_perc.meth_mean_Expres_percentile_mean_select <- data.frame(df_GR_CpGs_All_perc.meth_mean_Expres_percentile_mean$perc.meth_mean, 
                                                                          df_GR_CpGs_All_perc.meth_mean_Expres_percentile_mean$Hatchlings_TPM)
head(df_GR_CpGs_All_perc.meth_mean_Expres_percentile_mean_select)
colnames(df_GR_CpGs_All_perc.meth_mean_Expres_percentile_mean_select) <- c("perc.meth","TPM")

library("ggplot2")
dir.create("./Results/ggplots")
dir.create("./Results/ggplots/Meth_Expres")

pdf("./Results/ggplots/Meth_Expres/GR_CpGs_All_perc.meth_mean_Expres_percentile_JitterPlot_blackDot.pdf", onefile = FALSE, paper = "special", width = 8.5, height = 7.5)
ggplot(df_GR_CpGs_All_perc.meth_mean_Expres_percentile_mean_select, aes(x = log2(TPM + 1), y = perc.meth)) + 
  geom_jitter() + geom_smooth(color="darkred", lwd=1.5) + 
  labs(x= "log2(TPM + 1)", y = "% of methylation") + 
  theme_classic()
dev.off()
########################################################################
##########################################################################################

# heatmap on expression TPMs of all tissues divided by meth.cluster (Fig.4E)
##########################################################################################

# subsetting extract tx_name by transcript_ID for each cluster 
########################################################################
All_TPM_tissue_meth.cluster1 <- subset(All_TPM_tissue_our_GeneID_commonNA, 
                                       tx_name %in% meth.cluster1$transcript_ID, 
                                       select=c(Testes_TPM, Retina_TPM, Subesophageal_brain_TPM, 
                                                Supraesophageal_TPM, Axial_nerve_cord_TPM, 
                                                Optic_lobe_TPM, Arm_B5_L1_2_TPM, Arm_B5_L1_3_TPM, 
                                                Arm_B5_L1_1_TPM, Junior_OBJ1_TPM, Skin_TPM, Suckers_TPM, 
                                                Viscera_TPM, Ova_TPM, Posterior_salivary_gland_TPM))

All_TPM_tissue_meth.cluster2 <- subset(All_TPM_tissue_our_GeneID_commonNA, 
                                       tx_name %in% meth.cluster2$transcript_ID, 
                                       select=c(Testes_TPM, Retina_TPM, Subesophageal_brain_TPM, 
                                                Supraesophageal_TPM, Axial_nerve_cord_TPM, 
                                                Optic_lobe_TPM, Arm_B5_L1_2_TPM, Arm_B5_L1_3_TPM, 
                                                Arm_B5_L1_1_TPM, Junior_OBJ1_TPM, Skin_TPM, Suckers_TPM, 
                                                Viscera_TPM, Ova_TPM, Posterior_salivary_gland_TPM))

All_TPM_tissue_meth.cluster3 <- subset(All_TPM_tissue_our_GeneID_commonNA, 
                                       tx_name %in% meth.cluster3$transcript_ID, 
                                       select=c(Testes_TPM, Retina_TPM, Subesophageal_brain_TPM, 
                                                Supraesophageal_TPM, Axial_nerve_cord_TPM, 
                                                Optic_lobe_TPM, Arm_B5_L1_2_TPM, Arm_B5_L1_3_TPM, 
                                                Arm_B5_L1_1_TPM, Junior_OBJ1_TPM, Skin_TPM, Suckers_TPM, 
                                                Viscera_TPM, Ova_TPM, Posterior_salivary_gland_TPM))

All_TPM_tissue_meth.cluster4 <- subset(All_TPM_tissue_our_GeneID_commonNA, 
                                       tx_name %in% meth.cluster4$transcript_ID, 
                                       select=c(Testes_TPM, Retina_TPM, Subesophageal_brain_TPM, 
                                                Supraesophageal_TPM, Axial_nerve_cord_TPM, 
                                                Optic_lobe_TPM, Arm_B5_L1_2_TPM, Arm_B5_L1_3_TPM, 
                                                Arm_B5_L1_1_TPM, Junior_OBJ1_TPM, Skin_TPM, Suckers_TPM, 
                                                Viscera_TPM, Ova_TPM, Posterior_salivary_gland_TPM))
########################################################################

All_TPM_tissue_meth.cluster1_mat <- as.matrix(All_TPM_tissue_meth.cluster1)
All_TPM_tissue_meth.cluster2_mat <- as.matrix(All_TPM_tissue_meth.cluster2)
All_TPM_tissue_meth.cluster3_mat <- as.matrix(All_TPM_tissue_meth.cluster3)
All_TPM_tissue_meth.cluster4_mat <- as.matrix(All_TPM_tissue_meth.cluster4)


# Removing row with NAs where present
########################################################################
# in this case only meth.cluster4 has to begin with some NAs
# since there are transcripts not covered in our RNA-seq
All_TPM_tissue_meth.cluster4_mat_filtered <- All_TPM_tissue_meth.cluster4_mat[-which(is.na(All_TPM_tissue_meth.cluster4_mat)), ]
# removing only few transcripts we got 17,557 from 17,586

# log2 transformation +1 will generate a data normalized
All_TPM_tissue_meth.cluster1_mat_log2 <- log2(All_TPM_tissue_meth.cluster1_mat + 1)
All_TPM_tissue_meth.cluster2_mat_log2 <- log2(All_TPM_tissue_meth.cluster2_mat + 1)
All_TPM_tissue_meth.cluster3_mat_log2 <- log2(All_TPM_tissue_meth.cluster3_mat + 1)
All_TPM_tissue_meth.cluster4_mat_log2 <- log2(All_TPM_tissue_meth.cluster4_mat_filtered + 1)
########################################################################

# heatmapping
########################################################################
library("pheatmap")
dir.create("./Results/pheatmap")

# Modify ordering of the clusters using clustering callback option
# with - in from of sv it reverse the ordering weight as decrescent 
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = -sv)
  as.hclust(dend)
}

png("./Results/pheatmap/All_TPM_tissue_meth.cluster1_mat_log2_ordered.png", width= 5, height= 3.5, units="in", res=300)
pheatmap(All_TPM_tissue_meth.cluster1_mat_log2, cluster_cols=F, show_rownames=F, cellwidth=15, clustering_callback = callback)
dev.off()
png("./Results/pheatmap/All_TPM_tissue_meth.cluster2_mat_log2_ordered.png", width= 5, height= 3.5, units="in", res=300)
pheatmap(All_TPM_tissue_meth.cluster2_mat_log2, cluster_cols=F, show_rownames=F, cellwidth=15, clustering_callback = callback)
dev.off()
png("./Results/pheatmap/All_TPM_tissue_meth.cluster3_mat_log2_ordered.png", width= 5, height= 3.5, units="in", res=300)
pheatmap(All_TPM_tissue_meth.cluster3_mat_log2, cluster_cols=F, show_rownames=F, cellwidth=15, clustering_callback = callback)
dev.off()
png("./Results/pheatmap/All_TPM_tissue_meth.cluster4_mat_log2_ordered.png", width= 5, height= 3.5, units="in", res=300)
pheatmap(All_TPM_tissue_meth.cluster4_mat_log2, cluster_cols=F, show_rownames=F, cellwidth=15, clustering_callback = callback)
dev.off()
########################################################################
##########################################################################################


# Universal enrichment analysis by GO over-representation test in enricher function (Fig.4F)
##########################################################################################

# creating a unique plot for all the cluster
########################################################################
library("enrichplot")
# create a list with all transcript IDs
geneClusterMeth_list = list(meth.cluster1=meth.cluster1$transcript_ID,
                            meth.cluster2=meth.cluster2$transcript_ID,
                            meth.cluster3=meth.cluster3$transcript_ID,
                            meth.cluster4=meth.cluster4$transcript_ID
)

# GO_biological_process
#######################################
library("clusterProfiler")
compareClusterMeth_BioProc <- compareCluster(geneCluster = geneClusterMeth_list, fun = "enricher", TERM2GENE = term2gene_BioProc, TERM2NAME = term2name_BioProc)
#head(as.data.frame(compareClusterMeth_BioProc))
head(compareClusterMeth_BioProc)
tail(compareClusterMeth_BioProc)

library("enrichplot")
library("ggplot2")
postscript("./Results/clusterProfiler/ClusterMeth/compareClusterMeth_BioProc_15top_scaled.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 13, height = 19)
dotplot(compareClusterMeth_BioProc, showCategory=15) +
  geom_text(aes(label=Count, size=0.002, fontface='bold'), colour='white', show.legend = F) +
  scale_size(range = c(3, 9)) +
  theme(legend.key.height = unit(30,'points'),
        axis.text.y = element_text(size=22),
        axis.text.x = element_text(size=18, angle=45, hjust=1),
        text = element_text(size=18)) + 
  ggtitle("GO_biological_process for all clusters")
dev.off() # change in postscript width, height for a good aspect ratio
########################################################################
##########################################################################################
