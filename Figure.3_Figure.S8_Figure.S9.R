##########################################################################################
### Figure 3B - S8 ###
##########################################################################################

# call CpGs from BAM files on Hatchling RRBS
##########################################################################################
library("methylKit")
# create methylRawList.obj
objs=processBismarkAln(location="./Input_Data/Obi_Hatchlings_RRBS_sorted.bam",
                       sample.id="Obi_Hatchlings_RRBS",
                       assembly="PRJNA270931",
                       save.folder="./Output_Data/methylKit_Obimaculoides_280",save.context="CpG",read.context="CpG",
                       nolap=FALSE,mincov=10,minqual=20,phred64=FALSE)

# importing data exported by processBismarkAln as txt file
# containing "freqC" that is essentially the perc.meth
##########################################################################################
obj_perc.meth <- read.delim("./Output_Data/methylKit_Obimaculoides_280/Obi_Hatchlings_RRBS_CpG.txt")

# selecting the column we need for converting in GRanges object
library("dplyr")
library("tibble")
CpGs_All <- obj_perc.meth %>% select(chr, base, strand, freqC) %>%
  add_column(obj_perc.meth$base, .after = "base")

colnames(CpGs_All) <- c("chr", "start", "end", "strand", "perc.meth")

CpGs_All$strand <- gsub("F", "+", CpGs_All$strand)
CpGs_All$strand <- gsub("R", "-", CpGs_All$strand)

View(CpGs_All)

# subsetting by perc.meth
CpGs_All_LowMeth <- subset(CpGs_All, CpGs_All$perc.meth<20)
CpGs_All_HighMeth <- subset(CpGs_All, CpGs_All$perc.meth>80)

# converting to GRanges
GR_CpGs_All <- as(CpGs_All, "GRanges")
save(GR_CpGs_All, file="./R_Objects/GR_CpGs_All.RData")

GR_CpGs_All_LowMeth <- as(CpGs_All_LowMeth, "GRanges")
save(GR_CpGs_All_LowMeth, file="./R_Objects/GR_CpGs_All_LowMeth.RData")

GR_CpGs_All_HighMeth <- as(CpGs_All_HighMeth, "GRanges")
save(GR_CpGs_All_HighMeth, file="./R_Objects/GR_CpGs_All_HighMeth.RData")
##########################################################################################

# import CGmap files from Supra and Sub Esophageal Brain WGBS
##########################################################################################
Octopus_WGBS_Brain_SupraE <- read.table("./Input_Data/GSM4209498_mC_brain_octopus_supraesophageal_onlyCG.CGmap", sep="\t", header = F, quote = "")
Octopus_WGBS_Brain_SubE <- read.table("./Input_Data/GSM4209499_mC_brain_octopus_subesophageal_onlyCG.CGmap", sep="\t", header = F, quote = "")

View(Octopus_WGBS_Brain_SupraE)
View(Octopus_WGBS_Brain_SubE)

# naming column as per standard CGmap file
colnames(Octopus_WGBS_Brain_SupraE) <- c("chr", "strand", "start", "context", "Dinucleotide", "perc.meth", "count_mC", "count_C")
colnames(Octopus_WGBS_Brain_SubE) <- c("chr", "strand", "start", "context", "Dinucleotide", "perc.meth", "count_mC", "count_C")

# select only needed column
library("dplyr")
Octopus_WGBS_Brain_SupraE_select <- Octopus_WGBS_Brain_SupraE %>% select(chr, strand, start, perc.meth) 
View(Octopus_WGBS_Brain_SupraE_select)

Octopus_WGBS_Brain_SubE_select <- Octopus_WGBS_Brain_SubE %>% select(chr, strand, start, perc.meth) 
View(Octopus_WGBS_Brain_SubE_select)

# changing text to make it homogeneous to RRBS nomenclature
Octopus_WGBS_Brain_SupraE_select$strand <- gsub("C", "+", Octopus_WGBS_Brain_SupraE_select$strand)
Octopus_WGBS_Brain_SupraE_select$strand <- gsub("G", "-", Octopus_WGBS_Brain_SupraE_select$strand)

Octopus_WGBS_Brain_SubE_select$strand <- gsub("C", "+", Octopus_WGBS_Brain_SubE_select$strand)
Octopus_WGBS_Brain_SubE_select$strand <- gsub("G", "-", Octopus_WGBS_Brain_SubE_select$strand)

# transforming methylation values in 0 to 100 to make it homogeneous to RRBS
Octopus_WGBS_Brain_SupraE_select$perc.meth <- Octopus_WGBS_Brain_SupraE_select$perc.meth * 100
Octopus_WGBS_Brain_SubE_select$perc.meth <- Octopus_WGBS_Brain_SubE_select$perc.meth * 100

# Only Scaffold chromosomes (remove chrL) and add "end" base pair position column for GRanges conversion
library("dplyr")
library("tibble")
# Supra esophageal Brain
Octopus_WGBS_Brain_SupraE_select <- Octopus_WGBS_Brain_SupraE_select %>% 
  slice(6625:64852310) %>%
  add_column(end= Octopus_WGBS_Brain_SupraE_select$start, .after = "start")
View(Octopus_WGBS_Brain_SupraE_select)

Octopus_WGBS_Brain_SupraE_select_meth <- Octopus_WGBS_Brain_SupraE_select %>% 
  filter(perc.meth >80)
View(Octopus_WGBS_Brain_SupraE_select_meth)
Octopus_WGBS_Brain_SupraE_select_Below80 <- Octopus_WGBS_Brain_SupraE_select %>% 
  filter(perc.meth <80)
View(Octopus_WGBS_Brain_SupraE_select_Below80)

Octopus_WGBS_Brain_SupraE_select_UNmeth <- Octopus_WGBS_Brain_SupraE_select %>% 
  filter(perc.meth <20)
View(Octopus_WGBS_Brain_SupraE_select_UNmeth)

# Sub esophageal Brain
Octopus_WGBS_Brain_SubE_select <- Octopus_WGBS_Brain_SubE_select %>% 
  slice(6625:47060128) %>%
  add_column(end= Octopus_WGBS_Brain_SubE_select$start, .after = "start")
View(Octopus_WGBS_Brain_SubE_select)

Octopus_WGBS_Brain_SubE_select_meth <- Octopus_WGBS_Brain_SubE_select %>% 
  filter(perc.meth >80)
View(Octopus_WGBS_Brain_SubE_select_meth)

Octopus_WGBS_Brain_SubE_select_UNmeth <- Octopus_WGBS_Brain_SubE_select %>% 
  filter(perc.meth <20)
View(Octopus_WGBS_Brain_SubE_select_UNmeth)

# converting to GRanges
GR_Octopus_WGBS_Brain_SupraE_select <- as(Octopus_WGBS_Brain_SupraE_select, "GRanges")
GR_Octopus_WGBS_Brain_SupraE_select_meth <- as(Octopus_WGBS_Brain_SupraE_select_meth, "GRanges")
GR_Octopus_WGBS_Brain_SupraE_select_UNmeth <- as(Octopus_WGBS_Brain_SupraE_select_UNmeth, "GRanges")

GR_Octopus_WGBS_Brain_SubE_select <- as(Octopus_WGBS_Brain_SubE_select, "GRanges")
GR_Octopus_WGBS_Brain_SubE_select_meth <- as(Octopus_WGBS_Brain_SubE_select_meth, "GRanges")
GR_Octopus_WGBS_Brain_SubE_select_UNmeth <- as(Octopus_WGBS_Brain_SubE_select_UNmeth, "GRanges")

save(GR_Octopus_WGBS_Brain_SupraE_select, file="./R_Objects/GR_Octopus_WGBS_Brain_SupraE_select.RData")
save(GR_Octopus_WGBS_Brain_SupraE_select_meth, file="./R_Objects/GR_Octopus_WGBS_Brain_SupraE_select_meth.RData")
save(GR_Octopus_WGBS_Brain_SupraE_select_UNmeth, file="./R_Objects/GR_Octopus_WGBS_Brain_SupraE_select_UNmeth.RData")

save(GR_Octopus_WGBS_Brain_SubE_select, file="./R_Objects/GR_Octopus_WGBS_Brain_SubE_select.RData")
save(GR_Octopus_WGBS_Brain_SubE_select_meth, file="./R_Objects/GR_Octopus_WGBS_Brain_SubE_select_meth.RData")
save(GR_Octopus_WGBS_Brain_SubE_select_UNmeth, file="./R_Objects/GR_Octopus_WGBS_Brain_SubE_select_UNmeth.RData")
##########################################################################################

# All CpGs per Tissues in RRBS and WGBSs (Fig.S8B)
##########################################################################################
sq <- seq(max(length(GR_CpGs_All),
              length(GR_Octopus_WGBS_Brain_SupraE_select),
              length(GR_Octopus_WGBS_Brain_SubE_select)))
# this is forcing to keep the same lenght for all columns and put NAs
df_perc.meth_All_CpGs_Hatchings_SupraE_SubE <- data.frame(GR_CpGs_All$perc.meth[sq],
                                                          GR_Octopus_WGBS_Brain_SupraE_select$perc.meth[sq],
                                                          GR_Octopus_WGBS_Brain_SubE_select$perc.meth[sq])
head(df_perc.meth_All_CpGs_Hatchings_SupraE_SubE)
tail(df_perc.meth_All_CpGs_Hatchings_SupraE_SubE)
colnames(df_perc.meth_All_CpGs_Hatchings_SupraE_SubE) <- c("RRBS_Hatchlings",
                                                           "WGBS_Brain_SupraE",
                                                           "WGBS_Brain_SubE")
View(df_perc.meth_All_CpGs_Hatchings_SupraE_SubE)

dir.create("./Analysis")
write.table(df_perc.meth_All_CpGs_Hatchings_SupraE_SubE, file="./Analysis/perc.meth_All_CpGs_Hatchings_SupraE_SubE.txt", quote = FALSE, sep="\t", row.names = FALSE)

library("reshape2")
df_perc.meth_All_CpGs_Hatchings_SupraE_SubE_melt <- melt(df_perc.meth_All_CpGs_Hatchings_SupraE_SubE)
head(df_perc.meth_All_CpGs_Hatchings_SupraE_SubE_melt)
tail(df_perc.meth_All_CpGs_Hatchings_SupraE_SubE_melt)

# violin plot 
###############################################################
library("ggplot2")
dir.create("./Results/ggplots")
dir.create("./Results/ggplots/CpGs_All")

png("./Results/ggplots/CpGs_All/CpGs_All_CpGs_Hatchings_SupraE_SubE_ViolinPlot_Width.png", width= 10, height= 7.5, units="in", res=300)
ggplot(df_perc.meth_All_CpGs_Hatchings_SupraE_SubE_melt, aes(x = variable, y = value, fill=variable)) + 
  geom_violin(trim = TRUE, scale = "width") + scale_y_continuous(breaks=seq(0,100,25)) + 
  labs(x= "Density of CpGs per each sample", y = "% of methylation") + 
  theme_classic() + theme(axis.text.x = element_text(angle=30, hjust = 1)) + 
  scale_fill_manual(values = c("#009091","#7DBC75","#336633")) + 
  geom_boxplot(width=0.05, fill="white", notch = FALSE)
dev.off()
##########################################################################################

# Common CpGs to all three tissues across RRBS and WGBSs (Fig.3B and Fig.S8C-D)
##########################################################################################
library("GenomicRanges")
# merge by overlap is creating a dataset with both the input maintained in the output
RRBS_WGBS_Brain_SupraE_merged <- mergeByOverlaps(GR_CpGs_All, GR_Octopus_WGBS_Brain_SupraE_select)
RRBS_WGBS_Brain_SubE_merged <- mergeByOverlaps(GR_CpGs_All, GR_Octopus_WGBS_Brain_SubE_select)

# extract the column needed and recreate a GR for the second merge
###############################################################
df_RRBS_WGBS_Brain_SupraE_merged <- data.frame(seqnames = RRBS_WGBS_Brain_SupraE_merged@listData$GR_CpGs_All@seqnames, 
                                               ranges = RRBS_WGBS_Brain_SupraE_merged@listData$GR_CpGs_All@ranges, 
                                               strand = RRBS_WGBS_Brain_SupraE_merged@listData$GR_CpGs_All@strand, 
                                               perc.meth_Hatchlings = RRBS_WGBS_Brain_SupraE_merged@listData$GR_CpGs_All$perc.meth, 
                                               perc.meth_SupraE = RRBS_WGBS_Brain_SupraE_merged@listData$GR_Octopus_WGBS_Brain_SupraE_select$perc.meth)
GR_df_RRBS_WGBS_Brain_SupraE_merged <- as(df_RRBS_WGBS_Brain_SupraE_merged, "GRanges")

df_RRBS_WGBS_Brain_SubE_merged <- data.frame(seqnames = RRBS_WGBS_Brain_SubE_merged@listData$GR_CpGs_All@seqnames, 
                                             ranges = RRBS_WGBS_Brain_SubE_merged@listData$GR_CpGs_All@ranges, 
                                             strand = RRBS_WGBS_Brain_SubE_merged@listData$GR_CpGs_All@strand, 
                                             perc.meth_Hatchlings = RRBS_WGBS_Brain_SubE_merged@listData$GR_CpGs_All$perc.meth, 
                                             perc.meth_SubE = RRBS_WGBS_Brain_SubE_merged@listData$GR_Octopus_WGBS_Brain_SubE_select$perc.meth)
GR_df_RRBS_WGBS_Brain_SubE_merged <- as(df_RRBS_WGBS_Brain_SubE_merged, "GRanges")

# second merge between GR and extraction of column needed
###############################################################
# here the CpGs for hatcling will be reduntant so take only one
RRBS_WGBSs_commonCpGs_merged <- mergeByOverlaps(GR_df_RRBS_WGBS_Brain_SupraE_merged, GR_df_RRBS_WGBS_Brain_SubE_merged)

# extract CpGs perc.meth only for ggplot
df_RRBS_WGBSs_commonCpGs_merged <- data.frame(perc.meth_Hatchlings = RRBS_WGBSs_commonCpGs_merged@listData$GR_df_RRBS_WGBS_Brain_SupraE_merged$perc.meth_Hatchlings, 
                                              perc.meth_SupraE = RRBS_WGBSs_commonCpGs_merged@listData$GR_df_RRBS_WGBS_Brain_SupraE_merged$perc.meth_SupraE,
                                              perc.meth_SubE = RRBS_WGBSs_commonCpGs_merged@listData$GR_df_RRBS_WGBS_Brain_SubE_merged$perc.meth_SubE)
###############################################################

# scatter plot 
###############################################################
library("ggplot2")
dir.create("./Results/ggplots")
dir.create("./Results/ggplots/CpGs_common")

png("./Results/ggplots/CpGs_common/commonCpGs_SupraE.vs.SubE_perc.meth_Color.Hatchling_ScatterPlot_CountSupraE6.png", width= 10, height= 7.5, units="in", res=300)
ggplot(df_RRBS_WGBSs_commonCpGs_merged, aes(x = perc.meth_SupraE, y = perc.meth_SubE, color= perc.meth_Hatchlings)) + 
  geom_count(aes(size = after_stat(prop), group = perc.meth_SupraE)) +
  scale_size_area(max_size = 2.5) + geom_smooth(color="black", lwd=2) + 
  labs(x= "perc.meth_SupraE", y = "perc.meth_SubE") + 
  theme_classic() + scale_color_gradient(low="#0C7BDC", high="#FF950A")
dev.off()

png("./Results/ggplots/CpGs_common/commonCpGs_SupraE.vs.Hatchling_perc.meth_Color.SubE_ScatterPlot_CountSupraE6.png", width= 10, height= 7.5, units="in", res=300)
ggplot(df_RRBS_WGBSs_commonCpGs_merged, aes(x = perc.meth_SupraE, y = perc.meth_Hatchlings, color= perc.meth_SubE)) + 
  geom_count(aes(size = after_stat(prop), group = perc.meth_SupraE)) +
  scale_size_area(max_size = 2.5) + geom_smooth(color="black", lwd=2) + 
  labs(x= "perc.meth_SupraE", y = "perc.meth_Hatchlings") + 
  theme_classic() + scale_color_gradient(low="#0C7BDC", high="#FF950A")
dev.off()

png("./Results/ggplots/CpGs_common/commonCpGs_SubE.vs.Hatchling_perc.meth_Color.SupraE_ScatterPlot_CountSubE6.png", width= 10, height= 7.5, units="in", res=300)
ggplot(df_RRBS_WGBSs_commonCpGs_merged, aes(x = perc.meth_SubE, y = perc.meth_Hatchlings, color= perc.meth_SupraE)) + 
  geom_count(aes(size = after_stat(prop), group = perc.meth_SubE)) +
  scale_size_area(max_size = 2.5) + geom_smooth(color="black", lwd=2) + 
  labs(x= "perc.meth_SubE", y = "perc.meth_Hatchlings") + 
  theme_classic() + scale_color_gradient(low="#0C7BDC", high="#FF950A")
dev.off()
##########################################################################################



##########################################################################################
### Figure 3C - S9 ###
##########################################################################################

# generate annotation for full-lenght transcripts, gene.parts (promoter, exon, intron) and intergenic
##########################################################################################
library("GenomicFeatures")
gff3.file = "./Input_Data/Obimaculoides_280.gene.gff3"
txdb <- makeTxDbFromGFF(file=gff3.file, organism="Octopus bimaculoides")
GR_transcripts_OB <- transcripts(txdb)
save(GR_transcripts_OB, file="./R_Objects/GR_transcripts_OB.RData")

# export the txdb as a bed file to be re-imported and analyzed by genomation directly
###################################################
library("rtracklayer")
txdb_bed_out <- file.path("./Output_Data/bed", "txdb_out.bed")
export.bed(object=txdb, con=txdb_bed_out)

library("genomation")
bed.file = "./Output_Data/bed/txdb_out.bed"
gene.parts = readTranscriptFeatures(bed.file)

GR_TSS <- gene.parts$TSSes
GR_promoters <- gene.parts$promoters
GR_exons <- gene.parts$exons
GR_introns <- gene.parts$introns
save(GR_promoters, file="./R_Objects/GR_promoters.RData")
save(GR_exons, file="./R_Objects/GR_exons.RData")
save(GR_introns, file="./R_Objects/GR_introns.RData")

### get intergenic regions from full scaffold length from OB genome
###################################################
library("BSgenome.Obimaculoides.OIST.Obimaculoides.280") # this need to be build manually
GR_Scaffolds <- as(seqinfo(Obimaculoides), "GRanges")
GR_transcripts_OB_collapsed <- reduce(transcripts(txdb))
save(GR_Scaffolds, file="./R_Objects/GR_Scaffolds.RData")
save(GR_transcripts_OB_collapsed, file="./R_Objects/GR_transcripts_OB_collapsed.RData")

## modified strand to * or use unstrand(GR_transcripts_OB_collapsed)
strand(GR_transcripts_OB_collapsed) <- "*"
GR_intergenic <- GenomicRanges::setdiff(GR_Scaffolds, GR_transcripts_OB_collapsed)
save(GR_intergenic, file="./R_Objects/GR_intergenic.RData")
##########################################################################################

# generate annotation for whole genome CpGs
##########################################################################################
library("Biostrings")
library("GenomicRanges")
library("IRanges")

fastaFile <- readDNAStringSet("./Input_Data/PRJNA270931.fa")
seqnames = names(fastaFile)

CG_string <- DNAString("CG")
cgs <- vmatchPattern(CG_string, fastaFile)

cpgr <- as(cgs, "GRanges")
save(cpgr, file="./Output_Data/GR_CpGs_Obimaculoides_Obimaculoides_280.Rdata")
##########################################################################################

# Annotation by genomation package
##########################################################################################
library("genomation")
dir.create("./Results/Genomation")
dir.create("./Results/Genomation/AnnotationPlot")

### whole genome CpGs
###################################################
annot_cpgr_noStrand = annotateWithGeneParts(cpgr, gene.parts, strand=FALSE, intersect.chr=TRUE)
postscript("./Results/Genomation/AnnotationPlot/cpgr_noStrand_AnnotationPlot.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_cpgr_noStrand)
dev.off()
# take the numbers (with promoter > exon > intron precedence)
annot_cpgr_noStrand
###################################################

### Hatchling RRBS CpGs
###################################################
### All CpGs
annot_GR_CpGs_All_noStrand = annotateWithGeneParts(GR_CpGs_All, gene.parts, strand=FALSE, intersect.chr=TRUE)
postscript("./Results/Genomation/AnnotationPlot/GR_CpGs_All_noStrand_RRBS_AnnotationPlot.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_CpGs_All_noStrand)
dev.off()
# take the numbers (with promoter > exon > intron precedence)
annot_GR_CpGs_All_noStrand

### Low Meth CpGs
annot_GR_CpGs_All_LowMeth_noStrand = annotateWithGeneParts(GR_CpGs_All_LowMeth, gene.parts, strand=FALSE, intersect.chr=TRUE)
postscript("./Results/Genomation/AnnotationPlot/GR_CpGs_All_LowMeth_noStrand_RRBS_AnnotationPlot.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_CpGs_All_LowMeth_noStrand)
dev.off()
# take the numbers (with promoter > exon > intron precedence)
annot_GR_CpGs_All_LowMeth_noStrand

### High Meth CpGs
annot_GR_CpGs_All_HighMeth_noStrand = annotateWithGeneParts(GR_CpGs_All_HighMeth, gene.parts, strand=FALSE, intersect.chr=TRUE)
postscript("./Results/Genomation/AnnotationPlot/GR_CpGs_All_HighMeth_noStrand_RRBS_AnnotationPlot.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_CpGs_All_HighMeth_noStrand)
dev.off()
# take the numbers (with promoter > exon > intron precedence)
annot_GR_CpGs_All_HighMeth_noStrand
###################################################

### SupraE brain WGBS CpGs
###################################################
### All CpGs
annot_GR_WGBS_Brain_SupraE_noStrand = annotateWithGeneParts(GR_Octopus_WGBS_Brain_SupraE_select, gene.parts, strand=FALSE, intersect.chr=TRUE)
postscript("./Results/Genomation/AnnotationPlot/GR_WGBS_Brain_SupraE_noStrand_AnnotationPlot.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_WGBS_Brain_SupraE_noStrand)
dev.off()
# take the numbers (with promoter > exon > intron precedence)
annot_GR_WGBS_Brain_SupraE_noStrand
save(annot_GR_WGBS_Brain_SupraE_noStrand, file="./R_Objects/gene.parts_annot_GR_WGBS_Brain_SupraE_noStrand.RData")
rm(annot_GR_WGBS_Brain_SupraE_noStrand) # removing is suggested as this object is data heavy

### Low Meth CpGs
annot_GR_WGBS_Brain_SupraE_LowMeth_noStrand = annotateWithGeneParts(GR_Octopus_WGBS_Brain_SupraE_select_UNmeth, gene.parts, strand=FALSE, intersect.chr=TRUE)
postscript("./Results/Genomation/AnnotationPlot/GR_WGBS_Brain_SupraE_LowMeth_noStrand_AnnotationPlot.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_WGBS_Brain_SupraE_LowMeth_noStrand)
dev.off()
# take the numbers (with promoter > exon > intron precedence)
annot_GR_WGBS_Brain_SupraE_LowMeth_noStrand
save(annot_GR_WGBS_Brain_SupraE_LowMeth_noStrand, file="./R_Objects/gene.parts_annot_GR_WGBS_Brain_SupraE_LowMeth_noStrand.RData")
rm(annot_GR_WGBS_Brain_SupraE_LowMeth_noStrand) # removing is suggested as this object is data heavy

### High Meth CpGs
annot_GR_WGBS_Brain_SupraE_HighMeth_noStrand = annotateWithGeneParts(GR_Octopus_WGBS_Brain_SupraE_select_meth, gene.parts, strand=FALSE, intersect.chr=TRUE)
postscript("./Results/Genomation/AnnotationPlot/GR_WGBS_Brain_SupraE_HighMeth_noStrand_AnnotationPlot.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_WGBS_Brain_SupraE_HighMeth_noStrand)
dev.off()
# take the numbers (with promoter > exon > intron precedence)
annot_GR_WGBS_Brain_SupraE_HighMeth_noStrand
save(annot_GR_WGBS_Brain_SupraE_HighMeth_noStrand, file="./R_Objects/gene.parts_annot_GR_WGBS_Brain_SupraE_HighMeth_noStrand.RData")
rm(annot_GR_WGBS_Brain_SupraE_HighMeth_noStrand) # removing is suggested as this object is data heavy
###################################################

### SubE brain WGBS CpGs
###################################################
### All CpGs
annot_GR_WGBS_Brain_SubE_noStrand = annotateWithGeneParts(GR_Octopus_WGBS_Brain_SubE_select, gene.parts, strand=FALSE, intersect.chr=TRUE)
postscript("./Results/Genomation/AnnotationPlot/GR_WGBS_Brain_SubE_noStrand_AnnotationPlot.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_WGBS_Brain_SubE_noStrand)
dev.off()
# take the numbers (with promoter > exon > intron precedence)
annot_GR_WGBS_Brain_SubE_noStrand
save(annot_GR_WGBS_Brain_SubE_noStrand, file="./R_Objects/gene.parts_annot_GR_WGBS_Brain_SubE_noStrand.RData")
rm(annot_GR_WGBS_Brain_SubE_noStrand) # removing is suggested as this object is data heavy

### Low Meth CpGs
annot_GR_WGBS_Brain_SubE_LowMeth_noStrand = annotateWithGeneParts(GR_Octopus_WGBS_Brain_SubE_select_UNmeth, gene.parts, strand=FALSE, intersect.chr=TRUE)
postscript("./Results/Genomation/AnnotationPlot/GR_WGBS_Brain_SubE_LowMeth_noStrand_AnnotationPlot.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_WGBS_Brain_SubE_LowMeth_noStrand)
dev.off()
# take the numbers (with promoter > exon > intron precedence)
annot_GR_WGBS_Brain_SubE_LowMeth_noStrand
save(annot_GR_WGBS_Brain_SubE_LowMeth_noStrand, file="./R_Objects/gene.parts_annot_GR_WGBS_Brain_SubE_LowMeth_noStrand.RData")
rm(annot_GR_WGBS_Brain_SubE_LowMeth_noStrand) # removing is suggested as this object is data heavy

### High Meth CpGs
annot_GR_WGBS_Brain_SubE_HighMeth_noStrand = annotateWithGeneParts(GR_Octopus_WGBS_Brain_SubE_select_meth, gene.parts, strand=FALSE, intersect.chr=TRUE)
postscript("./Results/Genomation/AnnotationPlot/GR_WGBS_Brain_SubE_HighMeth_noStrand_AnnotationPlot.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotTargetAnnotation(annot_GR_WGBS_Brain_SubE_HighMeth_noStrand)
dev.off()
# take the numbers (with promoter > exon > intron precedence)
annot_GR_WGBS_Brain_SubE_HighMeth_noStrand
save(annot_GR_WGBS_Brain_SubE_HighMeth_noStrand, file="./R_Objects/gene.parts_annot_GR_WGBS_Brain_SubE_HighMeth_noStrand.RData")
rm(annot_GR_WGBS_Brain_SubE_HighMeth_noStrand) # removing is suggested as this object is data heavy
###################################################
##########################################################################################



##########################################################################################
### Figure 3D-E - S9 ###
##########################################################################################

# import Repeat Elements from gff file
##########################################################################################
library("rtracklayer")
my_file <- "./Input_Data/oct.m51.D6.c51.b28.5.p3.p3.p3.p10.fa.gff"
# read as data.frame
RE_df <- readGFF(my_file)

# reordering and renaming columns for export in BED
###################################################
library("dplyr")
RE_df <- RE_df %>% 
  select(seqid, start, end, 
         strand, score, Target) 

colnames(RE_df) <- c("chr", "start", "end", "strand", "score", "name")
View(RE_df)

RE_GR <- as (RE_df, "GRanges")
save(RE_GR, file="./R_Objects/RE_GR.RData")

library("rtracklayer")
RE_GR_bed_out <- file.path("./Output_Data/bed", "RE_GR_out.bed")
export.bed(object=RE_GR, con=RE_GR_bed_out)
###################################################

# call Transposon and Not_Trasposon 
###################################################
library("dplyr")
library("stringr")
RE_Transposons <- RE_df %>% 
  filter(str_detect(name, regex(c("SINE|LINE|RTE|Tx1|CR1|_R4_|L1|LTR|ERV|DNA|EnSpm|hAT|Kolobok|Helitron"), ignore_case = TRUE)))
View(RE_Transposons)

RE_Transposons_GR <- as(RE_Transposons, "GRanges")
save(RE_Transposons_GR, file="./R_Objects/RE_Transposons_GR.RData")

RE_Not_Transposons <- RE_df %>% 
  filter(str_detect(name, regex(c("Satellite|Simple_repeat|Low_complexity|Unknown|rRNA|tRNA-|snRNA"), ignore_case = TRUE)))
View(RE_Not_Transposons)

RE_Not_Transposons_GR <- as(RE_Not_Transposons, "GRanges")
save(RE_Not_Transposons_GR, file="./R_Objects/RE_Not_Transposons_GR.RData")

AllRepeats_GR <- merge(RE_Transposons_GR,
                       RE_Not_Transposons_GR,
                       all=TRUE)
save(AllRepeats_GR, file="./R_Objects/AllRepeats_GR.RData")
###################################################

# call Transposon and Not_Trasposon divided by family
###################################################
# using betools software on HPC to intersect CpGs and REs by chromosome position 
# this job is to heavy to be run in local computer or in R

################################
##### run on linux teminal #####
################################
# reorder the chr names by 1 2 3 instead of 1 10 11 with 2 col ordered by number position etc.

#sort -k1,1V -k2,2n RE_GR_out.bed > RE_GR_out.versionsorted.bed

# intersect reporting both the dataset and -s for strand.aware
# -sorted can give an error if the dataset do not start or contain the same chr

#bedtools intersect -wa -wb \
#-a GR_CpGs_All_score.bed \
#-b RE_GR_out.versionsorted.bed \
#-s -sorted > GR_CpGs_All_RE_GR_bothAandB.bed

################################
##### stop run on linux teminal #####
################################

# importing the intersected file as data.frame for Hatchling RRBS
###################################################
CpGs_All_RE_GR_bothAandB <- read.delim("./Output_Data/bed/GR_CpGs_All_RE_GR_bothAandB.bed", header = FALSE)
colnames(CpGs_All_RE_GR_bothAandB) <- c("chr", "start", "end", "name", "score", "strand",
                                        "chr_RE", "start_RE", "end_RE", "name_RE", "score_RE", "strand_RE")
View(CpGs_All_RE_GR_bothAandB)

# Filter by a string pattern specific REs family for Hatchling RRBS
###################################################
library("dplyr")
library("stringr")
# SINE
############################################################
CpGs_All_RE_GR_bothAandB_SINE <- CpGs_All_RE_GR_bothAandB %>% 
  filter(str_detect(name_RE, regex(c("SINE"), ignore_case = TRUE)))
View(CpGs_All_RE_GR_bothAandB_SINE)

CpGs_All_RE_SINE <- CpGs_All_RE_GR_bothAandB_SINE %>% 
  select(chr, start, end, strand, score) 
View(CpGs_All_RE_SINE)

CpGs_All_RE_SINE_GR <- as(CpGs_All_RE_SINE, "GRanges")
############################################################

# LINE
############################################################
CpGs_All_RE_GR_bothAandB_LINE <- CpGs_All_RE_GR_bothAandB %>% 
  filter(str_detect(name_RE, regex(c("LINE|RTE|Tx1|CR1|_R4_|L1"), ignore_case = TRUE)))
View(CpGs_All_RE_GR_bothAandB_LINE)

CpGs_All_RE_LINE <- CpGs_All_RE_GR_bothAandB_LINE %>% 
  select(chr, start, end, strand, score) 
View(CpGs_All_RE_LINE)

CpGs_All_RE_LINE_GR <- as(CpGs_All_RE_LINE, "GRanges")
############################################################

# LTR
############################################################
CpGs_All_RE_GR_bothAandB_LTR <- CpGs_All_RE_GR_bothAandB %>% 
  filter(str_detect(name_RE, regex(c("LTR|ERV"), ignore_case = TRUE)))
View(CpGs_All_RE_GR_bothAandB_LTR)

CpGs_All_RE_LTR <- CpGs_All_RE_GR_bothAandB_LTR %>% 
  select(chr, start, end, strand, score) 
View(CpGs_All_RE_LTR)

CpGs_All_RE_LTR_GR <- as(CpGs_All_RE_LTR, "GRanges")
############################################################

# DNA
############################################################
CpGs_All_RE_GR_bothAandB_DNA <- CpGs_All_RE_GR_bothAandB %>% 
  filter(str_detect(name_RE, regex(c("DNA|EnSpm|hAT|Kolobok|Helitron"), ignore_case = TRUE)))
View(CpGs_All_RE_GR_bothAandB_DNA)

CpGs_All_RE_DNA <- CpGs_All_RE_GR_bothAandB_DNA %>% 
  select(chr, start, end, strand, score) 
View(CpGs_All_RE_DNA)

CpGs_All_RE_DNA_GR <- as(CpGs_All_RE_DNA, "GRanges")
############################################################

# Satellite
############################################################
CpGs_All_RE_GR_bothAandB_Satellite <- CpGs_All_RE_GR_bothAandB %>% 
  filter(str_detect(name_RE, regex(c("Satellite"), ignore_case = TRUE)))
View(CpGs_All_RE_GR_bothAandB_Satellite)

CpGs_All_RE_Satellite <- CpGs_All_RE_GR_bothAandB_Satellite %>% 
  select(chr, start, end, strand, score) 
View(CpGs_All_RE_Satellite)

CpGs_All_RE_Satellite_GR <- as(CpGs_All_RE_Satellite, "GRanges")
############################################################

# MicroSatellite = Simple Repeat
############################################################
CpGs_All_RE_GR_bothAandB_MicroSatellite <- CpGs_All_RE_GR_bothAandB %>% 
  filter(str_detect(name_RE, regex(c("Simple_repeat"), ignore_case = TRUE)))
View(CpGs_All_RE_GR_bothAandB_MicroSatellite)

CpGs_All_RE_MicroSatellite <- CpGs_All_RE_GR_bothAandB_MicroSatellite %>% 
  select(chr, start, end, strand, score) 
View(CpGs_All_RE_MicroSatellite)

CpGs_All_RE_MicroSatellite_GR <- as(CpGs_All_RE_MicroSatellite, "GRanges")
############################################################

# Others
############################################################
CpGs_All_RE_GR_bothAandB_Others <- CpGs_All_RE_GR_bothAandB %>% 
  filter(str_detect(name_RE, regex(c("Low_complexity|Unknown|rRNA|tRNA-|snRNA"), ignore_case = TRUE)))
View(CpGs_All_RE_GR_bothAandB_Others)

CpGs_All_RE_Others <- CpGs_All_RE_GR_bothAandB_Others %>% 
  select(chr, start, end, strand, score) 
View(CpGs_All_RE_Others)

CpGs_All_RE_Others_GR <- as(CpGs_All_RE_Others, "GRanges")
############################################################

# merge CpGs diveded in each family by Transposon ans not-Transposon for Hatchling RRBS
###################################################
library("GenomicRanges")
CpGs_All_Trasposon_GR <- merge(CpGs_All_RE_DNA_GR,
                               CpGs_All_RE_LTR_GR,
                               CpGs_All_RE_LINE_GR,
                               CpGs_All_RE_SINE_GR,
                               all=TRUE)

CpGs_All_Not_Trasposon_GR <-  merge(CpGs_All_RE_Satellite_GR,
                                    CpGs_All_RE_MicroSatellite_GR,
                                    CpGs_All_RE_Others_GR,
                                    all=TRUE)

CpGs_All_AllRepeats_GR <- merge(CpGs_All_Trasposon_GR,
                                CpGs_All_Not_Trasposon_GR,
                                all=TRUE)
###################################################


# calculating matrix and plot the line-plot for all Repeat Elements (Fig.3D)
##########################################################################################
# Transposons All for Hatchling RRBS
########################################################################
# CpGs_All_Trasposon_GR has score column instead of perc.meth
ScoreBin_CpG_All_RE_GR_Transposons_All = ScoreMatrixBin(target=CpGs_All_Trasposon_GR,
                                                        windows=RE_Transposons_GR,
                                                        bin.num=15, strand.aware=TRUE,
                                                        weight.col="score",
                                                        is.noCovNA = TRUE)

pdf("./Results/Genomation/MetaPlot/RE/ScoreBin_CpG_All_RE_GR_Transposons_All.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_Transposons_All,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_Transposons_All), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################

# Transposons All for SupraE WGBS
########################################################################
ScoreBin_CpG_All_RE_GR_Transposons_All = ScoreMatrixBin(target=GR_Octopus_WGBS_Brain_SupraE_select,
                                                        windows=RE_Transposons_GR,
                                                        bin.num=15, strand.aware=TRUE,
                                                        weight.col="perc.meth",
                                                        is.noCovNA = TRUE)

pdf("./Results/Genomation/MetaPlot/RE/ScoreBin_CpG_All_RE_GR_Transposons_All.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_Transposons_All,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_Transposons_All), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################

# Transposons All for SubE WGBS
########################################################################
ScoreBin_CpG_All_RE_GR_Transposons_All = ScoreMatrixBin(target=GR_Octopus_WGBS_Brain_SubE_select,
                                                        windows=RE_Transposons_GR,
                                                        bin.num=15, strand.aware=TRUE,
                                                        weight.col="perc.meth",
                                                        is.noCovNA = TRUE)

pdf("./Results/Genomation/MetaPlot/RE/ScoreBin_CpG_All_RE_GR_Transposons_All.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_Transposons_All,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_Transposons_All), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################

# CpG density for whole genome CpGs
########################################################################
ScoreBin_WG_CpG_RE_GR_Transposons_All = ScoreMatrixBin(target=cpgr,
                                                       windows=RE_Transposons_GR,
                                                       bin.num=15, strand.aware=TRUE
)

pdf("./Results/Genomation/MetaPlot/RE/ScoreBin_WG_CpG_RE_GR_Transposons_All.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_WG_CpG_RE_GR_Transposons_All,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 0.06), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_WG_CpG_RE_GR_Transposons_All), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "CpG density (Avg. score)", xlab = "binned base")
dev.off()
########################################################################
##########################################################################################


# calculating matrix and plot the line-plot for all Transcripts  (Fig.3E)
##########################################################################################
### Hatchling RRBS CpGs
###################################################
library("genomation")
ScoreBin_CpG_All_GR_transcripts_OB = ScoreMatrixBin(target=GR_CpGs_All,
                                                    windows=GR_transcripts_OB,
                                                    bin.num=30, strand.aware=TRUE,
                                                    weight.col="perc.meth",
                                                    is.noCovNA = TRUE)
dir.create("./Results/Genomation/MetaPlot")
pdf("./Results/Genomation/MetaPlot/ScoreBin_CpG_All_GR_transcripts_OB.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_GR_transcripts_OB,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_GR_transcripts_OB), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
###################################################

### SupraE brain WGBS CpGs
###################################################
ScoreBin_WGBS_Brain_SupraE_GR_transcripts_OB = ScoreMatrixBin(target=GR_Octopus_WGBS_Brain_SupraE_select,
                                                              windows=GR_transcripts_OB,
                                                              bin.num=30, strand.aware=TRUE,
                                                              weight.col="perc.meth",
                                                              is.noCovNA = TRUE)

pdf("./Results/Genomation/MetaPlot/ScoreBin_WGBS_Brain_SupraE_GR_transcripts_OB.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_WGBS_Brain_SupraE_GR_transcripts_OB,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_WGBS_Brain_SupraE_GR_transcripts_OB), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
###################################################

### SubE brain WGBS CpGs
###################################################
ScoreBin_WGBS_Brain_SubE_GR_transcripts_OB = ScoreMatrixBin(target=GR_Octopus_WGBS_Brain_SubE_select,
                                                            windows=GR_transcripts_OB,
                                                            bin.num=30, strand.aware=TRUE,
                                                            weight.col="perc.meth",
                                                            is.noCovNA = TRUE)

pdf("./Results/Genomation/MetaPlot/ScoreBin_WGBS_Brain_SubE_GR_transcripts_OB.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_WGBS_Brain_SubE_GR_transcripts_OB,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_WGBS_Brain_SubE_GR_transcripts_OB), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
###################################################
##########################################################################################

# boxplot for All RE divided by family (Fig.S9A-B) for Hatchling RRBS
##########################################################################################
# All CpGs in each family
############################################################
sq <- seq(max(length(CpGs_All_RE_DNA_GR),
              length(CpGs_All_RE_LTR_GR),
              length(CpGs_All_RE_LINE_GR),
              length(CpGs_All_RE_SINE_GR),
              length(CpGs_All_RE_Satellite_GR),
              length(CpGs_All_RE_MicroSatellite_GR),
              length(CpGs_All_RE_Others_GR)))
# this is forcing to keep the same lenght for all columns and put NAs
df_perc.meth_All_RE_groups_Hatchlings <- data.frame(CpGs_All_RE_DNA_GR$score[sq],
                                                    CpGs_All_RE_LTR_GR$score[sq],
                                                    CpGs_All_RE_LINE_GR$score[sq],
                                                    CpGs_All_RE_SINE_GR$score[sq],
                                                    CpGs_All_RE_Satellite_GR$score[sq],
                                                    CpGs_All_RE_MicroSatellite_GR$score[sq],
                                                    CpGs_All_RE_Others_GR$score[sq])
head(df_perc.meth_All_RE_groups_Hatchlings)
tail(df_perc.meth_All_RE_groups_Hatchlings)
colnames(df_perc.meth_All_RE_groups_Hatchlings) <- c("RE-DNA",
                                                     "RE-LTR",
                                                     "RE-LINE",
                                                     "RE-SINE",
                                                     "RE-Satellite",
                                                     "RE-MicroSatellite",
                                                     "RE-Others")
View(df_perc.meth_All_RE_groups_Hatchlings)

dir.create("./Analysis")
write.table(df_perc.meth_All_RE_groups_Hatchlings, file="./Analysis/perc.meth_All_RE_groups_Hatchlings.txt", quote = FALSE, sep="\t", row.names = FALSE)

library("reshape2")
df_perc.meth_All_RE_groups_Hatchlings_melt <- melt(df_perc.meth_All_RE_groups_Hatchlings)
head(df_perc.meth_All_RE_groups_Hatchlings_melt)
tail(df_perc.meth_All_RE_groups_Hatchlings_melt)

library("ggplot2")
dir.create("./Results/ggplots")
dir.create("./Results/ggplots/CpGs_All_RE")

png("./Results/ggplots/CpGs_All_RE/CpGs_All_RE_groups_Hatchlings_BoxPlot.png", width= 10, height= 7.5, units="in", res=300)
ggplot(df_perc.meth_All_RE_groups_Hatchlings_melt, aes(x = variable, y = value, fill=variable)) + 
  geom_boxplot(width=0.5, notch = FALSE) + scale_y_continuous(breaks=seq(0,100,25)) + 
  labs(x= "All RE families", y = "% of methylation") + 
  theme_classic() + theme(axis.text.x = element_text(angle=30, hjust = 1)) + 
  scale_fill_manual(values = c("#1E141E", "#B980AA", "#8A4F8F", "#522A5C", "#D2D2D7", "#8F908F", "#737373"))
dev.off()
############################################################

# Random Sampling to have CpGs in each family of same length of the lowest populated
############################################################
length(CpGs_All_RE_DNA_GR)
# 9065
length(CpGs_All_RE_LTR_GR)
# 4751
length(CpGs_All_RE_LINE_GR)
# 6401
length(CpGs_All_RE_SINE_GR)
# 586521
length(CpGs_All_RE_Satellite_GR)
# 385
length(CpGs_All_RE_MicroSatellite_GR)
# 79157
length(CpGs_All_RE_Others_GR)
# 178258

# CpGs_All_RE_DNA_GR are 9065 so create the sampling or pit the number directly below
CpGs_All_RE_DNA_GR_sampling <- 1:9065

# since CpGs_All_RE_Satellite_GR are 385 and are the limiting number
# we downsample all the other to this number
sampling <- sample(CpGs_All_RE_DNA_GR_sampling, size=385)

# or directly put the number
sampling <- sample(1:9065, size=385)
CpGs_All_RE_DNA_GR_RandomSample <- CpGs_All_RE_DNA_GR[sampling, ] 

sampling <- sample(1:4751, size=385)
CpGs_All_RE_LTR_GR_RandomSample <- CpGs_All_RE_LTR_GR[sampling, ] 

sampling <- sample(1:6401, size=385)
CpGs_All_RE_LINE_GR_RandomSample <- CpGs_All_RE_LINE_GR[sampling, ] 

sampling <- sample(1:586521, size=385)
CpGs_All_RE_SINE_GR_RandomSample <- CpGs_All_RE_SINE_GR[sampling, ] 

sampling <- sample(1:79157, size=385)
CpGs_All_RE_MicroSatellite_GR_RandomSample <- CpGs_All_RE_MicroSatellite_GR[sampling, ] 

sampling <- sample(1:178258, size=385)
CpGs_All_RE_Others_GR_RandomSample <- CpGs_All_RE_Others_GR[sampling, ] 


df_perc.meth_All_RE_groups_Hatchlings_RandomSample <- data.frame(CpGs_All_RE_DNA_GR_RandomSample$score,
                                                                 CpGs_All_RE_LTR_GR_RandomSample$score,
                                                                 CpGs_All_RE_LINE_GR_RandomSample$score,
                                                                 CpGs_All_RE_SINE_GR_RandomSample$score,
                                                                 CpGs_All_RE_Satellite_GR$score,
                                                                 CpGs_All_RE_MicroSatellite_GR_RandomSample$score,
                                                                 CpGs_All_RE_Others_GR_RandomSample$score)
head(df_perc.meth_All_RE_groups_Hatchlings_RandomSample)
tail(df_perc.meth_All_RE_groups_Hatchlings_RandomSample)
colnames(df_perc.meth_All_RE_groups_Hatchlings_RandomSample) <- c("RE-DNA",
                                                                  "RE-LTR",
                                                                  "RE-LINE",
                                                                  "RE-SINE",
                                                                  "RE-Satellite",
                                                                  "RE-MicroSatellite",
                                                                  "RE-Others")
View(df_perc.meth_All_RE_groups_Hatchlings_RandomSample)

dir.create("./Analysis")
write.table(df_perc.meth_All_RE_groups_Hatchlings_RandomSample, file="./Analysis/perc.meth_All_RE_groups_Hatchlings_RandomSample.txt", quote = FALSE, sep="\t", row.names = FALSE)

library("reshape2")
df_perc.meth_All_RE_groups_Hatchlings_RandomSample_melt <- melt(df_perc.meth_All_RE_groups_Hatchlings_RandomSample)
head(df_perc.meth_All_RE_groups_Hatchlings_RandomSample_melt)
tail(df_perc.meth_All_RE_groups_Hatchlings_RandomSample_melt)

library("ggplot2")
dir.create("./Results/ggplots")
dir.create("./Results/ggplots/CpGs_All_RE")

png("./Results/ggplots/CpGs_All_RE/CpGs_All_RE_groups_Hatchlings_RandomSample_BoxPlot.png", width= 10, height= 7.5, units="in", res=300)
ggplot(df_perc.meth_All_RE_groups_Hatchlings_RandomSample_melt, aes(x = variable, y = value, fill=variable)) + 
  geom_boxplot(width=0.5, notch = FALSE) + scale_y_continuous(breaks=seq(0,100,25)) + 
  labs(x= "All RE families", y = "% of methylation") + 
  theme_classic() + theme(axis.text.x = element_text(angle=30, hjust = 1)) + 
  scale_fill_manual(values = c("#1E141E", "#B980AA", "#8A4F8F", "#522A5C", "#D2D2D7", "#8F908F", "#737373"))
dev.off()
############################################################
##########################################################################################


# Intergenic methylation in Transposon not-Transposon or in non-Repetive regions (Fig.S9C)
##########################################################################################
library("genomation")
library("GenomicRanges")
GR_intergenic_in_Trasposon <- subsetByOverlaps(GR_intergenic, RE_Transposons_GR) # Intergenic TEs
GR_intergenic_in_Not_Trasposon <- subsetByOverlaps(GR_intergenic, RE_Not_Transposons_GR) # Intergenic non-TEs Repeats
GR_intergenic_NOTin_AllRepeats <- subsetByOverlaps(GR_intergenic, AllRepeats_GR, invert=TRUE) # Intergenic non-repetitive

# Intergenic locus in Transposon for Hatchling RRBS
########################################################################
ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon = ScoreMatrixBin(target=GR_CpGs_All,
                                                                windows=resize(GR_intergenic_in_Trasposon, width=width(GR_intergenic_in_Trasposon)+1),
                                                                bin.num=15, strand.aware=TRUE,
                                                                weight.col="perc.meth",
                                                                is.noCovNA = TRUE)

dir.create("./Results/Genomation/MetaPlot/RE_intergenic")

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################

# Intergenic locus in Not_Transposon (all other RE) for Hatchling RRBS
########################################################################
ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon = ScoreMatrixBin(target=GR_CpGs_All,
                                                                    windows=resize(GR_intergenic_in_Not_Trasposon, width=width(GR_intergenic_in_Not_Trasposon)+1),
                                                                    bin.num=15, strand.aware=TRUE,
                                                                    weight.col="perc.meth",
                                                                    is.noCovNA = TRUE)

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################

# Intergenic locus NOT in AllRepeats for Hatchling RRBS
########################################################################
# CpGs_All_AllRepeats_GR has score column instead of perc.meth
ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats = ScoreMatrixBin(target=GR_CpGs_All,
                                                                    windows=resize(GR_intergenic_NOTin_AllRepeats, width=width(GR_intergenic_NOTin_AllRepeats)+1),
                                                                    bin.num=15, strand.aware=TRUE,
                                                                    weight.col="perc.meth",
                                                                    is.noCovNA = TRUE)

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################

# Intergenic locus in Transposon for SupraE WGBS
########################################################################
ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon = ScoreMatrixBin(target=GR_Octopus_WGBS_Brain_SupraE_select,
                                                                windows=resize(GR_intergenic_in_Trasposon, width=width(GR_intergenic_in_Trasposon)+1),
                                                                bin.num=15, strand.aware=TRUE,
                                                                weight.col="perc.meth",
                                                                is.noCovNA = TRUE)

dir.create("./Results/Genomation/MetaPlot/RE_intergenic")

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################

# Intergenic locus in Not_Transposon (all other RE) for SupraE WGBS
########################################################################
ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon = ScoreMatrixBin(target=GR_Octopus_WGBS_Brain_SupraE_select,
                                                                    windows=resize(GR_intergenic_in_Not_Trasposon, width=width(GR_intergenic_in_Not_Trasposon)+1),
                                                                    bin.num=15, strand.aware=TRUE,
                                                                    weight.col="perc.meth",
                                                                    is.noCovNA = TRUE)

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################

# Intergenic locus NOT in AllRepeats for SupraE WGBS
########################################################################
ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats = ScoreMatrixBin(target=GR_Octopus_WGBS_Brain_SupraE_select,
                                                                    windows=resize(GR_intergenic_NOTin_AllRepeats, width=width(GR_intergenic_NOTin_AllRepeats)+1),
                                                                    bin.num=15, strand.aware=TRUE,
                                                                    weight.col="perc.meth",
                                                                    is.noCovNA = TRUE)

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################

# Intergenic locus in Transposon for SubE WGBS
########################################################################
ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon = ScoreMatrixBin(target=GR_Octopus_WGBS_Brain_SubE_select,
                                                                windows=resize(GR_intergenic_in_Trasposon, width=width(GR_intergenic_in_Trasposon)+1),
                                                                bin.num=15, strand.aware=TRUE,
                                                                weight.col="perc.meth",
                                                                is.noCovNA = TRUE)

dir.create("./Results/Genomation/MetaPlot/RE_intergenic")

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_intergenic_in_Trasposon), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################

# Intergenic locus in Not_Transposon (all other RE) for SubE WGBS
########################################################################
ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon = ScoreMatrixBin(target=GR_Octopus_WGBS_Brain_SubE_select,
                                                                    windows=resize(GR_intergenic_in_Not_Trasposon, width=width(GR_intergenic_in_Not_Trasposon)+1),
                                                                    bin.num=15, strand.aware=TRUE,
                                                                    weight.col="perc.meth",
                                                                    is.noCovNA = TRUE)

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_intergenic_in_Not_Trasposon), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################

# Intergenic locus NOT in AllRepeats for SubE WGBS
########################################################################
ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats = ScoreMatrixBin(target=GR_Octopus_WGBS_Brain_SubE_select,
                                                                    windows=resize(GR_intergenic_NOTin_AllRepeats, width=width(GR_intergenic_NOTin_AllRepeats)+1),
                                                                    bin.num=15, strand.aware=TRUE,
                                                                    weight.col="perc.meth",
                                                                    is.noCovNA = TRUE)

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 100), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_CpG_All_RE_GR_intergenic_NOTin_AllRepeats), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "DNA methylation score", xlab = "binned base")
dev.off()
########################################################################


# Intergenic CpGs density in Transposon not-Transposon or in non-Repetive regions (Fig.S9D)
##########################################################################################
# Intergenic Transposons 
########################################################################
ScoreBin_WG_CpG_GR_intergenic_in_Trasposon = ScoreMatrixBin(target=cpgr,
                                                            windows=resize(GR_intergenic_in_Trasposon, width=width(GR_intergenic_in_Trasposon)+1),
                                                            bin.num=15, strand.aware=TRUE
)

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_WG_CpG_GR_intergenic_in_Trasposon.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_WG_CpG_GR_intergenic_in_Trasposon,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 0.06), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_WG_CpG_GR_intergenic_in_Trasposon), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "CpG density (Avg. score)", xlab = "binned base")
dev.off()
########################################################################

# Intergenic not-Transposons (all other Repeats) 
########################################################################
ScoreBin_WG_CpG_GR_intergenic_in_Not_Trasposon = ScoreMatrixBin(target=cpgr,
                                                                windows=resize(GR_intergenic_in_Not_Trasposon, width=width(GR_intergenic_in_Not_Trasposon)+1),
                                                                bin.num=15, strand.aware=TRUE
)

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_WG_CpG_GR_intergenic_in_Not_Trasposon.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_WG_CpG_GR_intergenic_in_Not_Trasposon,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 0.06), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_WG_CpG_GR_intergenic_in_Not_Trasposon), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "CpG density (Avg. score)", xlab = "binned base")
dev.off()
########################################################################

# Intergenic NOT in AllRepeats 
########################################################################
ScoreBin_WG_CpG_GR_intergenic_NOTin_AllRepeats = ScoreMatrixBin(target=cpgr,
                                                                windows=resize(GR_intergenic_NOTin_AllRepeats, width=width(GR_intergenic_NOTin_AllRepeats)+1),
                                                                bin.num=15, strand.aware=TRUE
)

pdf("./Results/Genomation/MetaPlot/RE_intergenic/ScoreBin_WG_CpG_GR_intergenic_NOTin_AllRepeats.pdf", onefile = FALSE, paper = "special", width = 10, height = 7.5)
plotMeta(ScoreBin_WG_CpG_GR_intergenic_NOTin_AllRepeats,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 0.06), dispersion="se",
         dispersion.col = rainbow(length(ScoreBin_WG_CpG_GR_intergenic_NOTin_AllRepeats), start = 4/6, alpha = 0.1),
         winsorize=c(0,99),line.col="black", lwd=2,
         ylab = "CpG density (Avg. score)", xlab = "binned base")
dev.off()
########################################################################

