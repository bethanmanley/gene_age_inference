
library('ggplot2')
library(stringr)
library("GenomicRanges")
library("GenomicFeatures")
library(rtracklayer)
library(chromPlot)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(tibble)
library(tidyr)

  
################################################################################
### Gene count distribution ( n of genes per phylorank)
################################################################################

PS_filt2 <- PS_filt 
PS_filt2$transcript_id <- gsub("-.*", "", PS_filt2$transcript_id)
PS_filt2 <- PS_filt2 %>% 
  group_by(transcript_id) %>% 
  filter(PS == min(PS)) %>% 
  distinct

#Summarise n of genes 
PS_sum <- PS_filt2 %>% 
  group_by(PS) %>% 
  summarise(Count = n())

#Plot
PS_sum %>% 
  ggplot(data=., aes(x=PS, y=Count, width=.5)) +
  geom_bar(stat="identity", fill="black") +
  scale_x_discrete(limits=1:9) +
  coord_flip() +
  theme_classic() +
  xlab("") +
  ylab("Number of genes") 



################################################################################
### Top 15 InterPro domains
################################################################################
library(dplyr)
library(tidytext)
library(textrank)
library(ggrepel)

GO <- read.table("Rhizophagus_irregularis_DAOM197198.annotations.txt", header = T, sep="\t", stringsAsFactors=FALSE, quote="")
GO <- GO %>% mutate_all(na_if,"")

#Add Gene annotation
names(GO)[2] <- "transcript_id"

PS_j_info <- left_join(PS_j, GO, by=c("transcript_id"))
PS_j_info <- PS_j_info[!is.na(PS_j_info$seqnames),]

#List of Glomeromycetes-restricted genes (Table S2)
PS_j_hd_info <- left_join(PS_j_hd, GO, by=c("transcript_id"))
PS_j_hd_info_PS6 <- PS_j_hd_info %>%
  filter(PS == 6)
PS_j_hd_info_PS6 <- PS_j_hd_info %>%
  filter(PS == 6)
PS_j_hd_info_PS6 <- PS_j_hd_info_PS6 %>%
  select(PFAM, InterPro, GO.Terms, Secreted, antiSMASH) %>%
  arrange(InterPro)
write.table(PS_j_hd_info_PS6, paste(out_dir,"TableS2_HDF_phylorank6_genes.tsv", sep="/"), col.names=T, quote=F, sep="\t", row.names=F)

head(PS_j_info)
colnames(PS_j_info)
PS_j_info$gDNA <- NULL
PS_j_info$mRNA <- NULL
PS_j_info$CDS.transcript <- NULL
PS_j_info$Translation <- NULL

#Extract words from InterPro column
Data <- PS_j_info %>% 
  group_by(PS) %>% 
  mutate(InterPro = paste0(InterPro, collapse = ";")) %>%
  select(PS, InterPro) %>% 
  distinct()

#Tokenise the InterPro words
myTokens <- Data %>% 
  unnest_tokens(word, InterPro) %>% 
  anti_join(stop_words, by = "word")

#Keep only words that start with 'ipr'
myTokens <- myTokens %>%
  filter(str_detect(word, "^ipr"))

#Count occurence of each IPR domain, grouped by PS
counts <- myTokens %>% 
  group_by(PS, word) %>% 
  summarize(count=n())

#Keep top 15 for each PS group
top15 <- counts %>%
  arrange(desc(count)) %>%
  group_by(PS) %>% 
  slice(1:15)

write.table(top15, paste(out_dir,"top15.PS.IPRs.txt", sep="/"), col.names=T, quote=F, sep="\t", row.names=F)



################################################################################
### Plot PHYLOSTRATA on chromosomes (LINE PLOT)
################################################################################

#Make cytoband dfs
PS_cyto <- data.frame(PS_j$seqnames, PS_j$start, PS_j$end, PS_j$gene_id, PS_j$PS)
colnames(PS_cyto) <- c("Chrom", "Start", "End", "Name", "Group") 
PS_cyto <- PS_cyto[!is.na(PS_cyto$Chrom), ]
PS_cyto <- PS_cyto[!is.na(PS_cyto$Group), ]
PS_cyto$Group <- as.numeric(PS_cyto$Group)
# Line plot with colour-scaled line
library(gam)
library(mgcv)
library(broom)
library(scales)
sapply(PS_cyto, class)
#Define model from smooth fit
model <- mgcv::gam(Group ~ s(Start, bs = "cs", by=Chrom), data=PS_cyto)
#Predict fitted values
predicted.values <- data.frame(predict(model, by=Chrom, type = "response"))
#Add modelled values to original dataframe
PS_cyto$gam <- predict(model, by=Chrom, type = "response")
sapply(PS_cyto,class)
PS_cyto$gam <- as.numeric(PS_cyto$gam)

#Plot with colours gradient
ggplot(data = PS_cyto, aes(x=Start, y=gam)) +
  geom_line(aes(color = gam), lineend = "round") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_gradientn(colours = c("#003f5c", "#2f4b7c", "#665191", "#a05195", "#d45087", "#f95d6a", "#ff7c43", "#ffa600"), 
                        values = rescale(x = c(1, 2, 3, 4, 5, 6, 7, 8,9), from = c(1, 9))) +
  facet_wrap(~Chrom, ncol=6)

PS_cyto$value <- NULL
write.table(PS_cyto, paste(out_dir,"geneages.txt", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

# The mating locus is at chr.30-:2,120,246-2,121,783
#KU128850.1 Rhizophagus irregularis isolate DAOM-197198 H+ ATPase gene, partial cds
#unnamed	30+	100.00%	1,538	100.00%	0	0	1	1,538	2,120,246	2,121,783	0	2841

#Make df where we have 1000 random samplings of gene age
test <- reduce(seq(0, 1000, by = 1), .init = PS_cyto,
               ~ mutate(.x, !!paste0('x', .y) := sample(Group)))
PS_cyto$mean_shuffle = rowMeans(select(test, x0:x1000))
#Plot the mean of a 1000 permutations: this is what the distributions would look like of gene age was randomly distributed.
ggplot(data = PS_cyto, aes(x=Start, y=mean_shuffle)) +
  geom_line(color = "grey76", lineend = "round") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(3.5,5.5)+
  facet_wrap(~Chrom, ncol=6)

#Model 1000x shuffled gene ages
shuf_model <- mgcv::gam(mean_shuffle ~ s(Start, bs = "cs", by=Chrom), data=PS_cyto)
shuf_predicted.values <- data.frame(predict(shuf_model, by=Chrom, type = "response"))
PS_cyto$gam.shuf <- predict(shuf_model, by=Chrom, type = "response")
PS_cyto$gam.shuf <- as.numeric(PS_cyto$gam.shuf)
#Plot the FIT of a 1000 permutations
ggplot(data = PS_cyto, aes(x=Start, y=gam.shuf)) +
  geom_line(color = "grey76", lineend = "round") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(3.5,5.5)+
  facet_wrap(~Chrom, ncol=6)

#write.table(PS_cyto_shuffle, paste(out_dir,"shuffled_geneages.txt", sep="/"), col.names=F, quote=F, sep="\t", row.names=F)

#T-test to check if gene age distributions are significantly different from random
df1 <- PS_cyto[, c("Chrom", "Start", "gam")]
df2 <- PS_cyto[, c("Chrom", "Start", "gam.shuf")]
df2 <- PS_cyto[, c("Chrom", "Start", "mean_shuffle")]

names(df2)[3] <- "gam"
df1$value <- "obs"
df2$value <- "shuf"
hey <- rbind(df1, df2)

library(tidyverse)
library(ggpubr)
library(rstatix)
pwc <- hey %>%
  group_by(Chrom) %>%
  pairwise_t_test(gam ~ value, paired = TRUE, p.adjust.method = "bonferroni")
pwc



################################################################################
### Plot SMALL RNA loci on IDEOGRAM
################################################################################

setwd("/Users/alexandradallaire/genomics/Manley_2022/sRNA/ShortStack_1644398621")
results <- read.table("Results.txt", header = FALSE, check.names = FALSE)
sRNAgff <- import.gff("ShortStack_D.gff3")
sRNAgff <- data.frame(sRNAgff)
colnames(sRNAgff) <- c("Seqname", "Start", "End", "width",  "Strand", "Source", "Type", "Score", "Phase", "ID", "DicerCall") 
anno_file_sRNA <- data.frame(sRNAgff$Seqname, sRNAgff$Source, sRNAgff$Type, sRNAgff$Start, sRNAgff$End, sRNAgff$Strand, sRNAgff$Phase, sRNAgff$ID)
colnames(anno_file_sRNA) <- c("Chr", "Source", "Feature", "Start", "End", "Strand", "Frame", "ID") 
names(results)[2] <- "ID"
sRNAgff2 <- left_join(results, anno_file_sRNA, by=c("ID"))
#Make new gff that includes sRNA RPM (column V5)
sRNAgff3 <- data.frame(sRNAgff2$Chr, sRNAgff2$Source, sRNAgff2$Feature, sRNAgff2$Start, sRNAgff2$End, sRNAgff2$V5, sRNAgff2$Strand, sRNAgff2$Frame, sRNAgff2$ID)
colnames(sRNAgff3) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "ID") 
#anno_file_sRNA3$V5 = log10(anno_file_sRNA3$V5)
cytoband <- data.frame(sRNAgff2$Chr, sRNAgff2$Start, sRNAgff2$End, sRNAgff2$Chr)
colnames(cytoband) <- c("Chrom", "Start", "End", "Seqname") 
cytoband$Colors = "red"
#cytoband$End = cytoband$Start+5000

chromPlot(gaps=chr_file, bands=cytoband, chr=c(1:32), figCols=5)
# PS line + sRNA bands
chromPlot(gaps=chr_file, bands=cytoband, stat=PS_cyto, statCol="Group", statName="Mean PS", statTyp="l", statSumm="mean", scex=0.01, spty=1, chr=c(1:32), figCols=5)
# DS line + sRNA bands
chromPlot(gaps=chr_file, bands=cytoband, stat=DS_cyto, statCol="Group", statName="Mean DS", statTyp="l", statSumm="mean", scex=0.01, spty=1, chr=c(1:32), figCols=5)
View(cytoband)
sapply(cytoband,class)



########################################
# Plot small RNA coverage across chromosomes (normalised RPKM (LINEPLOT))
########################################
samtools sort ~/genomics/rirman22/sRNA/ShortStack_1644398621/Rhiir_sRNA_mix3.bam -o ~/genomics/rirman22/sRNA/ShortStack_1644398621/Rhiir_sRNA_mix3_sorted.bam
samtools index ~/genomics/rirman22/sRNA/ShortStack_1644398621/Rhiir_sRNA_mix3_sorted.bam
bamCoverage --bam ~/genomics/rirman22/sRNA/ShortStack_1644398621/Rhiir_sRNA_mix3_sorted.bam --outFileName ~/genomics/rirman22/sRNA/ShortStack_1644398621/Rhiir_sRNA_mix3_sorted_normcoverage.bedgraph --outFileFormat bedgraph --binSize 200 --ignoreDuplicates --normalizeUsing RPKM -p 12

setwd("/Users/alexandradallaire/genomics/Manley_2022/sRNA/ShortStack_1644398621")
bamCoverage <- read.table("Rhiir_sRNA_mix3_sorted_normcoverage.bedgraph", header = F, sep="\t", stringsAsFactors=FALSE, quote="")
colnames(bamCoverage) <- c("Chrom", "Start", "End", "Cov") 
bamCoverage$Chrom <- as.factor(bamCoverage$Chrom)

#Define model from smooth fit
model2 <- mgcv::gam(Cov ~ s(Start, bs = "cs", by=Chrom), unconditional=FALSE, data=bamCoverage)
#Add modelled values to original dataframe
bamCoverage$gam <- predict(model2, by=Chrom, type = "response")
sapply(bamCoverage,class)
bamCoverage$gam <- as.numeric(bamCoverage$gam)
#Plot with colours gradient
ggplot(data = bamCoverage, aes(x=Start, y=gam)) +
  geom_line(aes(color = gam), lineend = "round") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_gradientn(colours = c("#440154FF", "#46337EFF", "#365C8DFF", "#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF"), 
                        values = rescale(x = c(1, 2, 3, 4, 5, 6, 8), from = c(1, 8))) +
  facet_wrap(~Chrom, ncol=6)

# Correlation between genome-wide sRNA levels and PS
library(data.table)
names(bamCoverage)[5] <- "gamcov"
names(bamCoverage)[1] <- "Chrom"
names(bamCoverage)[2] <- "Start2"
# rolling join:  nearest RPKM value for each gene (based on Start position)
dt1 <- as.data.table(PS_cyto)
dt2 <- as.data.table(bamCoverage)
dt1 <- dt1 %>% arrange(Chrom, Start)
dt2 <- dt2 %>% arrange(Chrom, Start2)

bind2 <- dt2[, gam := gamcov][PS_cyto, on = .(Chrom, Start2 = Start), roll = "nearest"]
bind_subsample = bind2[seq(1, nrow(bind2), 50), ]
ggplot(bind_subsample, aes(x=i.gam, y=gamcov)) + 
  geom_point(alpha=0.5) +
  stat_smooth(method=lm, se = TRUE) +
  stat_cor(method = "pearson", label.x = 2, label.y = 2) +
  ggtitle("1:50nth subsample") +
  labs(color='Chromosome') +
  xlab("PS") +
  ylab("sRNAcov") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 8, by = 1))
  

#Next: Correlation between smallRNAseq reads and PS, at the GENE-LEVEL (instead of all-encompassing genomic coverage)
########################################
#Small RNA counts PER GENE (counting based on the gene annotation)
########################################

featureCounts -M \
--fraction \
-T 8 \
-F GTF \
-g gene_id  \
-t transcript \
-a ~/genomics/rirman22/protein/funct_annot_postupdate/annotate_results/Rhizophagus_irregularis_DAOM197198_Illumina+ONT_curated.gtf \
-o ~/genomics/rirman22/sRNA/ShortStack_1644398621/sRNA_mix3_FeatureCounts_genecounts.txt ~/genomics/rirman22/sRNA/ShortStack_1644398621/Rhiir_sRNA_mix3_sorted.bam

#Load small RNA featureCounts output
setwd("/Users/alexandradallaire/genomics/Manley_2022/sRNA/ShortStack_1644398621")
sRNAcounts <- read.table("sRNA_mix3_FeatureCounts_genecounts.txt", header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="")
colnames(sRNAcounts) <- c("Name", "Chrom", "Start", "End", "Strand", "Length", "Count") 
#Reads per million
sum(sRNAcounts$Count) # 14540513
sRNAcounts$RPM = sRNAcounts$Count/14540513*1000000
sRNAcounts$RPKM = sRNAcounts$RPM/sRNAcounts$Length*1000

#sRNA RPKM+1 PER GENE
sRNA_genecounts <- left_join(sRNAcounts, PS_cyto, by=c("Name"))
sRNA_genecounts %>% 
  filter(!is.na(Group)) %>% 
  ggplot(., aes(x=Group, group = Group, y=log(RPKM+1))) + 
  geom_violin(aes(fill = Group)) + 
  geom_boxplot(width = 0.3, alpha=0.2, outlier.size = 0.1) +
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  theme_classic() +
  ylim(0,11) +
  ylab("log(small RNA RPKM+1)") +
  xlab("Phylostratum")+
  theme(legend.position="none")
#reg line
ggplot(sRNA_genecounts, aes(x=Group, y=log(RPKM+1))) +
  stat_smooth(method=lm, se = TRUE) +
  stat_cor(method = "pearson", label.x = 4, label.y = 1) +
  ggtitle("sRNA RPKM, linear regression") +
  xlab("PS") +
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  ylim(0,11) +
  theme_classic() +
  theme(legend.position="none")



################################################################################
### Plot long-read RNAseq RPKM (LINEPLOT)
################################################################################

samtools index nano_cDNA.coordSorted.bam
bamCoverage --bam nano_cDNA.coordSorted.bam --outFileName nano_cDNA_coverage.bedgraph --outFileFormat bedgraph --binSize 200 --ignoreDuplicates --normalizeUsing RPKM -p 12
setwd("/Users/alexandradallaire/genomics/Manley_2022/transcript_assembly/genes")
ONT_bamCoverage <- read.table("nano_cDNA_coverage.bedgraph", header = F, sep="\t", stringsAsFactors=FALSE, quote="")
colnames(ONT_bamCoverage) <- c("Chrom", "Start", "End", "Cov") 
ONT_bamCoverage$Chrom <- as.factor(ONT_bamCoverage$Chrom)

#Define model
model3 <- mgcv::gam(Cov ~ s(Start, bs = "cs", by=Chrom), unconditional=FALSE, data=ONT_bamCoverage)
#Add modelled values to original dataframe
ONT_bamCoverage$gam <- predict(model3, by=Chrom, type = "response")
ONT_bamCoverage$gam <- as.numeric(ONT_bamCoverage$gam)
#Plot with colours gradient
ggplot(data = ONT_bamCoverage, aes(x=Start, y=gam)) +
  geom_line(aes(color = gam), lineend = "round") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_gradientn(colours = c("#440154FF", "#46337EFF", "#365C8DFF", "#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF"), 
                        values = rescale(x = c(1, 2, 3, 4, 5, 6, 8), from = c(1, 8))) +
  facet_wrap(~Chrom, ncol=6)

#Correlation between long RNAseq coverage and PS
library(data.table)
names(ONT_bamCoverage)[5] <- "gamcov"
names(ONT_bamCoverage)[1] <- "Chrom"
names(ONT_bamCoverage)[2] <- "Start2"

# rolling join: keep nearest coverage value for each gene (based on Start position)
dt1 <- as.data.table(PS_cyto)
dt2 <- as.data.table(ONT_bamCoverage)
dt1 <- dt1 %>% arrange(Chrom, Start)
dt2 <- dt2 %>% arrange(Chrom, Start2)

bind2 <- dt2[, gam := gamcov][PS_cyto, on = .(Chrom, Start2 = Start), roll = "nearest"]

#Autocorrelation function
#Set order
bind2 <- bind2 %>% arrange(Chrom, Start2)
acf(bind2$i.gam, lag.max = 10000, type = c("correlation"), plot = TRUE) 
acf(bind2$gamcov, lag.max = 10000, type = c("correlation"),plot = TRUE)


#Scatter plot fit values of PS and coverage for every gene, with linear regression line
ggplot(bind2, aes(x=i.gam, y=gamcov)) + 
  geom_point(alpha=0.07, mapping=aes(color=Chrom)) +
  stat_smooth(method=lm, se = TRUE) +
  stat_cor(method = "pearson", label.x = 2, label.y = 72) +
  ggtitle("Smooth fit values (method=gam, formula=y~s(x,bs=cs)") +
  labs(color='Chromosome') +
  xlab("PS") +
  ylab("ONT RNAseq coverage") +
  scale_y_continuous(breaks = seq(0, 80, by = 20)) +
  scale_x_continuous(breaks = seq(1, 8, by = 1)) +
  theme_classic()
#Same but sample 1 in 200th value (based on acf plots)
bind2_subsample = bind2[seq(1, nrow(bind2), 50), ]
ggplot(bind2_subsample, aes(x=i.gam, y=gamcov)) + 
  geom_point(alpha=0.4) +
  stat_smooth(method=lm, se = TRUE) +
  stat_cor(method = "pearson", label.x = 2, label.y = 72) +
  ggtitle("1:50th subsample") +
  xlab("PS") +
  ylab("ONT RNAseq coverage") +
  scale_y_continuous(breaks = seq(0, 80, by = 20)) +
  scale_x_continuous(breaks = seq(1, 8, by = 1)) +
  theme_classic()

#Bar plot modeled values of PS and coverage for every gene, with linear regression line
ggplot(bind2, aes(x = Group, y = gamcov, color = Group)) +
  geom_jitter(alpha=0.05) +
  geom_smooth(method = "lm", se =TRUE) +
  stat_cor(method = "pearson", label.x = 0.5, label.y = 70) +
  labs(color='PS') +
  xlab("PS") +
  ylab("ONT RNAseq coverage") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 80, by = 20)) +
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  theme(legend.position="none")

#plot violin (raw values, without regression)
bind2$Group <- as.factor(bind2$Group)
ggplot(bind2, aes(x=Group, y=log(Cov+1), group=Group)) + 
  geom_violin(aes(fill = Group)) + 
  geom_boxplot(width = 0.3, alpha=0.2, outlier.size = 0.1) +
  xlab("PS") +
 # ylim(0,100) +
  ylab("ONT RNAseq coverage") +
  theme_classic()
#Regression raw values
bind2$Cov_2 <- (log((bind2$Cov)+1))
bind2$Group <- as.numeric(bind2$Group)
ggplot(bind2, aes(x=Group, y=Cov_2)) + 
  stat_smooth(method=lm, se = TRUE) +
  stat_cor(method = "pearson", label.x = 0.5, label.y = 3)+
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  ylim(0,6) +
  xlab("PS") +
  ylab("ONT RNAseq") +
  theme_classic() 



################################################################################
### CG methylation
################################################################################

meth <- read.table("Rirr.manley.fast5s.CG.call_mods.frequency.tsv")
options("scipen"=100)
# Take 1 off every row that contains a "-". This way, we get 2 lines with the same value (coordinate) in column 2 for those that have 1bp difference (same meth site, symmetrical)
meth[meth$V3=="-", 2] = meth[meth$V3=="-", 2] -1
# Aggregate rows that have the same V1 and V2 value (coordinate),  mean V10 values (methylation score). This reduces the n of datapoints by half
values = aggregate(meth$V10, by = list(meth$V1, meth$V2), FUN=mean)
head(values)
#Plot sRNA RPM as dots and TEs as bands
cytoband.meth <- data.frame(values$Group.1, values$Group.2, values$Group.2, values$x)
colnames(cytoband.meth) <- c("Chrom", "Start", "End", "meth.score") 
View(cytoband.meth)
chromPlot(gaps=chr_file,stat=cytoband.meth, statCol="meth.score", statName="methylation score", statTyp="p", scex=0.01, spty=1, chr=c(1:32), figCols=5)

cytoband.highmeth <- cytoband.meth %>% filter(meth.score >= .90)
sapply(cytoband.highmeth,class)
cytoband.highmeth$Chrom=as.factor(cytoband.highmeth$Chrom)
sapply(cytoband.highmeth,class)
cytoband.highmeth$End = cytoband.highmeth$Start+50
cytoband.highmeth$meth.score = cytoband.highmeth$Chrom
colnames(cytoband.highmeth) <- c("Chrom", "Start", "End", "Seqname") 
cytoband.highmeth$Colors <- "blue"
cytoband.highmeth$Start=as.integer(cytoband.highmeth$Start)
cytoband.highmeth$End=as.integer(cytoband.highmeth$End)
View(cytoband.highmeth)
cytoband.highmeth <- cytoband.highmeth[
  with(cytoband.highmeth, order(Chrom, Start)),
]
chromPlot(bands=cytoband.highmeth, chr=c(1:32), figCols=5)

#Make bed file for IGV
highly_meth <- meth %>% filter(V10 >= .80)
myvars <- c("V1", "V2", "V3", "V10")
highly_meth_f <- highly_meth[myvars]
head(highly_meth_f)
highly_meth_f$V3 = highly_meth_f$V2+15
head(highly_meth_f)

#File for IGV bundle
write.table(highly_meth_f, paste(out_dir,"high_meth.bed",sep="/"), col.names=F, quote=F, sep="\t", row.names=F)



################################################################################
### Phylstratigraphy of other Mucoromycota fungi
################################################################################

#Filtering and correcting phyloranks
PS <- read_tsv("*_phylostrata_assignation.tsv")

# Correcting phyloranks:
# PS rank at the “Fungi incertae sedis” level (basically a wastebasket group) is moved to the kingdom level (Fungi)
PS_filt <- PS %>% mutate(rank = replace(rank, rank %in% c(5), 4)) %>%
  mutate(rank = replace(rank, rank == 6,5)) %>%
  mutate(rank = replace(rank, rank == 7,6)) %>%
  mutate(rank = replace(rank, rank == 8,7)) %>%
  mutate(rank = replace(rank, rank == 9,8)) %>%
  mutate(rank = replace(rank, rank == 10,9)) %>%
  mutate(phylostratum = replace(phylostratum, phylostratum == "Fungi incertae sedis", "Fungi"))

# Strain-level rank moved to species level
PS_filt <- PS %>% mutate(rank = replace(V3, V3 == 9,8)) %>%
  mutate(phylostratum = replace(V2, V2 == "Linnemannia elongata AG-77", "Linnemannia elongata")) %>%
  select(V1, phylostratum, rank, V4) 

#Gene count distributions
#This is a file containing the number genes at each phylorank, for each species analysed
GC <- read.table(file = 'Gene_count_distribs.tsv', sep = '\t', header = T, fill=TRUE)
#Plot Figure S2A
ggplot(GC, aes(x = X.1, y = (X),
               group = X.2, color = X.2)) + 
  geom_point(size = 2, shape = 16) +
  geom_line(col = 'gray') +
  theme_classic() +  
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  facet_wrap(~X.2, ncol=3, scales='free_x') +
  ylab("log(number of genes)
       
       ") +
  theme(axis.text.x = element_text(angle = 50,  hjust = 1)) +
  theme(legend.position="none")




