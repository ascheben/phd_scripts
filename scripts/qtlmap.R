library(tidyverse)
library(reshape2)
library(LinkageMapView)
library(ASMap)
library(qtl)
library(ggplot2)
library(Hmisc)
library(corrplot)
library(tibble)
library(dplyr)
library(ggpubr)

# Function to set wd as script dir
set_wd <- function() {
  library(rstudioapi) # make sure you have it installed
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
  print( getwd() )
}
set_wd()

########################################
# Phenotype analysis                   #
########################################

# Fomat correlation matrix for plotting
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# Get phenotypes used for QTL mapping
qtl_pheno <- finalmap_f2$pheno
# Read in headered table of all individuals and phenotypes
mypheno <- read.csv("pheno_analysis_all.csv")
# Remove individual name column
pheno_obs <- mypheno[, 2:7]

# Comput matrix of Pearson's rank correlation coefficients
pcc <- rcorr(as.matrix(pheno_obs), type = "pearson")

# Convert matrix to four coloumn pairwise correlation coefficients and p-values
flattenCorrMatrix(pcc$r, pcc$P)

# Create correlation plot with alpha=0.05
corrplot(pcc$r, type="upper", order="hclust", 
         p.mat = pcc$P, sig.level = 0.05, insig = "blank")

# Shapiro-Wilkes normality test
shapiro.test(qtl_pheno$Flower)

# Plot phenotype histograms with ggplot
# Insert arrow pointing to position of parental samples
# use breaks to define bins
# add result of shapiro test with geom_text or geom_label

# The Bud value for sample 146 is '19', which is an extreme outlier
# The value was made NA for plotting purposes
# Change row 203, column 6 to NA
mypheno[203, 6] = NA

# Define vector of parent individuals to label histograms
parents <- c("S1","S2","S3","S4","W1","W2","W3")

# Germination
pdf(file = "Germination_hist_qtl.pdf")
shapiro_res <- shapiro.test(qtl_pheno$Germ)
test_label = paste("Shapiro-Wilkes p = ", round(shapiro_res$p.value, 4), sep = "")

p <- gghistogram(qtl_pheno, x = "Germ", 
            main = "Histogram of germination time",
            xlab = "Days to germination",
            ylab = "Individuals",
            fill = "lightgray", 
            bins = 30, 
            add_density = TRUE, 
            #label = "Number", 
            #label.select = parents, 
            xticks.by = 5,
            position = "stack"
            #label.rectangle = TRUE
            ) 
annotate_figure(p, fig.lab = test_label, fig.lab.pos = "bottom.right")
dev.off()  

# LeafD
pdf(file = "LeafD_hist_qtl.pdf")
shapiro_res <- shapiro.test(qtl_pheno$LeafD)
test_label = paste("Shapiro-Wilkes p = ", round(shapiro_res$p.value, 4), sep = "")

p <- gghistogram(qtl_pheno, x = "LeafD", 
                 main = "Histogram of leaf diameter",
                 xlab = "Leaf diameter (mm)",
                 ylab = "Individuals",
                 fill = "lightgray", 
                 bins = 30, 
                 add_density = TRUE, 
                 #label = "Number", 
                 #label.select = parents, 
                 xticks.by = 5,
                 position = "stack"
                 #label.rectangle = TRUE
) 
annotate_figure(p, fig.lab = test_label, fig.lab.pos = "bottom.right")
dev.off()  
# NLeaves
pdf(file = "NLeaves_hist_qtl.pdf")
shapiro_res <- shapiro.test(qtl_pheno$NLeaves)
test_label = paste("Shapiro-Wilkes p = ", round(shapiro_res$p.value, 4), sep = "")

p <- gghistogram(qtl_pheno, x = "NLeaves", 
                 main = "Histogram of leaf number at 5-leaf stage",
                 xlab = "Number of leaves",
                 ylab = "Individuals",
                 fill = "lightgray", 
                 bins = 30, 
                 add_density = TRUE, 
                 #label = "Number", 
                 #label.select = parents, 
                 xticks.by = 1,
                 position = "stack"
                 #label.rectangle = TRUE
) 
annotate_figure(p, fig.lab = test_label, fig.lab.pos = "bottom.right")
dev.off()  
# Hair
pdf(file = "Hair_hist_qtl.pdf")
shapiro_res <- shapiro.test(qtl_pheno$Hair)
test_label = paste("Shapiro-Wilkes p = ", round(shapiro_res$p.value, 4), sep = "")
p <- gghistogram(qtl_pheno, x = "Hair", 
                 main = "Histogram of hair",
                 xlab = "Hair presence or absence",
                 ylab = "Individuals",
                 fill = "lightgray", 
                 bins = 2, 
                 #label = "Number", 
                 #label.select = parents, 
                 xticks.by = 1,
                 position = "stack"
                 #binwidth = 5
                 #label.rectangle = TRUE
) 
annotate_figure(p, fig.lab = test_label, fig.lab.pos = "bottom.right")
dev.off()  
# Bud 
pdf(file = "Bud_hist_qtl.pdf")
shapiro_res <- shapiro.test(qtl_pheno$Bud)
test_label = paste("Shapiro-Wilkes p = ", round(shapiro_res$p.value, 4), sep = "")

p <- gghistogram(qtl_pheno, x = "Bud", 
                 main = "Histogram of time to bud",
                 xlab = "Days to first floral bud",
                 ylab = "Individuals",
                 fill = "lightgray", 
                 bins = 30, 
                 add_density = TRUE, 
                 #label = "Number", 
                 #label.select = parents, 
                 xticks.by = 2,
                 position = "stack"
                 #label.rectangle = TRUE
) 
annotate_figure(p, fig.lab = test_label, fig.lab.pos = "bottom.right")
dev.off()  

# Flower
pdf(file = "Flower_hist_qtl.pdf")
shapiro_res <- shapiro.test(qtl_pheno$Flower)
test_label = paste("Shapiro-Wilkes p = ", round(shapiro_res$p.value, 4), sep = "")

parents <- c("S1","S2","S3","S4","W1","W2","W3")
p <- gghistogram(qtl_pheno, x = "Flower", 
                 main = "Histogram of flowering time",
                 xlab = "Days to flower",
                 ylab = "Individuals",
                 fill = "lightgray", 
                 bins = 30, 
                 add_density = TRUE, 
                 #label = "Number", 
                 #label.select = parents, 
                 xticks.by = 5,
                 position = "stack"
                 #label.rectangle = TRUE
) 
annotate_figure(p, fig.lab = test_label, fig.lab.pos = "bottom.right")
dev.off()

# qqplot
ggqqplot(pheno_obs$Flower)


####################################
# Genotype error analysis          #
####################################

# Pairwise comparison of filtered SNPs was carried out for duplicate samples using bash
errfrq <- read.csv("Pairwise_genotype_congruence_duplicate_samples_merged.csv")

# Get mean values of all S samples and of all W samples
errfrqm_sw <- data.frame(Type=errfrq[,1], SMeans=rowMeans(errfrq[2:7]),WMeans=rowMeans(errfrq[8:10]))
errfrqm  <- data.frame(Type=errfrq[,1], Means=rowMeans(errfrq[,-1]))
# Melt table for plotting
errfrqm_sw <- melt(errfrqm_sw, id = "Type")
errfrq <- melt(errfrq, id = "Type")
# Set custom group order
errfrqm_sw$Type <- factor(errfrqm_sw$Type, levels = c("miss_sum","match_sum","./._./.","./._0/0","./._0/1","./._1/1","0/0_0/0","1/1_1/1","0/1_0/1","0/0_1/1","0/1_1/1","0/0_0/1"))
errfrqm$Type <- factor(errfrqm$Type, levels = c("miss_sum","match_sum","./._./.","./._0/0","./._0/1","./._1/1","0/0_0/0","1/1_1/1","0/1_0/1","0/0_1/1","0/1_1/1","0/0_0/1"))
errfrq$Type <- factor(errfrq$Type, levels = c("miss_sum","match_sum","./._./.","./._0/0","./._0/1","./._1/1","0/0_0/0","1/1_1/1","0/1_0/1","0/0_1/1","0/1_1/1","0/0_0/1"))
# Bar plot with all pairwise comparisons
ggplot(data = subset(errfrq,!(Type %in% c("miss_sum","match_sum"))), aes(x = Type, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'dodge') +  ylab("Frequency") + xlab("Genotype concordance type") +
  scale_fill_discrete(name = "Pairwise comparisons")
# Bar plot with all pairwise comparisons using proportion instead of frq
ggplot(data = subset(errfrq,!(Type %in% c("miss_sum","match_sum"))), aes(x = Type, y = value / 124804, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'dodge') +  ylab("Fraction of total genotypes") + xlab("Genotype concordance type") +
  scale_fill_discrete(name = "Pairwise comparisons")
# Plot using means
ggplot(data = subset(errfrqm,!(Type %in% c("miss_sum","match_sum"))), aes(x = Type, y = Means)) + 
  geom_bar(stat = 'identity', position = 'dodge') +  ylab("Frequency") + xlab("Genotype concordance type")
# Plot using means of S and W
ggplot(data = subset(errfrqm_sw,!(Type %in% c("miss_sum","match_sum"))), aes(x = Type, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'dodge') +  ylab("Frequency") + xlab("Genotype concordance type")
# Plot using means of S and W using fractions
ggplot(data = subset(errfrqm_sw,!(Type %in% c("miss_sum","match_sum"))), aes(x = Type, y = value / 124804, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'dodge') +  ylab("Fraction of total genotypes") + xlab("Genotype concordance type")
# Plot using means of S and W using fractions - simplified
ggplot(data = subset(errfrqm_sw,!(Type %in% c("./._./.","./._0/0","./._0/1","./._1/1","0/0_0/0","1/1_1/1","0/1_0/1"))), aes(x = Type, y = value / 124804, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'dodge') +  ylab("Fraction of total genotypes") + xlab("Genotype concordance type")

# Plot excluding types with missing
ggplot(data = subset(errfrq, !(Type %in% c("./._./.","./._0/0","./._0/1","./._1/1"))), aes(x = Type, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'dodge') +  ylab("Frequency") + xlab("Genotype concordance type") +
  scale_fill_discrete(name = "Pairwise comparisons")
ggplot(data = subset(errfrqm, !(Type %in% c("./._./.","./._0/0","./._0/1","./._1/1"))), aes(x = Type, y = Means)) + 
  geom_bar(stat = 'identity', position = 'dodge') +  ylab("Frequency") + xlab("Genotype concordance type")
ggplot(data = subset(errfrqm_sw, !(Type %in% c("./._./.","./._0/0","./._0/1","./._1/1"))), aes(x = Type, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'dodge') +  ylab("Frequency") + xlab("Genotype concordance type")

# Heterozygosity was calculated using vcftools 

hetstats <- read.table("Pop3_ddrad_raw_variants_bwa_v81_mim09_biallelic_minDP5_mm08_maf005.recode_mim05.het", sep = "\t", header = TRUE)
# Get fraction homozygous and add group column
hetfrq <- data.frame(INDV=hetstats[,1], FractionHom=hetstats$O.HOM./hetstats$N_SITES, Group="F2")
# Change parental samples in group column based on individual name
hetfrq$Group <- ifelse(grepl('S', hetfrq$INDV), 'Parent_S',
                       ifelse(grepl('W', hetfrq$INDV), 'Parent_W', 'F2'))

ggplot(data = hetfrq, aes(x = INDV, y = FractionHom,fill = Group)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=3)) +  ylab("Fraction homozygous genotypes") + xlab("Individuals")
ggsave("FractionHomozygous_Pop3_ddrad_filt_bwa_v81_barplot.pdf")


# Genotype correction statistics
mycols <- c("black","coral","seagreen","cornflowerblue")
#ID <- c("AA","AB","AH","A-","BB","BA","BH","B-","HH","HA","HB","H-","-A","-B","-H","--")
ChangedTo <- c("A","B","H","-","B","A","H","-","H","A","B","-","A","B","H","-")
UncorrectedGenotype <- c("A","A","A","A","B","B","B","B","H","H","H","H","-","-","-","-")
# Old numbers (for incorrect data)
#Frequency <- c(468810,1135,37169,0,441338,1229,23405,0,935905,10777,14492,0,28120,34256,62918,6334)
gc_res <- data.frame(ChangedTo,UncorrectedGenotype,Frequency)
ggplot(gc_res, aes(UncorrectedGenotype,Frequency, fill = ChangedTo )) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = mycols) + xlab("Genotype groups") + ylab("Genotype frequencies in corrected/imputed genotype group")  + 
  guides(fill=guide_legend(title="Genotype after correction"))
write.table(gc_res, file = "genotype_correction_frequencies.txt")

############################################
# Import genotypes and phenotypes as cross #
############################################

#Cross, filtered for SNPs polymorphic between parents
#Seg dist 0.01
#less than 40% missing
mycross <- read.cross("csvsr", "/ddrad_final/",
                      "Pop3_ddrad_bwa_filt_parentalgeno_noUn.rqtl.gen.csv",
                      "Pop3_ddrad_bwa_filt_parentalgeno_noUn.rqtl.phe.csv",
                      genotypes=c("A","H","B"), 
                      alleles=c("A", "B"),
                      estimate.map=FALSE,
                      na.strings=c("-"),
                      bcsft,
                      F.gen = 2,
                      BC.gen = 0)
## Genotype Corrector
# Unplaced_v81 markers dropped
# G-C qchet removed ~50% of markers
# Corrected genotypes with G-C before linkage mapping (window 15, error rate 0.03 for both)
corcross <-read.cross("csvsr", "/ddrad_final/",
                      "Pop3_ddrad_bwa_filt_parentalgeno_noUn.qchet.win15.cor.csv.rqtl.gen.csv",
                      "Pop3_ddrad_bwa_filt_parentalgeno_noUn.qchet.win15.cor.csv.rqtl.phe.csv",
                      genotypes=c("A","H","B"), 
                      alleles=c("A", "B"),
                      estimate.map=FALSE,
                      na.strings=c("-"),
                      bcsft,
                      F.gen = 2,
                      BC.gen = 0)

#########################################################
# Summary statistics on raw and corrected cross objects #
#########################################################

summary(mycross)
#Generate plot to show allele distribution per sample
g <- pull.geno(mycross)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3) 
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))
#Triangle plot
par(mar=rep(8,8), pty="s")
triplot(labels=c("AA","AB","BB"))
tripoints(gfreq, cex=0.8)
tripoints(c(0.25, 0.5, 0.25), col="red", lwd=2, cex=1, pch=4)


############################################
# Run MSTMap algorithm on uncorrected cross#
############################################
pg <- profileGen(corcross, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id =
                   "ID", xo.lambda = 25, layout = c(1, 3), lty = 2,cex = 0.3)
#xo <- pg$stat$xo
#xodf <- as.data.frame(xo)
# Individuals with over 2000 recombinations
#bad_ind <- paste("-",rownames(subset(xodf, xo>2000)),sep="")
# Exclude bad inviduals
#mycross<- subset(mycross, ind=bad_ind)
# Do the same thing with samples displaying high missingness
#miss <- pg$stat$miss
#missdf <- as.data.frame(miss)
#bad_ind_miss <- paste("-",rownames(subset(missdf, miss>6000)),sep="")
#mycross<- subset(mycross, ind=bad_ind_miss)

###################################
#Run MSTMap algorithm on G-C cross#
###################################
options(max.print=100000)

mymapcor0 <- mstmap(corcross,id = "ID", bychr = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                    p.value = 1e-12, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                    miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
save(mymapcor0, file="mymap_cor_p12_chr.RData") 
mymapcor1 <- mstmap(corcross,id = "ID", bychr = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                   p.value = 1e-14, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                   miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
save(mymapcor1, file="mymap_cor_p14.RData") 
mymapcor2 <- mstmap(corcross,id = "ID", bychr = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                   p.value = 1e-16, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                   miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
save(mymapcor2, file="mymap_cor_p16.RData") 
mymapcor3 <- mstmap(corcross,id = "ID", bychr = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                   p.value = 1e-18, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                   miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
save(mymapcor3, file="mymap_cor_p18.RData") 
mymapcor4 <- mstmap(corcross,id = "ID", bychr = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                    p.value = 1e-15, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                    miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
save(mymapcor4, file="mymap_cor_p15.RData") 
mymapcor5 <- mstmap(corcross,id = "ID", bychr = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                    p.value = 1e-13, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                    miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
save(mymapcor5, file="mymap_cor_p13.RData") 
mymapcor6 <- mstmap(corcross,id = "ID", bychr = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                    p.value = 1e-19, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                    miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
save(mymapcor6, file="mymap_cor_p19.RData") 
mymapcor7 <- mstmap(corcross,id = "ID", bychr = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                    p.value = 1e-20, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                    miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
save(mymapcor7, file="mymap_cor_p20.RData") 
mymapcor8 <- mstmap(corcross,id = "ID", bychr = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                    p.value = 1e-21, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                    miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
save(mymapcor8, file="mymap_cor_p21.RData") 
mymapcor9 <- mstmap(corcross,id = "ID", bychr = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                    p.value = 1e-22, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                    miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
save(mymapcor9, file="mymap_cor_p22.RData") 
mymapcor10 <- mstmap(corcross,id = "ID", bychr = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                    p.value = 1e-23, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                    miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
save(mymapcor10, file="mymap_cor_p23.RData") 
load("mymap_cor_p23.RData")

summary(mymapcor10)



########################################################
# Generate final chromosomal LG for map with p = 1e-23 #
########################################################
# Manually check markers on each LG
#markernames(mymapcor7, c("L.1"))
# Unlinked markers on LG with <6 markers are listed
unlinked_markers_cor10 <- c("A01_28857283","A01_28937663","A02_28958077","A10_19543962","C03_147242","C06_44555369",
                          "C05_44878834","C05_44860920","C05_44861671","C05_44875352","C05_44877855","C05_44905493",
                         "C05_44912325","C05_44987014","C05_45042083","C05_44857924","C05_44854138")
# Drop unlinked markers on LG with <6 markers
# Exception: L.2 contains 12 unlinked C01 and A01 markers and they were all dropped as they could not be assigned
mymapcor10_chr <- drop.markers(mymapcor10, unlinked_markers_cor10)

#L22 is C05 satellite but was dropped as unlinked
markernames(mymapcor10_chr, c("A10"))
# Rename LG to chromosomes based on physical marker positions
nam <- names(mymapcor10_chr$geno)
nam[nam=="L.1"] <- "A01"
nam[nam=="L.4"] <- "A02"
nam[nam=="L.7"] <- "A03"
nam[nam=="L.8"] <- "A04"
nam[nam=="L.10"] <- "A05"
nam[nam=="L.11"] <- "A06"
nam[nam=="L.13"] <- "A07"
nam[nam=="L.14"] <- "A08"
nam[nam=="L.15"] <- "A09"
nam[nam=="L.16"] <- "A10"
nam[nam=="L.18"] <- "C01"
nam[nam=="L.5"] <- "C02"
nam[nam=="L.20"] <- "C03"
nam[nam=="L.21"] <- "C05"
nam[nam=="L.9"] <- "C04"
nam[nam=="L.21"] <- "C05"
nam[nam=="L.23"] <- "C06"
nam[nam=="L.25"] <- "C07"
nam[nam=="L.26"] <- "C08"
nam[nam=="L.12"] <- "C09" ##contains some A06 (small physical region, so is prob real)
# Rename
names(mymapcor10_chr$geno)<- nam

mymapcor10_chr <- mstmap(mymapcor10_chr,id = "ID", bychr = TRUE, dist.fun = "kosambi", objective.fun = "COUNT",
                     p.value = 2, noMap.dist = 15, noMap.size = 0, trace = TRUE,
                     miss.thresh = 1, mvest.bc = FALSE, detectBadData = TRUE)
summary(mymapcor10_chr)
plotMap(mymapcor10_chr)

# Heatmaps
# Heatmap of recombination fraction for chromosomes A1 and A2      
#This takes a long time for all markers
chroms <- names(mymapcor10_chr$geno)
chroms <- "A01"
for (i in chroms) {
  outpdf <- paste(i,"_heatmap_lmax50_mymapcor_p23_chr.pdf",sep = "")
  pdf(outpdf) 
  heatMap(mymapcor10_chr, chr = c(i),lmax = 50)
  dev.off()
}


png("myhmap1_plots_300dpi_cex01.png",width = 4, height = 4, units = 'in', res = 300)
myhmap1 <- heatMap(mymapcor7_chr_cor, chr = c("A01" ,"A02" ,"A03", "A04","A05", "C01", "C02", "C03", "C04"),lmax = 50,cex = 0.1)
dev.off()
png("myhmap2_plots.png")
myhmap2 <- heatMap(mymapcor7_chr, chr = c("A06" ,"A07", "A08","A09","A10", "C05", "C06", "C07", "C08", "C09"),lmax = 50)
dev.off()
# decreasing point size will decrease label size
# Vertical/horizontal line artifacts can occur
# These can be avoided by trial and error changes to width and height
# Even small changes can prevent these artifacts
png("hmap_adjusted_h_w.png",width = 1850, height = 1400, pointsize = 5, res = 300,type="cairo")
plotRF(finalmap, chr = c("A01","A08","C01"), zmax = 50)
dev.off()

# Visualise recombination fraction
# This can take a long time and requires large memory
#rf <- pull.rf(mymap)
#lod <- pull.rf(mapthis, what="lod")
#plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

# Save cross object for later import
save(mymapcor10_chr, file="mymap_cor10_p23_chr.RData") 

# This file can be loaded using:
load("mymap_cor10_p23_chr.RData") 
# The mymapcor7_chr object will then be available in the session
summary(mymapcor10_chr)
############################################################################################################
# Manually correct marker orders that show inconsistencies in recombination fraction based on heatmap plots#
############################################################################################################
# A10 has some clear problem markers
# Minor inconsistencies exist on: A1, A2, A3, A4, A6,A9?,C2 C04,C05

### Correct A01
# Bad markers identified manually via an exported  rf table
A01_rf <- pull.rf(mymapcor10_chr,what = c("rf","lod"), chr = "A01")
write.table(A01_rf, file = "mymapcor10_chr_A01_rf.txt")
# Table was opened in excel and all values <0.1 RF marked, allowing identification of bad markers seen in the heatmap
bad_rf_markers <- c("A01_25686249"
                    ,"A01_25686832"
                    ,"A01_25687058"
                    ,"A01_25687943"
                    ,"A01_25691703"
                    ,"A01_25741778"
                    ,"A01_25747088"
                    ,"A01_25767080"
                    ,"A01_26750556"
                    ,"A01_26776582"
)
# Drop bad rf markers
mymapcor10_chr <- drop.markers(mymapcor10_chr, bad_rf_markers)

### Correct A02
A02_rf <- pull.rf(mymapcor10_chr,what = c("rf","lod"), chr = "A02")
write.table(A02_rf, file = "mymapcor10_chr_A02_rf.txt")
# Table was opened in excel and all values <0.1 RF marked, allowing identification of bad markers seen in the heatmap
bad_rf_markers <- c("A02_2460046"
                    ,"A02_2542529"
                    ,"A02_2549079"
                    ,"A02_2565542"
                    ,"A02_2566071"
                    ,"A02_2582608"
                    ,"A02_2589355"
                    ,"A02_3052977"
                    ,"A02_3060196"
                    ,"A02_3074458"
                    ,"A02_3074957"
)

mymapcor10_chr <- drop.markers(mymapcor10_chr, bad_rf_markers)
### Correct A04
A04_rf <- pull.rf(mymapcor10_chr,what = c("rf","lod"), chr = "A04")
write.table(A04_rf, file = "mymapcor10_chr_A04_rf.txt")
# Table was opened in excel and all values <0.1 RF marked, allowing identification of bad markers seen in the heatmap
bad_rf_markers <- c("A04_18725717",
                    "A04_18732733",
                    "A04_19007937"
)
mymapcor10_chr <- drop.markers(mymapcor10_chr, bad_rf_markers)

### Correct A06
A06_rf <- pull.rf(mymapcor10_chr,what = c("rf","lod"), chr = "A06")
write.table(A06_rf, file = "mymapcor10_chr_A06_rf.txt")
# Table was opened in excel and all values <0.1 RF marked, allowing identification of bad markers seen in the heatmap
bad_rf_markers <- c("A06_26029231"
                    ,"A06_26029232"
                    ,"A06_26029219"
                    ,"A06_26029203"
                    ,"A06_26026503"
                    ,"A06_25807086"
                    ,"A06_25839339"
                    ,"A06_25861994"
                    ,"A06_25916355"
                    ,"A06_25917920"
                    ,"A06_25584967"
                    ,"A06_23210887"
                    ,"A06_23208194"
                    ,"A06_23187037"
                    ,"A06_23204979"
                    ,"A06_22831585"
                    ,"A06_22656315"
                    ,"A06_22642300"
                    ,"A06_22655022"
                    ,"A06_22609727"
                    ,"A06_22609968"
                    ,"A06_22386550"
                    ,"A06_22530341"
)

mymapcor10_chr <- drop.markers(mymapcor10_chr, bad_rf_markers)


### Correct A09
A09_rf <- pull.rf(mymapcor10_chr,what = c("rf","lod"), chr = "A09")
write.table(A09_rf, file = "mymapcor10_chr_A09_rf.txt")
# Table was opened in excel and all values <0.1 RF marked, allowing identification of bad markers seen in the heatmap
bad_rf_markers <- c("A09_9628562"
                    ,"A09_9628673"
                    ,"A09_9628712"
                    ,"A09_9628717"
                    ,"A09_9628718"
                    ,"A09_9628721"
                    ,"A09_9628734"
                    ,"A09_9628738"
                    ,"A09_9628745"
                    ,"A09_39780313"
                    ,"A09_39780322"
                    ,"A09_39780325"
                    ,"A09_39780328"
                    ,"A09_39780333"
                    ,"A09_39780337"
                    ,"A09_39780338"
                    ,"A09_39780358"
                    ,"A09_39780329"
)

mymapcor10_chr <- drop.markers(mymapcor10_chr, bad_rf_markers)


### Correct A10
A10_rf <- pull.rf(mymapcor10_chr,what = c("rf","lod"), chr = "A10")
write.table(A10_rf, file = "mymapcor10_chr_A10_rf.txt")
# Table was opened in excel and all values <0.1 RF marked, allowing identification of bad markers seen in the heatmap
bad_rf_markers <- c("A10_10151144"
)

mymapcor10_chr <- drop.markers(mymapcor10_chr, bad_rf_markers)

### Correct C02
# First bin of markers is in reverse order C02_60028 to C02_1364424
C02_rf <- pull.rf(mymapcor10_chr,what = c("rf","lod"), chr = "C02")
write.table(C02_rf, file = "mymapcor10_chr_C02_rf.txt")
# Table was opened in excel and all values <0.1 RF marked, allowing identification of bad markers seen in the heatmap
bad_rf_markers <- c("C02_45655"
                    ,"C02_63297"
                    ,"C02_45163672"
                    ,"C02_45163705"
                    ,"C02_45163677"
                    ,"C02_45163679"
                    ,"C02_45163698"
                    ,"C02_45163699"
                    ,"C02_45163731"
                    ,"C02_45163734"
                    ,"C02_45163740"
                    ,"C02_45163707"
)

mymapcor10_chr <- drop.markers(mymapcor10_chr, bad_rf_markers)


### Correct C04
# First bin of markers is in reverse order C02_60028 to C02_1364424
C04_rf <- pull.rf(mymapcor10_chr,what = c("rf","lod"), chr = "C04")
write.table(C04_rf, file = "mymapcor10_chr_C04_rf.txt")
# Table was opened in excel and all values <0.1 RF marked, allowing identification of bad markers seen in the heatmap
bad_rf_markers <- c("C04_44580361"
                    ,"C04_44580009"
                    ,"C04_44580012"
                    ,"C04_44580072"
                    ,"C04_44580165"
                    ,"C04_44252729"
                    ,"C04_44256870"
                    ,"C04_44259508"
                    ,"C04_44264364"
                    ,"C04_44052775"
                    ,"C04_44080746"
                    ,"C04_44104561"
                    ,"C04_44110803"
)

mymapcor10_chr <- drop.markers(mymapcor10_chr, bad_rf_markers)
### Correct C05
C05_rf <- pull.rf(mymapcor10_chr,what = c("rf","lod"), chr = "C05")
write.table(C05_rf, file = "mymapcor10_chr_C05_rf.txt")
# Table was opened in excel and all values <0.1 RF marked, allowing identification of bad markers seen in the heatmap
bad_rf_markers <- c("C05_45061850"
                    ,"C05_45231162"
                    ,"C05_45247384"
                    ,"C05_45284044"
                    ,"C05_45352868"
                    
)
mymapcor10_chr <- drop.markers(mymapcor10_chr, bad_rf_markers)
save(mymapcor10_chr, file="mymap_cor10_p23_chr_cor.RData") 

for (i in chroms) {
  outpdf <- paste(i,"_heatmap_lmax50_mymapcor_p29_chr_cor.pdf",sep = "")
  pdf(outpdf) 
  heatMap(mymapcor7_chr_cor, chr = c(i),lmax = 50)
  dev.off()
}

#########################################
# Flip marker order when LG is reversed #
#########################################

# Reverse order of markers on chr based on physical order
# Data marker order and map order were not consistent anymore on A03

# Manually check whether marker order needs to be flipped based on marker names containing physical positions
mymapcor10_chr$geno$C02$map

# Flip marker order on identified chromosomes
mymapcor10_chr_cor  <- flip.order(mymapcor10_chr, c("A03","A05","A06","A10","C04","C06","C07","C09"))

# Plot RF and LOD heatmaps per chr using corrected orientation
chroms <- names(mymapcor10_chr_cor$geno)

for (i in chroms) {
  outpdf <- paste(i,"_heatmap_lmax50_mymapcor_p29_chr_cor_oriented.pdf",sep = "")
  pdf(outpdf) 
  heatMap(mymapcor10_chr_cor, chr = c(i),lmax = 50)
  dev.off()
}

#####################################
# Visualise genetic maps            #
#####################################

#Change marker names to remove redundant "chr" in front of chromosome
#write.cross(finalmap, filestem = "tmp", format = "csvsr")
finalmap <-read.cross("csvsr", ".",
                      "tmp_gen.csv",
                      "tmp_phe.csv",
                     genotypes=c("AA","AB","BB"), 
                      alleles=c("A", "B"),
                      estimate.map=FALSE,
                      na.strings=c("-"),
                      bcsft,
                      F.gen = 2,
                      BC.gen = 0)
# Marker names on different chromosomes can't be identical
# So one marker was renamed to be unique: A10_16607778 A03_16607778
lmv.linkage.plot(finalmap,"finalmap_plain.pdf")
lmv.linkage.plot(finalmap,"finalmap_plain_A.pdf", mapthese = c("A01","A02","A03","A04","A05","A06","A07","A08", "A09","A10"))
lmv.linkage.plot(finalmap,"finalmap_plain_C.pdf",mapthese = c("C01","C02","C03","C04","C05","C06","C07","C08", "C09"))

# Density plot
## draw tickmarks at each cM from 0 to largest position of linkage groups to be drawn
#maxpos <- floor(max(finalmap$Position[finalmap$Group == "A01" | finalmap$Group == "A02"]))
maxpos <- 250
at.axis <- seq(0, maxpos)

## put labels on ruler at every 10 cM
axlab <- vector()
for (lab in 0:maxpos) {
  if (!lab %% 10) {
    axlab <- c(axlab, lab)
  }
  else {
    axlab <- c(axlab, NA)
  }
}

lmv.linkage.plot(finalmap,"finalmap_density.pdf",denmap=TRUE, cex.axis = 1, at.axis = at.axis, labels.axis = axlab)


# 6 Non-unique positions were found
#A01/10:1361665, C01/A10:15032425, A03/A10:16607778, A08/A03:17161918, A03/A09:17519533, A07/A03:20543703
# Density plot
## draw tickmarks at each cM from 0 to largest position of linkage groups to be drawn
#maxpos <- floor(max(finalmap$Position[finalmap$Group == "A01" | finalmap$Group == "A02"]))

#maxpos <- 2100
#at.axis <- seq(0, maxpos)

## put labels on ruler at every 10 cM
#axlab <- vector()
#for (lab in 0:maxpos) {
#  if (!lab %% 10) {
#    axlab <- c(axlab, lab)
#  }
#  else {
#    axlab <- c(axlab, NA)
#  }
#}
#lmv.linkage.plot(mycrossmap18_chr,"mycrossmap18_chr_density.pdf",denmap=TRUE, cex.axis = 1, at.axis = at.axis, labels.axis = axlab)

######################################
#Correlation of physical and genetic map#
######################################

# Spearman rank correlation test
# Where x and y are vectors of centimorgan positions and physical positions per chr
#cor.test(A01_order_df_t$genpos, A01_order_df_t$pos, method="spearman")

# Plot correlation and test results for each chromosome
chroms <- names(qtlmap$geno)
for (i in chroms) {
  # Get df of marker names and positions
  chr_order <- pull.map(qtlmap, chr=i)
  chr_order_df <- as.data.frame(as.list(get(i,chr_order)))
  chr_order_df_t <- data.frame(t(chr_order_df[-1]))
  # Rename genetic distance column
  colnames(chr_order_df_t)[1] <- "genpos" 
  chr_order_df_t$names <- rownames(chr_order_df_t)
  # Split the marker name into chr and pos and make pos numeric with convert
  chr_order_df_t <- chr_order_df_t %>%
    separate(names, c("chr", "pos"), "_", convert = TRUE)
  # 119 markers have a conflict between physical and genetic chromosome location
  # Only markers without conflict will be used
  chr_order_df_t <- subset(chr_order_df_t, chr == i)
  chr_order_df_t <- transform(chr_order_df_t, pos_mbp = pos / 1000000)
  maintitle <- paste("Correlation of genetic and physical map of ",i,sep = "")
  outpdf <- paste(i,"_spearman_correlation_phys_gen.pdf",sep = "")
  pdf(outpdf) 
  # carry out spearman correlation test and make scatterplot
  # ggplot call has to be in print statement when in loop
  print(ggscatter(chr_order_df_t, x = "pos_mbp", y = "genpos", 
            conf.int = TRUE, main = maintitle,
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "Physical pos (Mbp)", ylab = "Genetic pos (cM)"))
  dev.off()
}

#############################################
# Obtain statistics for final mapping data  #
#############################################

# Function to convert list of matrices to list of data frames
matrixList2dataframeList <- function(x)
{
  if(any(class(x)=="list"))
  {
    x_copy = x
    attrs = setdiff(names(attributes(x)),"names")
    x = lapply(x,replace_sub_dataframes)
    if(length(attrs)>0)
    {
      for(i in 1:length(attrs))
      {
        attr(x,attrs[i]) <- matrixList2dataframeList (attr(x_copy,attrs[i]))
      }
    }
    return(x)
  }
  else
  {
    if(any(class(x)=="matrix")) 
      return(as.data.frame(x)) 
    else 
      return(x)
  }
}



# Plot genotypes by individual and by marker to show data completeness
pdf("finalmap_ntyped_plots.pdf")
par(mfrow=c(1,2), las=1)
plot(ntyped(qtlmap), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(qtlmap, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
dev.off()

# Check for clonal samples that need to be merged
# After genotype correction and imputation, the similarity between some samples increases substantially
clonecheck <- genClones(qtlmap, tol = 0.9, id = "ID")
sink('qtlmap_clonecheck.txt')
clonecheck$cgd
sink()

png("A03_geno_before_after_correction.png", width = 1834, height = 1395, pointsize = 5, res = 300,type="cairo")
par(mfrow=c(1,2))
geno.image(mycross, chr = "A03", main = "Genotypes before correction (AA: red, AB: blue, BB:green)")
geno.image(corcross, chr = "A03", main = "Genotypes after correction (AA: red, AB: blue, BB:green)")
dev.off()

#png("genetic_map_before_after_correction.png", width = 1850, height = 1400, pointsize = 5, res = 300,type="cairo")
#par(mfrow=c(1,2))
#plotMap(mycross)
#plotMap(qtlmap)
#dev.off()

pdf("A03_genotypes_before_after_correction_mycross_corcross.pdf", pointsize = 2)
par(mfrow=c(1,2))
geno.image(mycross, chr = "A03",main = "Genotypes before correction (AA: red, AB: blue, BB:green)")
geno.image(corcross, chr = "A03", main = "Genotypes after correction (AA: red, AB: blue, BB:green)")
dev.off()

pdf("genotypes_before_after_correction_mycross_corcross.pdf", pointsize = 2)
par(mfrow=c(1,2))
geno.image(mycross, main = "Genotypes before correction (AA: red, AB: blue, BB:green)")
geno.image(corcross, main = "Genotypes after correction (AA: red, AB: blue, BB:green)")
dev.off()

# Check genotypes for a single sample to evaluate locateXO accuracy
qtlmap_tmp<- subset(qtlmap, ind=c("129", "150", "107","156"))
pdf("tmp.pdf")
geno.image(qtlmap_tmp,chr = "C09")
dev.off()

# Estimate the locations of crossovers (CO) for each individual on a given chromosome. 
# Loop through chromsomes and get positions and frequency of CO
datalist = list()
maplist = list()
chroms <- names(qtlmap$geno)
for (i in chroms) {
  outfile <- paste(i,"_crossover_per_ind.csv",sep = "")
  myxo <- locateXO(qtlmap, chr = i, full.info=TRUE)
  myxodf = lapply(myxo, matrixList2dataframeList)
  # Remove individuals with no CO
  myxodfc <- myxodf[sapply(myxodf, length) > 0]
  # Concatenate df and insert ID column with individual name
  myxo_all <- bind_rows(myxodfc, .id = "IND")
  # Add chromosome column
  #myxo_all['chr']=i
  # Add to list of chr df
  datalist[[i]] <- myxo_all
  # Get marker positions on chr
  chrgeno <- get(i,qtlmap$geno)
  chrmap <- chrgeno$map
  mapdf <- enframe(chrmap)
  maplist[[i]] <- mapdf
  # Write dataframe per chr for bash/python processing
  #write.csv(myxo_all,file = outfile, row.names=FALSE,quote=FALSE)
}
all_chr <- bind_rows(datalist, .id = "chr")
all_mark <- bind_rows(maplist, .id = "chr")
write.csv(all_chr,file = "Pop3_ddrad_qtlmap_locatexo_all_chr.csv", row.names=FALSE,quote=FALSE)
# Markers flanking a crossover can be assigned physical positions using this table
write.csv(all_mark,file = "Pop3_ddrad_qtlmap_marker_positions.csv", row.names=FALSE,quote=FALSE)


##################################################
# QTL analysis scanone for all traits            #
##################################################

#Loop through traits (except hairiness)
#Output LOD per trait as table
#Output permutation results as table
#Importantly this analysis assumes normal distribution of traits (not true for most traits)


write.cross(mymapcor10_chr_cor, filestem = "finalmap_p23_cor_hetcleanGC_dropbadrfmarkers", format = "csvsr")
# The finalmap cross is encoded as bcsft 2 with 0 backcross generation
# This was recommended for F2 in a manual
# However, the object can also be encoded as a plain F2, which seems neater and may impact some functions
qtlmap <-read.cross("csvsr", ".",
                 "finalmap_p23_cor_hetcleanGC_dropbadrfmarkers_gen.csv",
                 "finalmap_p23_cor_hetcleanGC_dropbadrfmarkers_phe.csv",
                 genotypes=c("AA","AB","BB"), 
                 alleles=c("A", "B"),
                 estimate.map=FALSE,
                 na.strings=c("-"))


qtlmap$pheno
traits <- c(2,3,4,6,7)
for (i in traits) {
  # Scan single QTL
  result_hk <- scanone(qtlmap, pheno.col=i, method= "hk")     #Scan the first column using Haley-Knott
  result_np <- scanone(qtlmap, pheno.col=i, model = "np")
  out.cim <- cim(qtlmap, pheno.col=i,n.marcovar=3) #Composite interval mapping

  # Use permutations to calculate thresholds
  result_cim <- cim(qtlmap, pheno.col=i, n.marcovar=3, n.perm = 1000)
  perm_hk <- scanone(qtlmap, pheno.col=i, method= "hk",n.perm= 1000)
  perm_np <- scanone(qtlmap, pheno.col=i, model = "np",n.perm= 1000)

  # Get suggestive and significant threshold to two decimal places
  lowerbound_cim <- format(round(summary(result_cim)[[2]], 2), nsmall = 2)
  upperbound_cim <- format(round(summary(result_cim)[[1]], 2), nsmall = 2)
  lowerbound_hk <- format(round(summary(perm_hk)[[2]], 2), nsmall = 2)
  upperbound_hk <- format(round(summary(perm_hk)[[1]], 2), nsmall = 2)
  lowerbound_np <- format(round(summary(perm_np)[[2]], 2), nsmall = 2)
  upperbound_np <- format(round(summary(perm_np)[[1]], 2), nsmall = 2)
  
  maintitle <- paste("Pop3 QTL Scanone Trait ",i,sep = "")
  outhk <- paste(i,"_scanone_1000perm","_hk",lowerbound_hk,"-",upperbound_hk,
                   sep = "")
  outnp <- paste(i,"_scanone_1000perm","_np",lowerbound_np,"-",upperbound_np,
                 sep = "")
  outnppdf <- paste(outnp,".pdf", sep = "")
  outpdf <- paste(outhk,".pdf", sep = "")
  out_lod <- paste(outhk,"_lod.csv", sep = "")
  total <- merge(perm_hk,perm_np, by="row.names",all.x=TRUE)
  total <- merge(total, result_cim,by="row.names",all.x=TRUE)
  colnames(total) <- c("r1", "r2", "lod_hk","lod_np", "lod_cim")
  total <- select(total, lod_hk, lod_np, lod_cim)
  write.csv(total, file = out_lod)
 
  pdf(outpdf) 
  plot(result_hk,lwd=c(4,3,1),col=c("black"),main = maintitle, cex = 0.5)
  abline(h=summary(perm_hk)[[1]],col="black",lty=2,lwd=2)
  abline(h=summary(perm_hk)[[2]],col="gray",lty=2,lwd=2)
  dev.off()
  
  pdf(outnppdf) 
  # carry out spearman correlation test and make scatterplot
  # ggplot call has to be in print statement when in loop
  plot(result_np,lwd=c(4,3,1),col=c("black"),main = maintitle, cex = 0.5)
  abline(h=summary(perm_np)[[1]],col="black",lty=2,lwd=2)
  abline(h=summary(perm_np)[[2]],col="gray",lty=2,lwd=2)
  dev.off()
  
  
  cimtitle <- paste("Pop3 QTL CIM vs HK Trait ",i,sep = "")
  cim_outpdf <- paste(i,"_1000perm_CIM","_hk",lowerbound_hk,"-",upperbound_hk,
                  "_cim", lowerbound_cim,"-",upperbound_cim, 
                  ".pdf", sep = "")
  pdf(cim_outpdf)
  # Plot CIM and scanone results to compare them
  plot(result_hk,out.cim, col =c("red","blue"), main = cimtitle, cex = 0.5)
  legend("topleft",c("hk","cim"),lwd=c(4,3,1),col=c("red","blue"))
  abline(h=summary(perm_hk)[[1]],col="red",lty=2,lwd=2)
  abline(h=summary(result_cim)[[1]],col="blue",lty=2,lwd=2)
  dev.off()
}

result_np <- scanone(qtlmap_noclone, pheno.col=7, model = "np")
pdf("tmp.pdf")
plot(result_np)
dev.off()
# Remove pseudomarkers from scanone result
# https://groups.google.com/forum/#!topic/rqtl-disc/YjXZY_xlEc8
result_np <- result_np[markernames(qtlmap),] 

#Result Bud time QTL LOD interval
#lodint(result_np, drop = 1.5, chr="C02")

#Result Flowering time QTL LOD interval
#lodint(result_np, drop = 1.5, chr="C02")

# Summary of the significant loci
summary(result_np, perms=operm.np, alpha=0.05, pvalues=TRUE)



#########################
# QTL analysis hairiness#
#########################

# Hairiness is a binary trait so the binary model is used

result_hk <- scanone(qtlmap, pheno.col=5, method= "hk", model = "binary")     #Scan the first column using Haley-Knott
perm_hk <- scanone(qtlmap, pheno.col=5, method= "hk",n.perm=1000,model = "binary")

lowerbound_hk <- format(round(summary(perm_hk)[[2]], 2), nsmall = 2)
upperbound_hk <- format(round(summary(perm_hk)[[1]], 2), nsmall = 2)

# Plot
outpdf <- paste("Hairiness_perm1000_scanone","_hk",lowerbound_hk,"-",upperbound_hk,
                ".pdf", sep = "")
pdf(outpdf) 
plot(result_hk,lwd=c(4,3),col=c("black"),main = "Pop3 QTL Scan Hairiness (binary)", cex = 0.5)
abline(h=summary(perm_hk)[[1]],col="black",lty=2,lwd=2)
abline(h=summary(perm_hk)[[2]],col="gray",lty=2,lwd=2)
dev.off()
out_lod <- paste("Hairiness_perm1000_scanone","_hk",lowerbound_hk,"-",upperbound_hk,
                 "_lod.txt", sep = "")
write.csv(perm_hk, file = out_lod)

########################
# Marker effect sizes  #
########################
# Plot Phenotype versus genotype at significant marker flowering
pdf("qtlmap_plotpxg_flower.pdf")
plotPXG(qtlmap,"C02_4655461", pheno.col = 7)
dev.off()
# Plot Phenotype versus genotype at significant marker bud
pdf("qtlmap_plotpxg_bud.pdf")
plotPXG(qtlmap,"chr_5604", pheno.col = 6)
dev.off()

# The function fitqtl was used for calculating percentages of variance of the significant QTL 
# by calculating the coefficient of determination for each single-QTL model obtained using scanone
sug <- calc.genoprob(qtlmap, step=1)
# Flowering time QTL effect
sink('qtl_effect_flowering_qtlmap_step1.txt')
qtl <- makeqtl(sug, chr="C02", pos=1, what="prob")
out.fq <- fitqtl(sug, qtl=qtl, method="hk", pheno.col=7)
summary(fitqtl(sug,  pheno.col=7, qtl=qtl, method="hk", get.ests=TRUE, dropone=FALSE))
sink()
# Bud QTL effect
sink('qtl_effect_bud_qtlmap_step1_C02.txt')
qtl <- makeqtl(sug, chr="C02", pos=1, what="prob")
out.fq <- fitqtl(sug, qtl=qtl, method="hk", pheno.col=6)
summary(fitqtl(sug,  pheno.col=6, qtl=qtl, method="hk", get.ests=TRUE, dropone=FALSE))
sink()

################################
# Visualising the QTL on C02   #
################################
# Highlight single locus on linkage group
outfile = file.path(tempdir(), "hyper_showonly.pdf")
lmv.linkage.plot(hyper,outfile,mapthese=c(1,4,6,15),lcol="green",lcex=2,lfont=2,
                 rcol="red",rcex=2,rfont=3,
                 showonly=c("D1Mit123","D4Mit80","D6Mit135","D15Mit156"))

# Show QTL
qtldf <- data.frame(
  chr = character(),
  qtl = character(),
  so = numeric(),
  si = numeric(),
  ei = numeric(),
  eo = numeric(),
  col = character(),
  stringsAsFactors = FALSE
)

qtldf <- rbind(qtldf,
               data.frame(
                 chr = c("C02","C02","C02","C02","C02","C02","C02"),
                 qtl = c("Gene1","Gene2", "Gene3","Gene4","Gene5","Gene6","Gene7"),
                 so = c(98.1810039,124.5663771,123.4299626,123.1458656,122.2935537,122.0094547,94.08807),
                 si = c(98.1810039,124.5663771,123.4299626,123.1458656,122.2935537,122.0094547,94.08807),
                 ei = c(99.9430977,125.1345844,123.9981698,123.1458666,122.2935547,122.0094557,125.98694),
                 eo = c(99.9430977,125.1345844,123.9981698,123.1458666,122.2935547,122.0094557,125.98694),
                 col=c("green","red","red","red","red","red","black")
               ))
# make a list to pass label options
flist <- list()
locus <- c("chrC02_4324", "chrC02_45329")
font  <- c(2)   #bold
flist[[1]] <- list(locus = locus, font = font)
locus <- c("F3H", "FLS1")
font  <- c(4)   #bold italic
flist[[2]] <- list(locus = locus, font = font)
locus <- c()
font  <- c(3)   #italic
col <- c("red")
flist[[3]] <- list(locus = locus, font = font, col = col)

# use showonly parameter to only show specific loci.
outfile = "finalmap_f2_QTL_candidate_gene_on_C02.pdf"
lmv.linkage.plot(
  finalmap_f2,
  mapthese=c("C02"),
  showonly = c("chrC02_214513", "chrC02_315123"),
  outfile = outfile,
  ruler = TRUE,
  lgtitle = c("C02"),
  maxnbrcolsfordups = 2,
  markerformatlist = flist,
  lg.col = "lightblue1",
  pdf.height = 16,
  pdf.width = 10,
  qtldf=qtldf
)







