library(NOISeq)
library(tidyr)
library(ggplot2)
library(SuperExactTest)
library(limma)
library(topGO)
library(gplots)
library(data.table)

## Part 1. DE in each of the three species

setwd("Tripedalia/")

counts <- read.table("counts_matrix.txt", header=T)

lmat <- data.frame(id=counts$Geneid, l=counts[,6])
dmat <- data.frame(counts[,7:10])
rownames(dmat) <- counts$Geneid

trip_dat <- readData(data=dmat, length=lmat, 
                     factors=data.frame(Tissue=c("Man", "Rhop", "Tent", "Body")))

set.seed('123')
myresults <- noiseq(trip_dat, factor = "Tissue", condition=c("Rhop", "Man"), lc = 1, replicates = "no")

resRT <- noiseq(trip_dat, factor = "Tissue", condition=c("Rhop", "Tent"), lc = 1, replicates = "no")
resTM <- noiseq(trip_dat, factor = "Tissue", condition=c("Man", "Tent"), lc = 1, replicates = "no")

resRM <- myresults@results[[1]]
resRT <- resRT@results[[1]]
resTM <- resTM@results[[1]]

write.csv(resRT, file="Trip_Noiseq_RT.csv", quote=F)
write.csv(resTM, file="Trip_Noiseq_TM.csv", quote=F)
write.csv(myresults@results[[1]], file="Trip_Noiseq_res.csv", quote=F)

trip_counts <- counts
x <- trip_counts[,7:10] / trip_counts$Length
trip_TPM <- t( t(x) * 1e6 / colSums(x) )
trip_TPM10k <- trip_TPM * nrow(trip_counts) / 10000
rownames(trip_TPM10k) <- trip_counts$Geneid

write.csv(trip_TPM10k, file="Trip_TPM10K.csv", quote=F)

##Aurelia DESeq
setwd('../Aurelia/')

Aurelia_counts <-read.table("counts_matrix.txt", header=T)
aurelia_samples <- read.table('samples.tsv', header=T)

counts <- Aurelia_counts[,7:15]
names(counts) <- c("MJ1", "MJ23", "MJ26", "MJ27", "MJ29", "MJ32", "MJ35",
                   "MJ4", "MJ8") #get into the order of metadata

counts <- counts[,c(1,8,9,2:7)] 
names(counts)==aurelia_samples$sample_name
rownames(counts) <- Aurelia_counts$Geneid

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=aurelia_samples, 
                              design=~condition)

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
summary(res) 
head(res)

write.csv(res, file="Aurelia_DESeq_RM.csv", quote=F)

x <- counts / Aurelia_counts$Length
Aurelia_TPM <- t( t(x) * 1e6 / colSums(x) )
Aurelia_TPM10k <- Aurelia_TPM * nrow(Aurelia_counts) / 10000

write.csv(Aurelia_TPM10k, file="Aurelia_TPM10K.csv", quote=F)

##Sarsia DESeq
setwd("../Sarsia/")

Sarsia_counts <-read.table("counts_matrix.txt", header=T)
Sarsia_samples <- read.table('../Sarsia/samples.tsv', header=T)

counts <- Sarsia_counts[,7:12]
names(counts) <- c("SA01", "SA02", "SA03",
                   "SA04", "SA05", "SA06") #get into the order of metadata

names(counts)==Sarsia_samples$sample_name
rownames(counts) <- Sarsia_counts$Geneid

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=Sarsia_samples, 
                              design=~condition)

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, contrast = c("condition", "tentacle_bulbs", "manubrium"))
summary(res) 
head(res)

write.csv(res, file="Sarsia_DESeq_TBM.csv", quote=F)

res_TBT <- results(dds, alpha = 0.05, contrast = c("condition", "tentacle_bulbs", "tentacles"))
write.csv(res, file="Sarsia_DESeq_TBT.csv", quote=F)

res_TM <- results(dds, alpha = 0.05, contrast = c("condition", "tentacles", "manubrium"))
write.csv(res, file="Sarsia_DESeq_TM.csv", quote=F)

x <- counts / Sarsia_counts$Length
Sarsia_TPM <- t( t(x) * 1e6 / colSums(x) )
Sarsia_TPM10k <- Sarsia_TPM * nrow(Sarsia_counts) / 10000

write.csv(Sarsia_TPM10k, file="Sarsia_TPM10K.csv", quote=F)

##Part 2. Gene expression visualizations
tissue_cols <- c("#94d10e", "#40b693", "#e47cbc", "#fb7d40", "#768ec2")
species_cols <- c("#08A54F", "#0A75BC", "#D5003E")

setwd("../Aurelia/")
Aurelia_TPM10k <- read.csv("Aurelia_TPM10K.csv", row.names = 1)

pr <- prcomp(t(log2(Aurelia_TPM10k+0.001)))

to_plot <- data.frame(pr$x)
to_plot$Tissue <- aurelia_samples[match(rownames(to_plot), aurelia_samples$sample_name),]$condition

summary(pr)$importance[2,]

png(res=400, width=6, height=4, units='in', file="Aurelia_PCA.png")
ggplot(to_plot, aes(x=PC1, y=-PC2, color=Tissue)) + geom_point(size=3) + theme_bw() + 
  scale_color_manual(values=tissue_cols[2:3]) + xlab("PC1 (40.5%)") + ylab("PC2 (17.3%)") +
  theme(text=element_text(size=14))
dev.off()

#Trip
setwd("../Tripedalia/")
Trip_TPM10k <- read.csv("Trip_TPM10K.csv", row.names = 1)

pr <- prcomp(t(log2(Trip_TPM10k+0.001)))

to_plot <- data.frame(pr$x)
to_plot$Tissue <- c("Man", "Rhop", "Tent", "Body")

png(res=400, width=6, height=4, units='in', file="Trip_PCA.png")
ggplot(to_plot, aes(x=PC1, y=-PC2, color=Tissue)) + geom_point(size=3) + theme_bw() + 
  scale_color_manual(values=tissue_cols[c(1:3,5)]) + xlab("PC1 (46.5%)") + ylab("PC2 (30.7%)")  +
  theme(text=element_text(size=14))
dev.off()

#Sarsia
setwd("./Sarsia/")
Sarsia_TPM10k <- read.csv("Sarsia_TPM10K.csv", row.names = 1)

pr <- prcomp(t(log2(Sarsia_TPM10k+0.001)))

to_plot <- data.frame(pr$x)
to_plot$Tissue <- Sarsia_samples[match(rownames(to_plot), Sarsia_samples$sample_name),]$condition

png(res=400, width=6, height=4, units='in', file="Sarsia_PCA.png")
ggplot(to_plot, aes(x=PC1, y=PC2, color=Tissue)) + geom_point(size=3) + theme_bw() + 
  scale_color_manual(values=tissue_cols[c(2,4,5)]) + xlab("PC1 (34.6%)") + ylab("PC2 (26.1%)")  +
  theme(text=element_text(size=14))
dev.off()


#Cross-species (Orthogroups)
setwd("..")

OGs <- read.table("N2.tsv", fill=T, header=T, sep="\t")

OGs[OGs==""] = NA

OGs <- OGs %>% dplyr::select(HOG, Aurelia, Sarsia, Tripedalia) %>% 
  dplyr::filter(!(is.na(Aurelia) & is.na(Sarsia) & is.na(Tripedalia)))

All_sc <- read.csv("../All_SC_expression.csv")

for_pca <- All_sc %>% dplyr::select(-Tissue, -Species, -Gene) %>% pivot_wider(names_from = "Library", values_from="Expr")
for_pca[,2:20] <- log2(for_pca[,c(2:20)]+0.001)
pr <- prcomp(t(for_pca[,c(2:20)]))

to_plot <- data.frame(pr$x)
to_plot$Tissue <- c(aurelia_samples$condition, Sarsia_samples$condition, trip_samples$condition)
to_plot$Species <- c(rep("Aurelia", 9), rep("Sarsia", 6), rep("Tripedalia", 4))

summary(pr)$importance[2,]

png(res=400, width=6, height=4, units='in', file="All_PCA.png")
ggplot(to_plot, aes(x=PC1, y=-PC2, color=Tissue, shape=Species)) + geom_point(size=3) + theme_bw() + 
  scale_color_manual(values=tissue_cols) + xlab("PC1 (57.8%)") + ylab("PC2 (14.4%)")  +
  theme(text=element_text(size=14))
dev.off()

quant_mat <- for_pca
quant_mat[,2:20] <- normalizeQuantiles(as.matrix(quant_mat[,c(2:20)]))

pr <- prcomp(t(quant_mat[,c(2:20)]))

to_plot <- data.frame(pr$x)
to_plot$Tissue <- c(aurelia_samples$condition, Sarsia_samples$condition, trip_samples$condition)
to_plot$Species <- c(rep("Aurelia", 9), rep("Sarsia", 6), rep("Tripedalia", 4))

summary(pr)$importance[2,]

png(res=400, width=6, height=4, units='in', file="All_PCA_quant.png")
ggplot(to_plot, aes(x=PC1, y=-PC2, color=Tissue, shape=Species)) + geom_point(size=3) + theme_bw() + 
  scale_color_manual(values=tissue_cols,labels=c("Body", "Manubrium", "Rhopalia", "Tentacle bulbs", "Tentacles")) +
  xlab("PC1 (49.2%)") + ylab("PC2 (18.5%)")  +
  theme(text=element_text(size=16)) 
dev.off()

#perform Z-scoring of genes within each species

Z_aurelia <- quant_mat[,2:10]
Z_sarsia <- quant_mat[,11:16]
Z_trip <- quant_mat[,17:20]

Z_aurelia <- (Z_aurelia-rowMeans(Z_aurelia)) / apply(Z_aurelia,1,sd)
Z_sarsia <- (Z_sarsia-rowMeans(Z_sarsia)) / apply(Z_sarsia,1,sd)
Z_trip <- (Z_trip-rowMeans(Z_trip)) / apply(Z_trip,1,sd)

#some of these have no variability and so show up as na 
Z_all <- cbind(Z_aurelia, Z_sarsia, Z_trip)

Z_all[is.na(Z_all)] <- 0
pr <- prcomp(t(Z_all))

to_plot <- data.frame(pr$x)
to_plot$Tissue <- c(aurelia_samples$condition, Sarsia_samples$condition, trip_samples$condition)
to_plot$Species <- c(rep("Aurelia", 9), rep("Sarsia", 6), rep("Tripedalia", 4))

summary(pr)$importance[2,]

png(res=400, width=6, height=4, units='in', file="All_PCA_Zscore.png")
ggplot(to_plot, aes(x=PC1, y=-PC2, color=Tissue, shape=Species)) + geom_point(size=3) + theme_bw() + 
  scale_color_manual(values=tissue_cols,labels=c("Body", "Manubrium", "Rhopalia", "Tentacle bulbs", "Tentacles")) +
  xlab("PC1 (23.9%)") + ylab("PC2 (11.6%)")  +
  theme(text=element_text(size=14))
dev.off()

#Pearson correlation matrix
cor_mat <- cor(Z_all, method = "pearson")    

png("Pearson_cor_sp.png",height=8,width=8, units='in', res=400)
heatmap.2(as.matrix(cor_mat),keysize=0.7, key.xlab = "Pearson Correlation", key.title = "none",key=T, 
          density.info="none",trace="none",
          col=colorpanel(100,"white","blue"),
          ColSideColors=tissue_cols[c(3,3,2,2,2,3,2,2,2,5,4,2,5,4,2,2,3,5,1)], RowSideColors=species_cols[c(rep(1,9), rep(3,6), rep(2,4))],
          margin=c(12, 12), main="Sample clustering, Pearson correlation, Z-scored")
legend("topright", legend=c("Body","Manubria", "Rhopalia", "Tentacle bulbs", "Tentacles"), col=tissue_cols,
       fill=tissue_cols, border=FALSE, cex=0.7)
dev.off()

##Part 3. GO
# Gene Ontology Analysis

# Read the GO annotation map
geneID2GO <- readMappings(file = "Tripedalia/Trip_0524.go", sep="\t")
geneUniverse<-names(geneID2GO)

## Make a factor object of the DEG genes that are upregulated
genesRM_UP <- subset(res_trip, M>0 & prob>0.99)

genesRM_UP <- as.character(genesRM_UP$X)
geneListRM_UP <- factor(as.integer(geneUniverse%in%genesRM_UP))
names(geneListRM_UP)<-geneUniverse

## Build topGO object
GOdataBP_RM_UP <- new('topGOdata',
                      description = "GO analysis of Trip tissues; Upregulated in Rhopalia versus Manubrium",
                      ontology = 'BP', # Biological Process
                      allGenes = geneListRM_UP,
                      nodeSize=10,
                      annot = annFUN.gene2GO,
                      gene2GO = geneID2GO)

## Run Fisher's exact test 
resultFisherBP_RM_UP <- runTest(GOdataBP_RM_UP, algorithm="weight01", statistic="fisher") 

## Summarize the results for all GO terms
GO_resultsBP_RM_UP <- GenTable(GOdataBP_RM_UP, classic = resultFisherBP_RM_UP,  
                               ranksOf = "classic", topNodes = length(resultFisherBP_RM_UP@score))


subset(GO_resultsBP_RM_UP, as.numeric(classic)<0.01) #80 in Trip

## Export the results to a file
write.csv(GO_resultsBP_RM_UP, file="Tripedalia/GO_resultsBP_RM_UP.csv")


## Sarsia GO
# Read the GO annotation map
geneID2GO <- readMappings(file = "Sarsia/Sarsia_0524.go", sep="\t")
geneUniverse<-names(geneID2GO)

## Make a factor object of the DEG genes that are upregulated
genesRM_UP <- subset(res_sarsia, log2FoldChange>0 & padj<0.05)

genesRM_UP <- as.character(genesRM_UP$X)
geneListRM_UP <- factor(as.integer(geneUniverse%in%genesRM_UP))
names(geneListRM_UP)<-geneUniverse

## Build topGO object
GOdataBP_RM_UP <- new('topGOdata',
                      description = "GO analysis of Trip tissues; Upregulated in Rhopalia versus Manubrium",
                      ontology = 'BP', # Biological Process
                      allGenes = geneListRM_UP,
                      nodeSize=10,
                      annot = annFUN.gene2GO,
                      gene2GO = geneID2GO)

## Run Fisher's exact test 
resultFisherBP_RM_UP <- runTest(GOdataBP_RM_UP, algorithm="weight01", statistic="fisher") 

## Summarize the results for all GO terms
GO_resultsBP_RM_UP <- GenTable(GOdataBP_RM_UP, classic = resultFisherBP_RM_UP,  
                               ranksOf = "classic", topNodes = length(resultFisherBP_RM_UP@score))

subset(GO_resultsBP_RM_UP, as.numeric(classic)<0.01) #208 in Sarsia

## Export the results to a file
write.csv(GO_resultsBP_RM_UP, file="Sarsia/GO_resultsBP_RM_UP.csv")

##Aurelia GO
#Aight let's compare GO terms among species...

Aurelia_GO <- read.csv("Aurelia/GO_resultsBP_RM_UP.csv")
Sarsia_GO <- read.csv("Sarsia/GO_resultsBP_TBM_UP.csv")
Trip_GO <- read.csv("Tripedalia/GO_resultsBP_RM_UP.csv")

Aurelia_sig <- Aurelia_GO[as.numeric(Aurelia_GO$classic)<0.01,] 
Sarsia_sig <- Sarsia_GO[as.numeric(Sarsia_GO$classic)<0.01,]
Trip_sig <- Trip_GO[as.numeric(Trip_GO$classic)<0.01,]

GO_all3<-intersect(intersect(Aurelia_sig$GO.ID,Sarsia_sig$GO.ID),Trip_sig$GO.ID) 
#11

GO_AT <- intersect(Aurelia_sig$GO.ID,Trip_sig$GO.ID)#26
GO_AS <- intersect(Aurelia_sig$GO.ID,Sarsia_sig$GO.ID) #43
GO_ST <- intersect(Sarsia_sig$GO.ID,Trip_sig$GO.ID) #22

#superexact Test for overlap between GO sets
ST_input <- list(Aurelia_sig$GO.ID, Sarsia_sig$GO.ID, Trip_sig$GO.ID)
#n = number of go terms that could be shared, AKA:
n=length(intersect(intersect(Aurelia_GO$GO.ID, Sarsia_GO$GO.ID), Trip_GO$GO.ID))
test <- supertest(ST_input, n=n)
test$P.value

###
# ##Calculate overlaps of DE OGs between taxa.

OrthoGroups <- read.csv(file="OGs_multicopy_parsed.csv")

res_aurelia <- na.omit(read.csv("Aurelia/Aurelia_DESeq_RM.csv"))
res_sarsia <- na.omit(read.csv("Sarsia/Sarsia_DESeq_TBM.csv"))
res_trip <- na.omit(read.csv("Tripedalia/Trip_Noiseq_res.csv"))

DEGs_a <- res_aurelia[res_aurelia$log2FoldChange>0 & res_aurelia$padj<0.05,]$X
DEGs_s <- res_sarsia[res_sarsia$log2FoldChange>0 & res_sarsia$padj<0.05,]$X
DEGs_t <- res_trip[res_trip$M>0 & res_trip$prob>0.99,]$X

OrthoGroups <- data.table(OrthoGroups)

Aurelia_sig_OGs <- unique(OrthoGroups[Gene %in% DEGs_a]$OG)
Sarsia_sig_OGs <- unique(OrthoGroups[Gene %in% DEGs_s]$OG)
Trip_sig_OGs <- unique(OrthoGroups[Gene %in% DEGs_t]$OG)

ST_input <- list(Aurelia_sig_OGs, Sarsia_sig_OGs, Trip_sig_OGs)
#n = number of OGs with at least one gene in all 3 species: 5966
n=5966
test <- supertest(ST_input, n=n)
p.adjust(test$P.value, method="BH")



