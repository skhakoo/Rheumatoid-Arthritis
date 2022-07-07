library(affy)
library(affxparser)
library(oligo)
library(limma)
library(dplyr)
library(pd.hta.2.0)
library(arrayQualityMetrics)
library(affycoretools)
library(caret)
library(ggplot2)
library(hta20transcriptcluster.db)


#importing CEL files
unzip("VEDERA_Raw.zip")
gc()
basedir<-"/mnt/jw01-aruk-home01/projects/moonshot/sophie/data/VEDERA_Raw/Cel_files"
setwd(basedir)
CD19 = VEDERA_CD19_Clinical$Sample_id
list(CD19)
data = read.celfiles(CD19)
save(data, file = "/mnt/jw01-aruk-home01/projects/moonshot/sophie/output/CD19/CD19_HTA_Feature_Set.Rdata")
head(data)


#Background correction, normalization, and summarization
#Using rma
core_genes = rma(data, target="core")
save(core_genes, file = "/mnt/jw01-aruk-home01/projects/moonshot/sophie/output/CD19/rmaCD19_core_genes.Rdata")
gc()

probe_set = rma(data, target = "probeset")
save(probe_set, file = "/mnt/jw01-aruk-home01/projects/moonshot/sophie/output/CD19/rmaCD19_probeset.Rdata")



#Removing control probes - from core genes
load(paste0(path.package("pd.hta.2.0"), "/extdata/netaffxTranscript.rda"))
annot <- pData(netaffxTranscript)
t_annot<-annot[(annot$category=="main"),]
main_clusterID <- t_annot$transcriptclusterid
transcriptID <- as.data.frame(main_clusterID)

expr <- exprs(core_genes)
expr <- cbind(main_clusterID = rownames(expr), expr)
rownames(expr) <- 1:nrow(expr)
core_expr <- merge(expr, transcriptID, by="main_clusterID")
rownames(core_expr) <- core_expr$main_clusterID
core_expr <- within(core_expr, rm(main_clusterID))

save(core_expr, file="CD19_core_expr.Rdata")
gc()

#Remove control probes from probeset genes
load(paste0(path.package("pd.hta.2.0"), "/extdata/netaffxProbeset.rda"))
annot <- pData(netaffxProbeset)
p_annot<-annot[(annot$probesettype=="main"),]
main_clusterID <- p_annot$probesetid
probesetID <- as.data.frame(main_clusterID)

expr <- exprs(probe_set)
expr <- cbind(main_clusterID = rownames(expr), expr)
rownames(expr) <- 1:nrow(expr)
probe_set <- merge(expr, probesetID, by="main_clusterID")
rownames(probe_set) <- probe_set$main_clusterID
probe_set <- within(probe_set, rm(main_clusterID))

save(probe_set, file="CD19_probe_set.Rdata")



#PCA using prcomp()
df <- VEDERA_CD19_Clinical
core_pca <- prcomp(t(core_expr), scale=TRUE)
head(core_pca$x[,1:5])

df$PC1 <- core_pca$x[,1]
df$PC2 <- core_pca$x[,2]
df$Responder_wk24[df$Responder_wk24 == 0] <- "poor"
df$Responder_wk24[df$Responder_wk24 == 1] <- "good"
df <- df[complete.cases(df[ ,14]),]

df_wk0<-df[!(df$Time=="wk24"),]
df_wk24<-df[!(df$Time=="wk00"),]

all_pca <- ggplot(df, aes(PC1, PC2)) +
  geom_point(aes(colour = Responder_wk24, shape=Time)) +
  ggtitle("CD19 - Response PCA") +
  theme_bw() +
  theme(panel.grid.major = element_blank())
all_pca
ggsave("CD_19_all_pca.pdf")
ggsave(all_pca, plot=lastplot())

wk0_pca <- ggplot(df_wk0, aes(PC1, PC2, label=Patient)) +
  geom_point(aes(colour = Responder_wk24)) +
  geom_text(hjust=1, vjust=1, size=2)+
  ggtitle("CD19 - Response Week 0") +
  theme_bw() +
  theme(panel.grid.major = element_blank())
wk0_pca
ggsave("CD19_wk0_pca.pdf")

wk24_pca <- ggplot(df_wk24, aes(PC1, PC2, label=Patient)) +
  geom_point(aes(colour = Responder_wk24)) +
  geom_text(hjust=1, vjust=1, size=2)+
  ggtitle("CD19 - Response Week 24") +
  theme_bw() +
  theme(panel.grid.major = element_blank())
wk24_pca
ggsave("CD19_wk24_pca.pdf")


#LIMMA model
df <- VEDERA_CD19_Clinical
rm(VEDERA_CD19_Clinical)
gc()

df$Time[df$Time== "wk00"] <- 0
df$Time[df$Time== "wk24"] <- 6
df$Responder_wk24[df$Responder_wk24 == 0] <- "NR"
df$Responder_wk24[df$Responder_wk24 == 1] <- "GR"
timeanalysis <- factor(paste(df$Responder_wk24, df$Time, sep="_"))

design <- model.matrix(~0+timeanalysis+df$Allocation+df$agebl+df$gender)
design


#array weights
aw <- arrayWeights(core_expr, design)
lm <- lmFit(core_expr, design, weights=aw)

#make contrasts
contrasts <- makeContrasts(Response_Baseline = GR_0 - NR_0,
                           Response_6months = GR_6 - NR_6,
                           DE_Poor = NR_6 - NR_0,
                           DE_Good = GR_6 - GR_0,
                           Overall_DE = ((GR_6 - GR_0)-(NR_6 - NR_0)), levels = design)


contr.fit <-contrasts.fit(lm, contrasts)
contr.fit <- eBayes(contr.fit, trend=TRUE)
save(contr.fit, file= "/mnt/jw01-aruk-home01/projects/moonshot/sophie/output/CD19_DGE_core.Rdata")

R_Baseline <- topTable(contr.fit, coef ="Response_Baseline", number = 10)
R_3months <- topTable(contr.fit, coef ="Response_3months", number = 10)
DE_Poor <- topTable(contr.fit, coef ="DE_Poor", number = 10)
DE_Good <- topTable(contr.fit, coef ="DE_Good", number = 10)
overall_DE <- topTable(contr.fit,coef ="Overall_DE", number=67528)

save(R_Baseline, file = "CD19_Baseline_DEG.Rdata")
save(R_3months, file = "CD19_3months_DEG.Rdata")
save(DE_poor, file = "CD19_Poor_DEG.Rdata")
save(DE_Good, file = "CD19_GOOD_DEG.Rdata")
save(overall_DE, file = "CD19_overallDGE.Rdata")


#Volcano plot
df$diffexpressed <- "NO"
df$diffexpressed[df$logFC > 0.3 & df$P.Value < 0.05] <- "UP"
df$diffexpressed[df$logFC < -0.3 & df$P.Value < 0.05] <- "DOWN"

df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- rownames(df)[df$diffexpressed != "NO"]


png("/mnt/jw01-aruk-home01/projects/moonshot/sophie/output/CD19_VolcanoPlot.png",  width = 15, height = 15, units = 'in', res = 300)
ggplot(data=df, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.3, 0.3), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
dev.off()


#Convert probe_ID to ensembl_ID and gene symbol)
df <- overall_DE
df <- cbind(probe_ID = rownames(df), df)
rownames(df) <- 1:nrow(df)
df$probe_ID<-gsub(".1","",as.character(df$probe_ID))

library(biomaRt)
listMarts()
ensembl <- useEnsembl('ensembl',dataset = 'hsapiens_gene_ensembl')
annot <- getBM(
  attributes = c(
    'affy_hta_2_0',
    'hgnc_symbol',
    'ensembl_gene_id',
    'entrezgene_id',
    'ensembl_transcript_id'),
  mart = ensembl)
head(annot)
colnames(annot)[which(names(annot) == "affy_hta_2_0")] <- "probe_ID"
DEGs <- merge(df, annot, by="probe_ID")
