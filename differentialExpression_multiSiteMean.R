
#### Tom Mitchell August 2020
#### Differential expression analysis comparing the mRNA from SDH renal with those from TCGA (KIRC/ KIRP/ KICH/ PGPL)
#### Aims to define important pathways upregulated, compare 
## setup on farm5 R-3.6
## now with a mean expression for our multisite sample to avoid distorting the effects

setwd("/lustre/scratch124/casm/team332/tjm/renal/SDH/RNA/diffExp")
library(DESeq2)
library(ggplot2)
library(tibble)
library(dplyr)

### load TCGA counts and manifest for kidney mRNA, remove paeds samples, and shorten sample and gene ENSG names
tcgaCounts = read.delim("/lustre/scratch124/casm/team332/tjm/renal/science1_scRNA/bulk_RNA/TCGA_data/all_HTseq_Counts.txt")
tcgaCounts = tcgaCounts[-1:-5,]  ### remove general feature information
rownames(tcgaCounts) = substr(rownames(tcgaCounts),1,15)
# tcgaCounts <- read.delim("/lustre/scratch116/casm/cgp/users/tjm/renal/science1_scRNA/bulk_RNA/TCGA_data/all_FPKM.txt")  
tcgaMani <- read.delim("/lustre/scratch124/casm/team332/tjm/renal/habitat/rna/cellSignalAnalysis/data/updated_histology.txt")
tcgaCounts <- tcgaCounts[,!grepl("TARGET", colnames(tcgaCounts))]
colnames(tcgaCounts) <- substr(colnames(tcgaCounts), 1, 16)
tcgaTypes <- tcgaMani$Additional.Pathology.Data[match(substr(colnames(tcgaCounts), 1, 12), gsub("-",".",tcgaMani$bcr_patient_barcode))]
tcgaCounts <- tcgaCounts[,!is.na(tcgaTypes)]
tcgaTypes <- tcgaTypes[!is.na(tcgaTypes)]

### load the count data from SDH renal sequencing into one dataframe
files = list.files("/lustre/scratch124/casm/team332/tjm/renal/SDH/RNA/counts/data", full.names = T)
sdhCounts <- read.delim(files[1], skip=6, row.names = 1)[,c(1,4)]
for (file in files) {
  sdhC <- read.delim(file, skip=6)[,6]
  sdhCounts <- cbind(sdhCounts, sdhC)
}
sdhCounts <- sdhCounts[,3:ncol(sdhCounts)]
colnames(sdhCounts) <- sub(".*data/","",sub(".con.*","",files))
multi <- round(rowMeans(sdhCounts[,grep("PR47450",colnames(sdhCounts))]))
sdhCounts <- sdhCounts[,!grepl("PR47450",colnames(sdhCounts))]
sdhCounts <- cbind("PR47450"= multi, sdhCounts)
sdhTypes <- rep("SDH renal", ncol(sdhCounts))
sdhTypes[grepl("PR47452a", colnames(sdhCounts))] <- "Clear cell RCC"
sdhTypes[grepl("PR47454a", colnames(sdhCounts))] <- "Paraganglioma"

### load the count data from previous clear cell sequencing into one dataframe
files = list.files("/lustre/scratch124/casm/team332/tjm/renal/habitat/rna/cellSignalAnalysis/bulkData/HABITATsamples", full.names = T)
ccCounts <- NULL
for (file in files) {
  ccC <- read.delim(file, row.names=1)[,2]
  ccCounts <- cbind(ccCounts, ccC)
}
colnames(ccCounts) <- sub(".*samples/","",sub(".tsv","",files))
rownames(ccCounts) <- read.delim(file[1])[,1]
ccTypes <- rep("Clear cell RCC", ncol(ccCounts))
ccTypes[grepl("b", colnames(ccCounts))] <- "Normal kidney"


# sdhbt <- tcgaCounts[grep("ENSG00000117118", rownames(tcgaCounts)),]
# sdhCounts[sdhCounts$gene %in% "SDHB",]

### merge the datasets after ordering the genes
sdhCounts <- sdhCounts[rownames(sdhCounts) %in% rownames(tcgaCounts),]
ccCounts <- ccCounts[rownames(ccCounts) %in% rownames(tcgaCounts),]
tcgaCounts <- tcgaCounts[rownames(tcgaCounts) %in% rownames(sdhCounts),]
allCounts <- cbind(tcgaCounts, sdhCounts, ccCounts)  ### merge the datasets
allCounts <- allCounts[,!duplicated(colnames(allCounts))]
allTypes <- c(as.character(tcgaTypes), sdhTypes, ccTypes)
sampleInfo <- data.frame("sample"=colnames(allCounts), "tumourType"= allTypes, stringsAsFactors = F)
sampleInfo$tumourType[grep("11A", sampleInfo$sample)] <- "Normal kidney"
sampleInfo$batch <- "sanger"
sampleInfo$batch[grep("TCGA", sampleInfo$sample)] <- "TCGA"
sampleInfo$isSDH <- sampleInfo$tumourType %in% "SDH renal"



#################################################################################################################################################
### differential expression analysis all samples
if (file.exists("DE.RData")) {load("DE.RData")}
if (!file.exists("DE.RData")) {
  sampleInfo$tumourType <- relevel(as.factor(sampleInfo$tumourType), "Normal kidney")  ### make normal kidney the first factor and therefore the intercept for the design matrix
  design = as.formula(~ isSDH + batch)
  modelMatrix <- model.matrix(design, sampleInfo)
  
  ddsObj.raw <- DESeqDataSetFromMatrix(countData = allCounts,
                                       colData = sampleInfo,
                                       design = design)
  
  vstcounts <- vst(ddsObj.raw, blind=TRUE)
  pcaData <- plotPCA(vstcounts, intgroup=c("tumourType"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  pdf("plots/PCA_all.pdf", width=7, height=5)
  ggplot(pcaData, aes(PC1, PC2, color=tumourType)) +
    geom_point(size=c(3, 0.4)[grepl("TCGA", rownames(pcaData))+1]) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    theme_classic()
  
  dev.off()
  
  
  ddsObj <- DESeq(ddsObj.raw)
  save(ddsObj, file="DE.RData")
  
  resultsNames(ddsObj)
  
}

resSDHall <- results(ddsObj, 
                    name="isSDHTRUE", 
                    alpha = 0.05)
topGenesSDHall <- as.data.frame(resSDHall) %>%
  rownames_to_column("GeneID") %>% 
  arrange(padj)
rownames(topGenesSDHall) <- topGenesSDHall$GeneID
geneMeta <- read.delim(list.files("/lustre/scratch116/casm/cgp/users/tjm/renal/SDH/RNA/counts/data", full.names = T)[1], skip=6, row.names = 1)[,c(1,2)]
topGenesSDHall <- merge(topGenesSDHall, geneMeta, by=0, sort=F)
topGenesSDHall <- topGenesSDHall[order(topGenesSDHall$log2FoldChange, decreasing = T),]
write.table(topGenesSDHall[,-1], "output/topGenesSDHSangerandTCGAAll.tsv", quote=F, sep="\t", row.names = F)
topGenesSDHall <- topGenesSDHall[topGenesSDHall$biotype %in% "protein_coding",]
write.table(topGenesSDHall[,-1], "output/topGenesSDHSangerandTCGA.tsv", quote=F, sep="\t", row.names = F)

topTCGAInformed <- topGenesSDHall[topGenesSDHall$gene %in% topGenesSDH$gene[1:2000],]
dat <- data.frame("tumourType" = as.character(sampleInfo[,2]), t(allCounts[topTCGAInformed$GeneID[1:20],]))
pdf("plots/markerGenesCombined4.pdf", width=5, height=6)
for (i in 1:20) {
  print(ggplot(dat, aes_string(x="tumourType", y=colnames(dat)[i+1])) +
          geom_violin() +
          theme_classic() +
          scale_y_continuous(trans='log10') +
          theme(axis.text.x=element_text(angle=90,hjust=1)) +
          geom_jitter(shape=16, position=position_jitter(0.2)) +
          ggtitle(topTCGAInformed$gene[i]))
}
dev.off()

#################################################################################################################################################
### differential expression analysis all samples compared to normal only
if (file.exists("DESDHNorm.RData")) {load("DESDHNorm.RData")}
if (!file.exists("DESDHNorm.RData")) {
  sampleInfo$tumourType <- relevel(as.factor(sampleInfo$tumourType), "Normal kidney")  ### make normal kidney the first factor and therefore the intercept for the design matrix
  sampleInfo$isNormal <- sampleInfo$tumourType %in% "Normal kidney"
  design = as.formula(~ tumourType + batch)
  modelMatrix <- model.matrix(design, sampleInfo)
  
  ddsObj.raw <- DESeqDataSetFromMatrix(countData = allCounts,
                                       colData = sampleInfo,
                                       design = design)
  
  ddsObj <- DESeq(ddsObj.raw)
  save(ddsObj, file="DESDHNorm.RData")
  
  resultsNames(ddsObj)
  resSDHNorm <- results(ddsObj, name =  "tumourType_SDH.renal_vs_Normal.kidney")
  
  topGenesSDHnorm <- as.data.frame(resSDHNorm) %>%
    rownames_to_column("GeneID") %>% 
    arrange(padj)
  rownames(topGenesSDHnorm) <- topGenesSDHnorm$GeneID
  topGenesSDHnorm <- merge(topGenesSDHnorm, geneMeta, by=0, sort=F)
  topGenesSDHnorm <- topGenesSDHnorm[order(topGenesSDHnorm$log2FoldChange, decreasing = T),]
  topGenesSDHnorm <- topGenesSDHnorm[topGenesSDHnorm$biotype %in% "protein_coding",]
  write.table(topGenesSDHnorm[,-1], "output/topGenesSDHSanger4andTCGAcompNormal.tsv", quote=F, sep="\t", row.names = F)
  
  dat <- data.frame("tumourType" = as.character(sampleInfo[,2]), t(allCounts[topGenesSDHnorm$GeneID[1:20],]))
  
  keys=factor(c("pRCC1","pRCC2","pRCCu","ccRCC","Normal","chRCC","MD-chRCC","SDH_renal","PCPG"), levels=c("SDH_renal", "ccRCC", "pRCC1", "pRCC2", "pRCCu","chRCC","MD-chRCC","Normal"))
  names(keys)=c("Type 1 Papillary RCC","Type 2 Papillary RCC","Unclassified Papillary RCC","Clear cell RCC","Normal kidney","ChRCC","Metabolically Divergent (MD-)ChRCC","SDH renal","Paraganglioma")
  dat$tumourType <- keys[dat$tumourType]
  pdf("plots/markerGenesCombinedSDHNormComp4.pdf", width=5, height=6)
  for (i in 1:20) {
    print(ggplot(dat, aes_string(x="tumourType", y=colnames(dat)[i+1], color="tumourType")) +
            geom_boxplot(outlier.shape = NA) +
            theme_classic() +
            scale_y_continuous(trans='log10') +
            theme(axis.text.x=element_text(angle=90,hjust=1)) +
            geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.2) +
            ggtitle(topGenesSDHnorm$gene[i]))
  }
  dev.off()
}



####### plot genes of interest
genes <- c("PAPPA2","CD274","TGFB1","VEGFA", "MTOR", "MET", "JUN", "PDCD1", "CTLA4", "RET", "HRAS", "MAPK1")
ens <- rownames(geneMeta)[match(genes, geneMeta$gene)]
dat <- data.frame("tumourType" = as.character(sampleInfo[,2]), t(allCounts[ens,]))
keys=factor(c("pRCC1","pRCC2","pRCCu","ccRCC","Normal","chRCC","MD-chRCC","SDH_renal","PCPG"), levels=c("SDH_renal", "ccRCC", "pRCC1", "pRCC2", "pRCCu","chRCC","MD-chRCC","Normal"))
names(keys)=c("Type 1 Papillary RCC","Type 2 Papillary RCC","Unclassified Papillary RCC","Clear cell RCC","Normal kidney","ChRCC","Metabolically Divergent (MD-)ChRCC","SDH renal","Paraganglioma")
colnames(dat)[2:ncol(dat)] <- genes
dat$tumourType <- keys[match(dat$tumourType,names(keys))]
dat <- dat[!is.na(dat$tumourType),]
write.table(dat, "output/markerGeneCounts.tsv", quote=F, sep="\t")
dat <- read.delim("C:/Users/tjm/Documents/renal/SDH/RNA/diffExp/output/markerGeneCounts.tsv")
dat$tumourType <- factor(dat$tumourType, levels=c("SDH_renal","ccRCC","pRCC1","pRCC2","pRCCu","chRCC","Normal"))
pdf("C:/Users/tjm/Documents/renal/SDH/RNA/diffExp/plots/manualGenes.pdf", width=3.5, height=3.5)
for (i in 1:12) {
  ens <- colnames(dat)[i+1]
  print(ggplot(dat, aes_string(x="tumourType", y=ens, color="tumourType")) +
          geom_boxplot(outlier.shape = NA) +
          theme_classic()  + scale_color_brewer(palette = "Dark2") +
          scale_y_continuous(trans='log10') +
          theme(axis.text.x=element_text(angle=90,hjust=1)) +
          geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.2, cex=0.5) +
          ggtitle(ens))
}
dev.off()


pdf("plots/manualGenes.pdf", width=3, height=3)
for (i in 1:length(genes)) {
  ens <- rownames(geneMeta)[geneMeta$gene %in% genes[i]]
  print(ggplot(dat, aes_string(x="tumourType", y=ens)) +
          geom_violin() +
          theme_classic() + scale_fill_brewer(palette = "Dark2") +
          scale_y_continuous(trans='log10') +
          theme(axis.text.x=element_text(angle=90,hjust=1)) +
          geom_jitter(shape=16, position=position_jitter(0.2)) +
          ggtitle(genes[i]))
}
dev.off()


