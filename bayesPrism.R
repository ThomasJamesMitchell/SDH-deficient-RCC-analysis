
####################################################################################
##################### BayesPrism deconvolution #####################
############## uses method from https://www.nature.com/articles/s43018-022-00356-3 ##########
####################################################################################

setwd("~/renal/SDH/RNA/deconvolve/")

library(TED)
library(Matrix)

#load("scripts/gbm.rdata")

#### which reference to use
ref = "YoungRef" #"CSF"   #

####  load bulk SDH sequencing data first
bulkRNASDH <- read.delim("data/bulkRNA/allSDH_renal.tsv")
colnames(bulkRNASDH) <- sub(".converted_counts", "", colnames(bulkRNASDH))
bulkRNASDH <- t(bulkRNASDH)

####  load bulk TCGA sequencing data first
bulkRNA <- read.delim("data/bulkRNA/allTCGA_renal.tsv")
bulkRNA <- t(bulkRNA)

####  merge SDH and TCGA datasets
bulkRNASDH <- bulkRNASDH[,match(colnames(bulkRNA), colnames(bulkRNASDH))]
bulkRNA <- bulkRNA[,match(colnames(bulkRNASDH), colnames(bulkRNA))]
bRNA <- rbind(bulkRNASDH, bulkRNA)

#### load reference scRNA seq data
if (ref %in% "YoungRef") {
  scRef <- readMM("data/scData/scSigsSymbol.mtx")
  scRef <- as.matrix(scRef)
  colnames(scRef) <- read.table("data/scData/scSigsSymbol_columnNames.tsv")[[1]]
  rownames(scRef) <- read.table("data/scData/scSigsSymbol_rowNames.tsv")[[1]]
}
if (ref %in% "CSF") {
  scRef <- readMM("data/scData/csfAll.mtx")
  scRef <- as.matrix(scRef)
  colnames(scRef) <- read.table("data/scData/csfAll_columnNames.tsv")[[1]]
  rownames(scRef) <- read.table("data/scData/csfAll_rowNames.tsv")[[1]]
}

#### make manifest and load signature file
mani = data.frame(sub(":.*","",colnames(scRef)))
rownames(mani) <- colnames(scRef)
if (ref %in% "YoungRef") {sigs <- read.delim("data/signatureKey2.txt")}
if (ref %in% "CSF") {sigs <- read.delim("data/signatureKeyCSF.txt")}
mani$cell.type <- sigs[match(mani[,1],sigs[,1]),2]
mani$cell.subtype <- sigs[match(mani[,1],sigs[,1]),3]
mani$patient.id <- sub("___.*","",sub(".*:","",rownames(mani)))

#### filter mito, sex and lowly expressed genes as recommended
ref.dat.filtered <- cleanup.genes(ref.dat=t(scRef),
              species = "hs",
              gene.type = c("RB","chrM","chrX","chrY"),
              input.type = "scRNA",
              exp.cells = 5)

rm(list=setdiff(ls(), c("ref.dat.filtered", "bRNA","mani","ref")))

#### run
tedOut = run.Ted(ref.dat=ref.dat.filtered,
        X=bRNA,
        cell.type.labels = as.character(mani$cell.type),
        cell.subtype.labels = as.character(mani$cell.subtype),
        input.type = "scRNA",
        n.cores=4,
        pdf.name="SDH_prism")

#### save
saveRDS(tedOut, paste0("output/prism/SDHTCGArenal",ref,".RDS"))
write.table(tedOut$res$first.gibbs.res$gibbs.theta, paste0("output/prism/SDHTCGArenal",ref,"Subtypes.tsv"), quote=F, sep="\t")
write.table(tedOut$res$final.gibbs.theta, paste0("output/prism/SDHTCGArenal",ref,"Types.tsv"), quote=F, sep="\t")

####################################################################################
#### plot and interrogate
####################################################################################
ref="YoungRef"#"CSF"#
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(reshape2)
library(ggsci)

setwd("~/renal/SDH/RNA/deconvolve/")

#### load bayesprism output
exposures <- t(read.delim(paste0("output/prism/SDHTCGArenal",ref,"Types.tsv")))

####################
# Load sample meta-data #
####################
mani <- read.delim("../../../RNAseqOutcomes/data/full_mainfest_all_tumours.txt")
mani <- mani[,c(23, 255, 126, 74, 108, 198, 194, 94, 95)]
histo <- read.delim("../../../RNAseqOutcomes/data/updated_histology.txt")   ##### This is the updated histology file
mani$type <- mani$gdc_cases.project.name
levels(mani$type) <- c(levels(mani$type), "Normal", "ccRCC", "pRCC", "chRCC")
mani$type[grep("Normal", mani$gdc_cases.samples.sample_type)] <- "Normal"
mani$type[mani$type %in% "Kidney Renal Clear Cell Carcinoma"] <- "ccRCC"
mani$type[mani$type %in% "Kidney Renal Papillary Cell Carcinoma"] <- "pRCC"
mani$type[mani$type %in% "Kidney Chromophobe"] <- "chRCC"
mani$update <- histo$Additional.Pathology.Data[match(mani$xml_bcr_patient_barcode, histo$bcr_patient_barcode)]  ## update the histology
CIMP <- c("TCGA-A4-7915", "TCGA-BQ-5879", "TCGA-BQ-5893", "TCGA-BQ-5894", "TCGA-F9-A4JJ", "TCGA-G7-6793", "TCGA-GL-7966", "TCGA-P4-A5E8", "TCGA-P4-A5EA")
HLRCC <- c("TCGA-BQ-5879","TCGA-BQ-5893","TCGA-BQ-5894","TCGA-F9-A4JJ")
mani$update[mani$xml_bcr_patient_barcode %in% CIMP] <- "CIMP" 
mani$update[mani$xml_bcr_patient_barcode %in% HLRCC] <- "HLRCC" 
mani$update[mani$gdc_cases.project.name %in% "Pheochromocytoma and Paraganglioma"] <- "PCPG"
levels(mani$update) <- c(levels(mani$update), "Normal", "ccRCC", "pRCC", "pRCC1", "pRCC2", "pRCCu", "chRCC", "MD-chRCC")
mani$update[mani$type %in% "Normal"] <- "Normal"
mani <- mani[!is.na(mani$update),]  ### for those without updated path we will use the original. Can remove later instead of this assumption
mani$update[mani$update %in% "Clear cell RCC"] <- "ccRCC"
mani$update[mani$update %in% "ChRCC"] <- "chRCC"
mani$update[mani$update %in% "Metabolically Divergent (MD-)ChRCC "] <- "MD-chRCC"
mani$update[mani$update %in% "Type 1 Papillary RCC"] <- "pRCC1"
mani$update[mani$update %in% "Type 2 Papillary RCC"] <- "pRCC2"
mani$update[mani$update %in% "Unclassified Papillary RCC"] <- "pRCCu"
mani$update[mani$update %in% "pRCC"] <- "pRCCu"
mani[,1][mani[,4] %in% "Pheochromocytoma and Paraganglioma"] <- mani[,2][mani[,4] %in% "Pheochromocytoma and Paraganglioma"]
mani[,1] <- gsub("-", ".", toupper(mani[,1]))

colnames(exposures) <- sub("TCGA_", "", colnames(exposures))
colnames(exposures)[grep("TCGA",colnames(exposures))] <- substr(colnames(exposures)[grep("TCGA",colnames(exposures))],1,12)
#colnames(exposures) <- paste0(mani$update[match(colnames(exposures), mani[,1])], "___", colnames(exposures))
colnames(exposures) <- mani$update[match(colnames(exposures), mani[,1])]
colnames(exposures)[1:10] <- c(rep("SDH_renal",6),"ccRCC","SDH_renal","PCPG","SDH_renal")
exposures <- exposures[,!is.na(colnames(exposures))]

table(colnames(exposures), useNA="always")

#### we won't use a heatmap of all celltypes as this is not looking at cell of origin and a comprehensive cancer cell single cell reference is not
#### available for all tumour types
colMain = colorRamp2(c(0,1),c('white','black'))
hmdata <- exposures
topAnno = HeatmapAnnotation(Histology=colnames(hmdata))
hm = Heatmap(hmdata,
             col=colMain,
             name='Exposure',
             top_annotation = topAnno,
             #             right_annotation = sideAnno,
             cluster_rows = FALSE,
             row_gap = unit(2.0,'mm'),
             #             row_order = order(sMani$over, sMani$Alias),
             #             row_labels=sMani$Alias,
             #             row_split = sMani$Key,
             column_split = paste0(colnames(hmdata)),
             column_gap = unit(3.0,'mm'),
             show_column_name=F,
             # column_names_max_height = unit(5,'in'),
             column_labels = F
)

if (ref %in% "YoungRef") {out <- melt(exposures[c("Leukocytes","Myeloid","Endothelial"),])
lims = c(0.024,0.15)}
if (ref %in% "CSF") {out <- melt(exposures[c("Lymphoid","Myeloid","Endothelial"),])
lims=c(c(0.08,0.15))}
colnames(out) <- c("cellType","Tumour", "Contribution")
out <- out[!out$Tumour %in% "PCPG",]
out$Tumour <- factor(out$Tumour, levels=c("SDH_renal", "ccRCC", "pRCC1", "pRCC2", "pRCCu", "CIMP", "HLRCC","chRCC","MD-chRCC","Normal"))
pdf(paste0("plots/bayesPrismFineAnnot",ref,".pdf"), width=3.5,height=3.5)
for (type in unique(out$cellType)) {
  dat <- out[out$cellType %in% type,]
  print(ggplot(dat, aes(x=Tumour, y=Contribution, color=Tumour)) +
          ggtitle(type)+
          geom_boxplot(outlier.shape = NA) +
          theme_classic()+ scale_color_brewer(palette = "Paired") +
          geom_jitter(position=position_jitter(0.2), alpha=0.2, cex=0.5) +
#          scale_y_continuous(trans = 'log10')
          coord_cartesian(xlim = NULL, ylim = c(0,lims[unique(out$cellType) %in% type]),expand = TRUE,default = FALSE,clip = "on") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  )
  print(pairwise.wilcox.test(dat$Contribution, dat$Tumour,
                             p.adjust.method = "BH"))
  
}
dev.off()

