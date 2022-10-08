######### Analysis of VAFs for SNV to determine clonality of samples
######### Updated analysis accounts for ploidy of tumour samples

setwd("~/renal/SDH/DNA/muts")

library(ggridges)
library(ggplot2)
library(dplyr)


file = "data/allMutsPassedFlags.txt"
strs <- readLines(file)
muts <- read.delim(file, skip=max(grep("##",strs)), header=T)
muts <- muts[muts$ASMD>=140 & muts$CLPM<1,]

#### collate purity estimates
samples = c("PD47450a", "PD47450c", "PD47450d", "PD47450e", "PD47451a", "PD47453a", "PD47454c")
purity = NULL
cn <- read.delim("../dpclust/data/CNsegments.txt")
muts$cn <- "diploid"
for (sample in samples) {
  purity <- rbind(purity, data.frame("Sample"=sample, "cellularity"=read.delim(paste0("../battenberg/data/",sample,".battenberg.rho_and_psi.txt"))[2,1]))
  cns <- cn[cn$Sample %in% sample,]
  for (i in 1:nrow(cns)) {
    muts$cn[muts$Sample %in% sample & muts$Chrom %in% cns[i,"chr"] & muts$Pos > cns[i,"lower"] & muts$Pos < cns[i,"upper"]] <- cns[i,"event"]
  }
}

#### restrict to SDH renal tumours only
muts <- muts[muts$Sample %in% samples,]

#### normalise VAF to account for purity
muts$cellularity <- purity$cellularity[match(muts$Sample, purity$Sample)]
muts$donor <- substr(muts$Sample, 1, 7)


muts$ccf = NULL
muts$ccf[muts$cn %in% "diploid"] =  as.numeric(muts$PM.Tum[muts$cn %in% "diploid"])/muts$cellularity[muts$cn %in% "diploid"]*2
muts$ccf[muts$cn %in% "loh"] =  as.numeric(muts$PM.Tum[muts$cn %in% "loh"])/muts$cellularity[muts$cn %in% "loh"]*1
muts$ccf[muts$cn %in% "gain"] =  as.numeric(muts$PM.Tum[muts$cn %in% "gain"])/muts$cellularity[muts$cn %in% "gain"]*3
muts$ccf[muts$cn %in% "gain"][muts$ccf[muts$cn %in% "gain"] > 1.5] <- muts$ccf[muts$cn %in% "gain"][muts$ccf[muts$cn %in% "gain"] > 1.5]/2
muts <- muts[!muts$Chrom %in% c("X","Y"),]

p <- ggplot(muts, aes(x = ccf, y=Sample, fill=donor)) +
  geom_density_ridges() + scale_fill_brewer(palette = "Dark2") +
  theme_minimal() + theme_classic()+xlab("Cancer cell fraction")+ylab("Density")+ theme(legend.position = "none") +coord_cartesian(clip = 'off') +
  scale_x_continuous(limits=c(0,3.5)) + 
  theme(axis.text.y=element_blank()) +
  annotate("text",label = paste0(unique(muts$Sample), ", purity = ",round(unique(muts$cellularity),2)), x=3, y=(1:7)+0.7)

pdf("plots/VAFHistogramsSDH.pdf", width=5, height=3)
p
dev.off()

