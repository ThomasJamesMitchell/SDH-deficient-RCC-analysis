#### Code to plot trinucleoptide context of samples to compare subclonal to clonal samples


### R-3.3.0 on cgpfoo
library("seqinr")
library("Rsamtools")
library("GenomicRanges")
genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"
bases <- c("A", "C", "G", "T")

setwd("/lustre/scratch116/casm/cgp/users/tjm/renal/SDH/DNA/tncontext/plots")

## 1A load mutations
patients <- substr(list.files("/lustre/scratch116/casm/cgp/users/tjm/renal/SDH/DNA/dpclust/data/output"), 1, 7)
muts <- read.delim("/lustre/scratch116/casm/cgp/users/tjm/renal/SDH/DNA/muts/data/allMutsPassedFlags.txt", skip=116, stringsAsFactors = F)
samples = c("PD47450a", "PD47450c", "PD47450d", "PD47450e", "PD47451a", "PD47453a", "PD47454c")
muts <- muts[muts$ASMD>=140 & muts$CLPM<1,]
muts <- muts[muts$Sample %in% samples,]
### count each mutation only once
muts <- muts[!duplicated(muts$Pos),]

purity = NULL
for (sample in samples) {
  purity <- rbind(purity, data.frame("Sample"=sample, "cellularity"=read.delim(paste0("../../battenberg/data/",sample,".battenberg.rho_and_psi.txt"))[2,1]))
}

#### normalise VAF to account for purity
muts$cellularity <- purity$cellularity[match(muts$Sample, purity$Sample)]
muts$ccf = as.numeric(muts$PM.Tum)/muts$cellularity*2
muts$donor <- substr(muts$Sample, 1, 7)

#### split into low and high VAF for clonal and subclonal. Based on histograms, a clonality split at 50% is reasonable
muts$subclonal <- muts$ccf <0.5

## 1. Loading the mutational signatures from Ludmil
##    This generates the 2 necessary variables:
##    - signatures: numeric matrix of the probs of each mutation type under each signature
##    - muttype: index vector to classify substitutions in the 96 channels
mutlist = paste0(rep(rep(bases, times=6), each=4),
                 rep(c("C", "T"), each=48),
                 rep(bases, times=24),">",rep(rep(bases, times=6), each=4),
                 rep(c("A", "G","T","A","C","G"), each=16),
                 rep(bases, times=24))

oppstrand = function(x) {    
  trin1 = paste(comp(rev(strsplit(substr(x,1,3),"")[[1]]), forceToLower=F), collapse="")
  trin2 = paste(comp(rev(strsplit(substr(x,5,7),"")[[1]]), forceToLower=F), collapse="")
  return(sprintf("%s>%s",trin1,trin2))
}
mutlist_oppstrand = sapply(mutlist, oppstrand)

# Indexing vector to classify mutations in the 96 channels
muttype = rep(1:96,2)
names(muttype) = c(mutlist,mutlist_oppstrand)



## 2. Loading the mutations
mutations = muts[,c("Chrom", "Pos", "Ref", "Alt","subclonal", "Sample")]

# Annotating the trinucleotide context
seqs = scanFa(genomeFile, GRanges(mutations$Chrom, IRanges(mutations$Pos-1, mutations$Pos+1)))
ref_trinuc = as.vector(seqs)
mutations$trinuc_sub = paste(ref_trinuc, paste(substr(ref_trinuc,1,1), mutations$Alt, substr(ref_trinuc,3,3),sep=""), sep=">")
mutations$channel = muttype[mutations$trinuc_sub]


## 3.  Plotting the trinucleotide context
cols <- rep(c("blue", "black", "red", "grey", "green", "pink"), each=32)
# only include sigs that contribute > x %

pdf(paste("clonalVsSubclonal_TN_context.pdf", sep=""), width = 12, height = 3.6, useDingbats=F)
barplot(rbind(tabulate(mutations$channel[!mutations$subclonal])/sum(!mutations$subclonal), tabulate(mutations$channel[mutations$subclonal])/sum(mutations$subclonal)),
        col=rep(c("blue", "black", "red", "grey", "green", "pink"), each=32), border=T, ylab="Propertion of mutations", xlab="", xaxt="n",
        main="", beside=T)
mtext(mutlist, side=1, las=2, at=(1:96)*3, cex=0.6, padj=-.1)
dev.off()

