
#### Code to time 1p loss (via 1q gain), and MRCA of multiregion sample

setwd("~/renal/SDH/DNA")

library(lme4)
patients <- unique(substr(list.files("dpclust/data/output"), 1, 7))

####### we need to define non-shared mutations between patients as we assume the shared muts are contamination/ artefact
muts <- read.delim("dpclust/data/SDHmutCount.txt", stringsAsFactors = F)
maxMuts <- NULL
for (patient in patients) {
  dat <- muts[,grep(patient, colnames(muts))]
  if (!is.null(dim(dat))) {max <- apply(dat, 1, max)}
  if (is.null(dim(dat))) {max <- dat}
  maxMuts <- cbind(maxMuts, max)
}
colnames(maxMuts) <- patients
rownames(maxMuts) <- muts[,1]
shared <- rowSums(maxMuts>0)>1
sharedPos <- as.integer(sub(".*_","",names(shared)[shared]))

####### load collated caveman output and do usual filtering before removing shared variants
###then calculate the number of mutations called treating each sample independantly 
data <- read.delim("muts/data/allMutsPassedFlags.txt", skip=116, stringsAsFactors = F)
data <- data[data$Type %in% "Sub" & as.numeric(data$CLPM) == 0 & as.numeric(data$ASMD) >140,]
data <- data[!data$Pos %in% sharedPos,]
nms <- data.frame(table(data$Sample))

###### load ages
ages <- read.delim("muts/data/ages.txt")
nms$age <- ages$Age[match(substr(nms[,1], 1, 7), ages$SangerID)]


###### for the 2 multi-regional samples take the dpclust output to calulate the true branch lengths
nmc <- NULL
for (patient in c("PD47450")) {
  no <- read.delim(paste0("dpclust/output/", patient, "/DP_output/", patient, "_optimaInfo_0.01.txt"),
                   stringsAsFactors = F)
  no <- no[no$Terminal %in% c("No", "Yes"),]
  no$id <- gsub(".*>","",no$Name)
  no <- no[order(nchar(no$Name)),]  ### reverse order
  no$parent <- gsub(".*>","",substr(no$Name, 1, nchar(no$Name)-3))
  if (sum(duplicated(no$Name))>0) {
    dupl <- which(duplicated(no$Name) | duplicated(no$Name, fromLast = T))
    addn <- sum(no$no.of.mutations.assigned[dupl])
    no$no.of.mutations.assigned[dupl[1]] <- addn
    no <- no[-dupl[2:length(dupl)],]
  }
  no$nmuts <- 0
  no$nmuts[1] <- no$no.of.mutations.assigned[1]
  for (i in 2:nrow(no)) {
    no$nmuts[i] <- no$no.of.mutations.assigned[i] + no$nmuts[no$id %in% no$parent[i]]
  }
  no$PDbranch <- substr(colnames(no)[grep("PD",colnames(no))][apply(no[,grep("PD",colnames(no))], 1, which.max)], 8, 20)
  nmc <- rbind(nmc, no[no$Terminal %in% "Yes", c("PDbranch", "nmuts")])
}
nmc$age <- ages$Age[match(substr(nmc[,1], 1, 7), ages$SangerID)]

###### make estimates of rate
nmsr <- nms[!nms[,1] %in% c("PD47452a", "PD47454a"),]
nmsr$patient <- substr(nmsr[,1], 1, 7)
muts.per.year.lmer <- lmer(Freq ~ age - 1 + (age - 1 | patient), data=nmsr, REML=FALSE)
summary.lmer <- summary(muts.per.year.lmer)

nmsr$Scaled.num.muts <- scale(nmsr$Freq, center=FALSE)
scale.factor <- attr(nmsr$Scaled.num.muts, which = "scaled:scale")
age.lmer <- lmer(age ~ Scaled.num.muts - 1 + (Scaled.num.muts - 1 | patient), data=nmsr)
nmsr$Pt.fitted.num.muts <- scale.factor / (ranef(age.lmer)$patient[nmsr$patient,] + fixef(age.lmer)["Scaled.num.muts"]) * nmsr$age


####### plot age vs n muts
pdf("dpclust/plots/agesMuts.pdf", height=5, width=5)
plot(NA, xlim=c(10,max(nms$age)), ylim=c(0,max(nms[,2])),
     xlab="Age, yrs", ylab="Mutations, n", bty="n")

polygon(x=c(0,70,70,0), y=c(0,70*(summary.lmer$coefficients[1]+1.96*summary.lmer$coefficients[2]), 70*(summary.lmer$coefficients[1]-1.96*summary.lmer$coefficients[2]),0),
        col="light blue", border=NA)
abline(0, summary.lmer$coefficients[1], lwd=3, col="black")
points(nms$age[!nms[,1] %in% c("PD47452a", "PD47454a")], nms[!nms[,1] %in% c("PD47452a", "PD47454a"),2], pch = 16)
points(nmc$age[!nmc[,1] %in% c("PD47452a", "PD47454a")], nmc[!nmc[,1] %in% c("PD47452a", "PD47454a"),2], pch = 17, bty="n")
points(nms$age[nms[,1] %in% "PD47452a"], nms[nms[,1] %in% "PD47452a",2], pch = 16, col="Green")
points(nms$age[nms[,1] %in% "PD47454a"], nms[nms[,1] %in% "PD47454a",2], pch = 16, col="Red")
points(nmc$age[nmc[,1] %in% "PD47454a"], nmc[nmc[,1] %in% "PD47454a",2], pch = 17, col="Red")
legend("bottomright", legend = c("Estimate from single region", "Estimate from mutational clusters", "SHD-renal", "ccRCC", "Paraganglioma"),
       pch=c(1, 2, 16, 16, 16), col = c(rep("black", 3), "Green", "Red"))
dev.off()



########## plot phylogenies
pdf(paste("dpclust/plots/phylos.pdf", sep=""), width = 3, height = 5, useDingbats=F)
for (patient in c("PD47450")) {
  no <- read.delim(paste0("dpclust/output/", patient, "/DP_output/", patient, "_optimaInfo_0.01.txt"),
                   stringsAsFactors = F)
  no <- no[no$Terminal %in% c("No", "Yes"),]
  no$id <- gsub(".*>","",no$Name)
  no <- no[order(nchar(no$Name), no$Name),]
  no$parent <- gsub(".*>","",substr(no$Name, 1, nchar(no$Name)-3))
  if (sum(duplicated(no$Name))>0) {
    dupl <- which(duplicated(no$Name) | duplicated(no$Name, fromLast = T))
    addn <- sum(no$no.of.mutations.assigned[dupl])
    no$no.of.mutations.assigned[dupl[1]] <- addn
    no <- no[-dupl[2:length(dupl)],]
  }
  no$nmuts <- 0
  no$nmuts[1] <- no$no.of.mutations.assigned[1]
  for (i in 2:nrow(no)) {
    no$nmuts[i] <- no$no.of.mutations.assigned[i] + no$nmuts[no$id %in% no$parent[i]]
  }
  no$PDbranch <- substr(colnames(no)[grep("PD",colnames(no))][apply(no[,grep("PD",colnames(no))], 1, which.max)], 8, 20)
  if (patient %in% "PD47450") {no$x <- c(0, -.05, 1, -.6, -.02, -1, -.6, -.2, .2)}
  if (patient %in% "PD47454") {no$x <- c(0, -1, .5, .1, .9)}
  
  plot(x=NA,y=NA, xaxt="n", ylab="Mutations, n",xlab="", xlim=c(-1, 1), cex=0.6, ylim = c(0,700), bty="n")
  segments(no$x[1],no$nmuts[1],0,0, lwd=6, col="grey")
  for (j in 2:nrow(no)) {
    segments(no$x[j],no$nmuts[j],no$x[no$id %in% no$parent[j]],no$nmuts[no$id %in% no$parent[j]], lwd=6, col="grey")
  }
  symbols(no$x, no$nmuts,circles=rep(0.1, nrow(no)), add=TRUE, inches=0.03,
          fg="dark grey", bg="dark grey")
  symbols(no$x[1],no$nmuts[1],circles=0.1, add=TRUE, inches=0.03,
          fg="red", bg="red")
  mtext(patient, 1)
  
}
dev.off()


###############################################################################################
###############################################################################################
######################################  MRCA timing ######################################
###############################################################################################
# Generate point estimates for the timing of the MRCA emergence
new.dat <- data.frame(Sample=substr(nmsr$Var1[1],1,7), Scaled.num.muts=68 / scale.factor, Age=nmsr$age[1])
new.dat$MRCA.pred.age <- predict(age.lmer, newdata = new.dat)

# Now generate 95% CIs on the predictions for MRCA timing
# Note that CIs for lme models are challenging - bootstrapping appears to be most robust
mySumm <- function(.) {
  predict(., newdata=new.dat, re.form=NULL)
}

####Collapse bootstrap into median, 95% PI
sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

boot1 <- bootMer(age.lmer, mySumm, nsim=1000, use.u=FALSE, type="parametric")

PI.boot1 <- sumBoot(boot1)
PI.boot1$Age <- new.dat$Age
PI.boot1$pred.fit <- new.dat$MRCA.pred.age
PI.boot1$Time.lag <- PI.boot1$Age - PI.boot1$fit
PI.boot1$Time.lag[PI.boot1$Time.lag < 0] <- PI.boot1$Age[PI.boot1$Time.lag < 0] - PI.boot1$pred.fit[PI.boot1$Time.lag < 0]
PI.boot1$lag.lwr <- PI.boot1$Age - PI.boot1$upr
PI.boot1$lag.lwr[PI.boot1$lag.lwr < 0] <- 0
PI.boot1$lag.upr <- PI.boot1$Age - PI.boot1$lwr
PI.boot1$ID <- new.dat$Sample
PI.boot1 <- PI.boot1[order(PI.boot1$Time.lag),]

plot(PI.boot1$Time.lag, 1:nrow(PI.boot1), pch=20, xlim=c(-10,max(PI.boot1$lag.upr)), axes=FALSE, xlab="Estimated time between MRCA and diagnosis (years)", ylab="")
segments(x0 = PI.boot1$lag.lwr, x1 = PI.boot1$lag.upr, y0 = 1:nrow(PI.boot1))
axis(side=1, at=(0:5)*10, labels=(0:5)*10)
text(x = PI.boot1$lag.lwr-1, y = 1:nrow(PI.boot1), labels = PI.boot1$LRIID, adj = 1, cex=0.75)

mySumm.muts.per.yr <- function(.) {
  unlist(fixef(.) + ranef(.)$patient[,1]) 
}

boot2 <- bootMer(muts.per.year.lmer, mySumm.muts.per.yr, nsim=1000, use.u=TRUE, type="parametric")
colnames(boot2$t) <- row.names(ranef(muts.per.year.lmer)$patient)

###############################################################################################
###############################################################################################
################### 1p loss by gain timing ######################################
###############################################################################################
# First we define a function that reads in the relevant mutations and decides whether a variant is clonal or subclonal; and if clonal whether it occurs pre- or post-duplication

sa.df <- read.table("../../data/summary/sample_summary.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
sa.df <- sa.df[sa.df$include == "yes",]

mutCount <-  read.delim("dpclust/data/SDHmutCount.txt", row.names = 1)
wtCount <-  read.delim("dpclust/data/SDHWTCount.txt", row.names = 1)
cellularity <- read.delim("dpclust/data/SDHcellularity.txt", header=F)[[1]]
vaf <- t(t(mutCount/(wtCount+mutCount))/cellularity)*3
vaf <- vaf[!grepl("X|Y", rownames(vaf)),]
vaf <- vaf[!apply(vaf, 1, function(x) any(is.na(x))),]

timing <- function(input) {
  sample = input$Sample
  chr = input$chr
  lower = input$lower
  upper = input$upper
  
  pt = substr(sample, 1, 7)
  k = substr(sample, 8, 8)
  
  if (pt %in% "PD47450") {
    # Get the mutation data for that sample
    muts <- data.frame("vaf"=vaf[as.numeric(sub("_.*","",rownames(vaf))) == chr & as.numeric(sub(".*_","",rownames(vaf))) >= lower & as.numeric(sub(".*_","",rownames(vaf))) <= upper
                                 & rowSums(vaf[,substr(colnames(vaf),1,7)%in%pt])>0,
                                 grep(sample,colnames(vaf))])}
  
  if (!pt %in% "PD47450") {
    muts <- data.frame("vaf"=vaf[as.numeric(sub("_.*","",rownames(vaf))) == chr & as.numeric(sub(".*_","",rownames(vaf))) >= lower & as.numeric(sub(".*_","",rownames(vaf))) <= upper
                               & vaf[,colnames(vaf)%in%sample]>0,
                               grep(sample,colnames(vaf))])}
  
  # From CCF, infer whether mutation occurred pre- or post-duplication in that sample
  muts[,"event"] <- rep("Uninformative", nrow(muts))
  muts[muts$vaf > 1.5, "event"] <- "Clonal, duplicated"
  muts[muts$vaf <= 1.5, "event"] <- "non-duplicated"
  print(table(muts$event))
  
  
  if (pt %in% "PD47450") {
    clusterInfo <- cbind(read.delim("dpclust/output/PD47450/loci.txt"), read.delim("dpclust/output/PD47450/DP_output/PD47450_DP_and cluster_info_0.01.txt"))
    assignment <- read.delim("dpclust/output/PD47450/DP_output/PD47450_optimaInfo_0.01.txt")
    cluster <- clusterInfo$most.likely.cluster[match(as.numeric(sub(".*_","",rownames(muts))), clusterInfo$pos)]
    assign <- rep("Uninformative", length(cluster)) 
    assign[cluster %in% assignment[assignment[,grep(paste0(pt,k),colnames(assignment))]>0.2,"cluster.no"]] <- "Subclonal, present in this sample"
    assign[cluster %in% assignment[assignment[,grep(paste0(pt,k),colnames(assignment))]<0.2,"cluster.no"]] <- "Subclonal, not present in this sample"
    assign[cluster %in% 1] <- "Uninformative"
    assign[assignment$Name[cluster] %in% c("T")] <- "Clonal"
    muts[,"event"] <- paste0(assign, ", ", muts[,"event"])
    
    point.col <- as.character(factor(muts$event, levels=c("Clonal, non-duplicated", "Clonal, duplicated", "Subclonal, present in this sample, non-duplicated", "Subclonal, not present in this sample, non-duplicated", "Uninformative, non-duplicated"), labels=c("blue", "black", "red", "plum3", "yellow")))
    plot(as.numeric(sub(".*_","",rownames(muts))), muts$vaf, pch=20, col=point.col, xlab=paste("Chr", chr, "position"), ylab="Number of mutation copies per cancer cell", las=1,
         main=paste0(pt,k), bty="n", ylim=c(0,2.5))
    legend(2.5E8,2.5, pch=20, col=c("blue", "black", "red", "plum3", "yellow"), legend=c("Clonal, non-duplicated", "Clonal, duplicated", "Subclonal, present in this sample", "Subclonal, absent in this sample", "Uncertain"))
  }
  
  if (!pt %in% "PD47450") {
    point.col <- as.character(factor(muts$event, levels=c("non-duplicated", "Clonal, duplicated"), labels=c("blue", "black")))
    plot(as.numeric(sub(".*_","",rownames(muts))), muts$vaf, pch=20, col=point.col, xlab=paste("Chr", chr, "position"), ylab="Number of mutation copies per cancer cell", las=1,
         main=paste0(pt,k), bty="n", ylim=c(0,2.5))
    legend(2.5E8,2.5, pch=20, col=c("blue", "black"), legend=c("Non-duplicated", "Duplicated"))
    
  }
  return(muts)
}


# The second method estimates the timing directly from the estimated mutation rate and number of mutations acquired before duplication
triploid.direct.est <- function(ploidy.2, bs, ID, dup.region.size.in.Gb, callable.genome.size=5.32) 
{
  # bs is bootstrap output from above
  iter <- bs$R
  mut.rate.bs <- bs$t[,ID] / callable.genome.size
  resamps <- rpois(n=iter, lambda = ploidy.2)
  ests.bs <- resamps / dup.region.size.in.Gb / mut.rate.bs 
  
  return(quantile(ests.bs, c(0.5,0.025,0.975)))
}


# Set up the data-frame for the results
tloc.ids.len <- nrow(nmsr)
timingOut <- data.frame(nmsr, triploid.lwr = rep(0, tloc.ids.len), triploid.est = rep(0, tloc.ids.len), triploid.upr = rep(0, tloc.ids.len), stringsAsFactors = FALSE)
timingOut$Est.mut.rate.per.Gb.per.year <- (timingOut$Pt.fitted.num.muts/timingOut$age) / 6
row.names(timingOut) <- timingOut$Var1

# Run the patients through the functions in turn ([params/timing_inputs.csv](file:///params/timing_inputs.csv));
time.inputs <- read.table( "dpclust/data/timingParams.txt", header=T, stringsAsFactors=F)
pdf("dpclust/plots/timingDotPlots.pdf", width=8, height=10)
par(mfrow=c(7,1), mar=c(4,4,2,15), xpd=T)
for (i in c(1:6,8)) {
  timing(time.inputs[i,])
}
dev.off()

for (i in 1:nrow(time.inputs)) {
  timing(time.inputs[i,]) -> a
  sum(grepl("Clonal, duplicated",a$event))
  timingOut[i, c("triploid.est", "triploid.lwr", "triploid.upr")] <- 
    triploid.direct.est(sum(grepl("Clonal, duplicated",a$event)),
                        boot2, substr(time.inputs[,1],1,7), (time.inputs[i,4]-time.inputs[i,3])/1E9)
}

### as there are no pre-duplication mutations in PD47450 to estimate the timing, we can use the MRCA estimate instead
finalEstimate <- timingOut[!duplicated(substr(rownames(timingOut),1,7)),]
finalEstimate[1,c("triploid.lwr", "triploid.est", "triploid.upr")] <- PI.boot1[,c("lwr","fit","upr")]
finalEstimate$type <- c("MRCA", "1p loss", "1p loss", "1p loss")
finalEstimate$triploid.lwr[finalEstimate$Var1 %in% "PD47453a"] <- NA   ### clearly looking at the VAF data for PD47453 we see that due to low cellularity, the winner's curse is extensive and true non-duplicated mutations cannot be called 

pdf(paste("dpclust/plots/timingSummary.pdf", sep=""), width = 5, height = 5, useDingbats=F)
cols <- c("blue", "red")[factor(finalEstimate$type)]
colsSD <- c("light blue", "pink")[factor(finalEstimate$type)]
plot(NA, bty="n", yaxt="n", xlim=c(0,65), ylab="", xlab="Years", ylim=c(0.4,4.5))
rect(0, (1:4)-0.1, finalEstimate$age, (1:4)+0.1, border=NA, lwd=2, col="grey")
segments(finalEstimate$age,(1:4)-0.1, finalEstimate$age, (1:4)+0.1, lwd=4, lend=1)
mtext(substr(rownames(finalEstimate),1,7), side=2, at=1:4, las=2)
rect(finalEstimate$triploid.lwr, (1:4)-0.1, finalEstimate$triploid.upr, (1:4)+0.1, col=colsSD, lwd=2, border=NA)
segments(finalEstimate$triploid.est, (1:4)-0.1, finalEstimate$triploid.est, (1:4)+0.1, lwd=4, lend=1, col=cols)
legend("topright", c("Age at surgery", "Estimated age of 1p loss", "Estimated age of MRCA"), col=c("black","blue","red"), pch=15)
dev.off()



################# revision blurb to include only single PD47450 sample (mean)
# 
# temp <- nmsr[c(1,6:8),]
# temp[1,"Freq"] <- mean(nmsr[1:5,"Freq"])
# lm(Freq ~ age - 1, temp)
