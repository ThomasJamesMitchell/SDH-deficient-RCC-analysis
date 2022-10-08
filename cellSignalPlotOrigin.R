############### Code to infer cell of origin from bulk RNA and appropriate single cell reference.
############### Uses output from https://www.nature.com/articles/s41467-021-23949-5


library(ggsci)
library(reshape2)


setwd("C:/Users/tjm/Documents/renal/SDH/RNA/deconvolve")

studies <- c("kidneyNormalAdultAndFoetalAdrenal", "GSE114530", "csfAll", "distalMouse")
names(studies) <- c("Young _et al.", "Hochane _et al.", "Li _et al._ (unpublished)", "Chen _et al.")

pdf("plots/exposurePlot_RefComparison.pdf", width=12, height=4)
par(mfrow=c(1,4))
for (study in studies) {
  exposures = read.delim(paste0('./output/', study,'_fitExposures.tsv'), header=T, row.names = 1)
  
  ### remove these lines as fit/ obs counts do not tally
  exposures <- exposures[!rownames(exposures) %in% c("pR2", "fitCount", "obsCount"),]
  
  exposures <- apply(exposures, 2, function(x) x/sum(x))
  
  ### only include the SDH_renal tumours (remove the ccRCC and the paraganglioma)
  exposures <- exposures[,c(1:6,8,10)]
  colnames(exposures) <- rep(names(studies)[studies %in% study],8)

  include <- rownames(exposures)[apply(exposures, 1, mean)>0.02]
  include <- rownames(exposures)[apply(exposures, 1, max)>0.1]
  print(length(include))
  exposures <- exposures[include,]
  exposurePlot(exposures, groupBySigs = F, labCex = 0.6)
}
dev.off()

##### full plot for all SDH mouse exposures
studies <- "distalMouse"
names(studies) <- "Chen _et al."
pdf("plots/exposurePlot_RefSDHdistalMouseAll.pdf", width=10, height=6)
for (study in studies) {
  exposures = read.delim(paste0('./output/', study,'_fitExposures.tsv'), header=T, row.names = 1)
  
  ### remove these lines as fit/ obs counts do not tally
  exposures <- exposures[!rownames(exposures) %in% c("pR2", "fitCount", "obsCount"),]
  
  exposures <- apply(exposures, 2, function(x) x/sum(x))
  
  ### only include the SDH_renal tumours (remove the ccRCC and the paraganglioma)
  exposures <- exposures[,c(1:6,8,10)]
  colnames(exposures) <- rep(names(studies)[studies %in% study],8)
  exposures <- exposures[c("PT","CTAL???1","CTAL???2","DCT1-a","DCT1-b",
                           "DCT1-c","DCT1-d","DCT1-pro","DCT1-inf","DCT2","MD",
                           "CNT","A-IC","B-IC", "Intercept"),]
  
  # include <- rownames(exposures)[apply(exposures, 1, mean)>0.02]
  # include <- rownames(exposures)[apply(exposures, 1, max)>0.1]
  # print(length(include))
  # exposures <- exposures[include,]
  exposurePlot(exposures, groupBySigs = F, labCex = 0.6)
}
dev.off()

##### full plot for all SDH mouse exposures
studies <- "stewart"
names(studies) <- "Stewart _et al."
pdf("plots/exposurePlot_RefSDHStewart.pdf", width=6, height=5)
par(mar=c(5,4,4,2))
for (study in studies) {
  exposures = read.delim(paste0('./output/', study,'_fitExposures.tsv'), header=T, row.names = 1)
  
  ### remove these lines as fit/ obs counts do not tally
  exposures <- exposures[!rownames(exposures) %in% c("pR2", "fitCount", "obsCount"),]
  
  exposures <- apply(exposures, 2, function(x) x/sum(x))
  
  ### only include the SDH_renal tumours (remove the ccRCC and the paraganglioma)
  exposures <- exposures[,c(1:6,8,10)]
  colnames(exposures) <- rep(names(studies)[studies %in% study],8)
  exposures <- exposures[c("Proximal tubule", "Epithelial progenitor cell",
                           "Thick ascending limb of Loop of Henle","Connecting tubule",
                           "Principal cell","Intercalated cell","Pelvic epithelium",
                           "Transitional urothelium","Intercept"),]
  
  # include <- rownames(exposures)[apply(exposures, 1, mean)>0.02]
  # include <- rownames(exposures)[apply(exposures, 1, max)>0.1]
  # print(length(include))
  # exposures <- exposures[include,]
  exposurePlot(exposures, groupBySigs = F, labCex = 0.6)
}
dev.off()


defColFuns = list(pal_npg(),
                  pal_igv(),
                  pal_futurama(),
                  pal_rickandmorty(),
                  pal_simpsons(),
                  pal_uchicago(),
                  pal_jama()
)


#' Boxplots with points overlaid
#'
#' Groups data either first by tumour type, then signatures or the other way around.  Only boxplot is show for control samples and with decreased width.  Colours are used to indicate the second split.  Points are either shown using jitter or by binning and creating a dot plot.
#'
#' @param dat Matrix with rows representing signatures, columns samples named in the format <TumourType>___<SampleID> and entries the exposures.
#' @param isControl A named vector indicating which of the tumour types are being used as controls.
#' @param groupBySigs Main grouping is by tumour type by default, if set to TRUE, do main grouping by Signatures.
#' @param dotPlot Default behaviour is to jitter points.  This instead produces a dot plot.
#' @param greyScale Don't use colours, use greyscale.
#' @param controlWidth The width, relative to 1 for a normal group, devoted to each control.
#' @param colPal Colour palette to use.  Can also be named list of colours.
#' @param jitterSDs How many standard deviations to use to jitter x.  Higher number means less jitter
#' @param splitLims Region within [0,1] to use for actually plotting in the first split.
#' @param pointLims Region within [0,1] to use for plotting data within the second split.
#' @param nBins Number of vertical bins for dot plot.
#' @param xPtWid How much space per point in x direction.  If NULL set so that the biggest bin in the plot is filled.
#' @param yPtWid How much space per point in y direction.
#' @param fixNumX Have this many points per row, overrides xPtWid.
#' @param fixNumY Have this many points per col, overrides yPtWid.
#' @param scaleToBox Scale number of points so the fullest bin in each box fills that row.
#' @param packDense Err on the side of more points?
#' @param maxY y-limit for plot
#' @param ptAlpha Alpha value for points.  If a value greater than one, alpha set dynamically such that alpha*numPoints = ptAlpha
#' @param pch Shape of points.
#' @param cex Size of points.
#' @param labCex Size of labels.
#' @param boxSqueeze Scale boxplot width relative to dot plotting by this factor.
#' @param yLabel Y-axis label.
#' @param controlColChange Shift all the colour to the boxes.  Tweaks everything to make this work.
#' @param abb Abbreviation for split labels.  If NULL calculated automatically.
exposurePlot = function(dat,isControl=NULL,groupBySigs=FALSE,dotPlot=TRUE,greyScale=FALSE,controlWidth=0.2,colPal = defColFuns[[1]],jitterSDs=5,splitLims=c(0.1,0.9),pointLims=c(0.05,0.95),nBins=100,xPtWid=0.01,yPtWid=0.01,fixNumX=1,fixNumY=1,scaleToBox=TRUE,packDense=TRUE,maxY=ceiling(max(dat)),ptAlpha=0.2,pch=19,cex=0.8,labCex=3.0,boxSqueeze=0.8,yLabel='Exposure',controlColChange=FALSE,abb=NULL){
  ##############
  # Preprocess #
  ##############
  #Stick the controls on the right
#  dat = dat[,order(isControl[gsub('___.*','',colnames(dat))]),drop=FALSE]
  splitWidth = splitLims[2]-splitLims[1]
  df = melt(as.matrix(dat),as.is=TRUE,varnames=c('Signatures','Sample'))
  df$class = gsub('___.*','',df$Sample)
  if(is.null(isControl))
    isControl = setNames(rep(FALSE,length(unique(df$class))),unique(df$class))
  #Designate controls
  df$isControl = isControl[df$class]
  #Define order
  df$Signatures = factor(df$Signatures,levels=rownames(dat))
  df$class = factor(df$class,levels=unique(df$class[match(colnames(dat),df$Sample)]))
  #Check that controls are all together 
  if(length(which(isControl[levels(df$class)]))>1 & any(diff(which(isControl[levels(df$class)]))!=1))
    stop("Controls must form a continuous block.  Please re-order.")
  if(groupBySigs){
    tmp = colnames(df)
    tmp[match(c('Signatures','class'),tmp)]=c('class','Signatures')
    colnames(df) = tmp
    xBreaks = seq(nlevels(df$class)+1)-0.5
    defWidths = rep(splitWidth/nlevels(df$Signatures),nlevels(df$Signatures))
    defWidths = defWidths*ifelse(isControl[levels(df$Signatures)],controlWidth,1)
    defWidths = defWidths/sum(defWidths)*splitWidth
    x2Breaks = c(0.1,0.1+cumsum(defWidths))
  }else{
    #Boundaries of each class
    xBreaks = c(0.5,0.5+cumsum(ifelse(isControl[levels(df$class)],controlWidth,1)))
    #Boundaries of each signature, assuming class area has width 1
    x2Breaks = seq(splitLims[1],splitLims[2],length.out=nlevels(df$Signatures)+1)
  }
  #Boundaries of plot area within signature, assuming signature area has width 1
  x3Breaks = pointLims
  xWidths = diff(xBreaks)
  #y-breakpoints
  yBreaks = seq(floor(min(dat)),maxY,length.out=nBins+1)
  #############################
  # Calculate point locations #
  #############################
  xx = rep(NA,nrow(df))
  yy = rep(NA,nrow(df))
  #Coordinates for extra annotation
  boxVals = list()
  if(is.null(xPtWid)){
    #Determine the number of points in the largest bin that we will plot
    bigBinNo = df[!df$isControl,]
    bigBinNo = with(bigBinNo,split(value,interaction(Signatures,class)))
    bigBinNo = lapply(bigBinNo,cut,yBreaks,include.lowest=TRUE)
    bigBinNo = sapply(lapply(bigBinNo,table),max)
    bigBinNo = max(bigBinNo)
    #Set point spacing so the fullest bin would be just filled
    xPtWid = max(xWidths[!isControl[names(xWidths)]])*max(diff(x2Breaks))/bigBinNo
  }
  for(i in levels(df$class)){
    boxVals[[i]]=list()
    for(j in levels(df$Signatures)){
      #Get the data
      w = which(df$class==i & df$Signatures==j)
      #Define region
      xLeft = xBreaks[match(i,levels(df$class))] + x2Breaks[match(j,levels(df$Signatures))]*xWidths[match(i,levels(df$class))]
      xRight = xBreaks[match(i,levels(df$class))] + x2Breaks[match(j,levels(df$Signatures))+1]*xWidths[match(i,levels(df$class))]
      #Tweak by final level of breaks
      xMid = (xLeft+xRight)*0.5
      xLeft = xMid + (x3Breaks[1]-0.5)*(xRight-xLeft)
      xRight = xMid + (x3Breaks[2]-0.5)*(xRight-xLeft)
      #Get number of points in the x-direction
      if(is.null(fixNumX)){
        nX = (xRight-xLeft)/xPtWid
        if(packDense){
          nX = ceiling(nX)
        }else{
          nX = max(1,floor(nX))
        }
      }else{
        nX = fixNumX
      }
      #Get quantiles
      boxVals[[i]][[j]]=list(quant=quantile(df$value[w]),xLeft=xLeft,xRight=xRight,dat=df$value[w])
      #Give plotting order for each point in each bin
      if(dotPlot){
        #tmp = beeswarm(df$value[w],method='hex',breaks=seq(0,1,length.out=nBins+1),corral='wrap',do.plot=FALSE,cex=cex)
        #yy[w] = tmp$y
        ##Centre around 0
        #tmp = (tmp$x -1)
        ##Boundaries should be at -0.5/0.5, scale them so they're at xLeft/xRight
        #tmp = tmp*(xRight-xLeft)
        ##Shift centre
        #tmp = tmp+0.5*(xLeft+xRight)
        #xx[w] = tmp 
        #Split into bin
        cc = cut(df$value[w],yBreaks,include.lowest=TRUE)
        if(scaleToBox)
          nX = max(table(cc))
        for(k in levels(cc)){
          #Process each y-bin at a time
          ww = which(cc==k)
          o = seq_along(ww)-1
          #Get number of points in the y-direction
          ii = match(k,levels(cc))
          if(is.null(fixNumY)){
            nY = (yBreaks[ii+1]-yBreaks[ii])/yPtWid
            if(packDense){
              nY = ceiling(nY)
            }else{
              nY = max(1,floor(nY))
            }
          }else{
            nY = fixNumY
          }
          #Line stacker
          #nY = 1
          #nX = length(ww)
          #Convert plotting order to positions relative to middle of bin
          #Place in y-direction first
          a = o %% nY
          #Oscillate around middle
          a = ceiling(a/2)*ifelse(a%%2,1,-1)
          #Then in x-direction
          b = (o %/% nY) %% nX
          b = ceiling(b/2)*ifelse(b%%2,1,-1)
          #Define the x point positions for each
          xPos = seq(xLeft,xRight,length.out=nX+1)
          xPos = xPos[-length(xPos)]+diff(xPos)/2
          #And thy y point positions, which depends on the bin number
          yPos = seq(yBreaks[ii],yBreaks[ii+1],length.out=nY+1)
          yPos = yPos[-length(yPos)]+diff(yPos)/2
          #Finally construct the positions
          xPos = xPos[b+ceiling(nX/2)]
          yPos = yPos[a+ceiling(nY/2)]
          #Centre it?
          xPos = xPos - mean(xPos) + (xLeft+xRight)/2
          yPos = yPos - mean(yPos) + (yBreaks[ii]+yBreaks[ii+1])/2
          xx[w[ww]] = xPos
          yy[w[ww]] = yPos
        }
      }else{
        #Jitter instead
        xx[w] = rnorm(length(w),mean=0.5*(xLeft+xRight),sd = (xRight-xLeft)/jitterSDs) 
        yy[w] = df$value[w]
      }
    }
  }
  ###############
  # Do the plot #
  ###############
  #Define colours
  if(greyScale){
    cols = rep('grey',nlevels(df$Signatures))
  }else{
    if(is.function(colPal)){
      cols = colPal(nlevels(df$Signatures))
    }else{
      cols = colPal[levels(df$Signatures)]
    }
    #Repeat if needed
    if(any(is.na(cols))){
      nCols = max(which(!is.na(cols)))
      cols = cols[((seq_along(cols)-1) %% nCols)+1]
    }
  }
  names(cols) = levels(df$Signatures)
  par(mar=c(7.1, 4.1, 4.1, 9.1))
  #Define the empty plot area
  plot(0,0,type='n',xlab='',ylab=yLabel,ylim=c(floor(min(dat)),maxY),xlim=range(xx),xaxt='n',frame.plot=FALSE,las=1,cex.axis=labCex)
  #Dynamically determine alpha
  if(ptAlpha>1){
    dynAlpha = ptAlpha/table(df$Signatures,df$class)
    m = cbind(match(df$Signatures,rownames(dynAlpha)),match(df$class,colnames(dynAlpha)))
    alpha = pmin(1,dynAlpha[m])
  }else{
    alpha = ptAlpha
  }
  w = which(!controlColChange | !df$isControl)
  if(controlColChange){
    ptCols = colAlpha('#606060',alpha)
  }else{
    ptCols = colAlpha(cols[df$Signatures],alpha)[w] 
  }
  points(xx[w],yy[w],
       col=ptCols,
       pch=pch,
       cex=cex,
       xaxt='n',
       yaxt='n'
       )
  #Plot boxplots
  boxDat = unlist(lapply(boxVals,function(e) lapply(e,function(ee) ee$dat)),recursive=FALSE)
  boxPos = unlist(lapply(boxVals,function(e) lapply(e,function(ee) 0.5*(ee$xLeft+ee$xRight))))
  boxWidth = unlist(lapply(boxVals,function(e) lapply(e,function(ee) (ee$xRight-ee$xLeft))))
  #Boxplot
  if(!groupBySigs){
    sigNoms = gsub('^(.+)\\.(.+)$','\\1',names(boxDat))
    classNoms = gsub('^(.+)\\.(.+)$','\\2',names(boxDat))
  }else{
    #This is a bit weird and hacky, but is what is needed to make the colour selection work for boxplots
    sigNoms = gsub('^(.+)\\.(.+)$','\\2',names(boxDat))
    classNoms = gsub('^(.+)\\.(.+)$','\\2',names(boxDat))
  }
  boxplot(boxDat,
          at=boxPos,
          width=boxWidth*boxSqueeze,
          pars=list(boxwex=max(boxWidth)*boxSqueeze,#Make widths absolute, not relative
                    staplewex = 0, #Turn off things at top
                    whisklty=1), #Solid whiskers
          outline=FALSE,
          boxlwd = ifelse(controlColChange,3,1),
          whisklwd = ifelse(controlColChange,3,1),
          #border = ifelse(controlColChange & isControl[sigNoms],
          #                cols[classNoms],
          #                'black'
          #                ),
          border = cols[classNoms],
          xaxt='n',
          yaxt='n',
          col =colAlpha('white',0.4),
          frame.plot=FALSE,
          add=TRUE
          )
  #Define an abbreviation for the inner split.  Increase string length until unique
  if(is.null(abb)){
    #Increase each name independently
    ll=rep(NA,nlevels(df$Signatures))
    for(i in seq_along(ll)){
      j=1
      while(j<=max(nchar(levels(df$Signatures)))){
        abb = substr(levels(df$Signatures),1,j)
        if(sum(abb==abb[i])==1)
          break
        j = j +1
      }
      ll[i]=j
    }
    abb = setNames(substr(levels(df$Signatures),1,ll),levels(df$Signatures))
  }
  #Create a series of axes labelling things
  for(i in seq_along(xBreaks[-1])){
    off = (x2Breaks[-length(x2Breaks)]+diff(x2Breaks)/2)*(xBreaks[i+1]-xBreaks[i])
    #Create the deeper split labels
    axis(1,at=xBreaks[i]+off,labels = abb,cex.axis=labCex)
    #Now add the larger axis
    axis(1,at=c(xBreaks[i]+off[1],0.5*(xBreaks[i]+xBreaks[i+1]),xBreaks[i]+off[length(off)]),labels = c(NA,levels(df$class)[i],NA),line=2.5,lwd.ticks=0,cex.axis=labCex)
  }
  #Finally add control label if needed
  if(!groupBySigs){
    tmp = isControl[levels(df$class)]
    #Do we need to draw it?
    if(any(tmp)){
      low = xBreaks[min(which(tmp))]
      high = xBreaks[max(which(tmp))+1]
      axis(1,at=c(low,0.5*(low+high),high),labels=c(NA,'Controls',NA),line=5,lwd.ticks=0,cex.axis=labCex)
    }
  }
  #Add axis label
  #title(xlab=ifelse(groupBySigs,'Signatures','TumourType'),line=5)
  #Add legends
  if(!groupBySigs){
    #Add legend
    legend(x=max(xx)+(max(xx)-min(xx))*0.05,
           y=max(yy),
           legend=levels(df$Signatures),
           col = cols,
           pch = pch,
           bty='n',
           title='Signatures',
           cex=cex,
           xpd=TRUE)
  }else{
    #Add extra heading to legend to show controls
    labs = levels(df$Signatures)
    cc = cols
    #Order so controls at the bottom
    o = order(isControl[labs])
    labs = labs[o]
    cc = cc[o]
    #Add in extra label
    cc = c(cc[seq(sum(!isControl))],NA,cc[-seq(sum(!isControl))])
    labs = c(labs[seq(sum(!isControl))],'Controls',labs[-seq(sum(!isControl))])
    #Add legend
    legend(x=max(xx)+(max(xx)-min(xx))*0.05,
           y=max(yy),
           legend=labs,
           col = cc,
           pch = pch,
           bty='n',
           title='TumourTypes',
           cex=cex,
           xpd=TRUE)
  }
}


#' Adds transparency to colour
#'
#' @param cols Vector of colours.
#' @param alphas Single value or vector of alphas
#' @param ... Passed to rgb
#' @return rgb colours with transparency set.
colAlpha = function(cols,alphas,...) {
  if(length(alphas)==1)
    alphas = rep(alphas,length(cols))
  tmp = col2rgb(cols)
  sapply(seq_len(ncol(tmp)),function(e) rgb(tmp[1,e],tmp[2,e],tmp[3,e],alphas[e]*255,maxColorValue=255,...))
}

