#' Segments the data using the multipcf algorithm and visualizes DNA methylation data and generative cumulative plots using conumee 2.0
#' @param mSetsAnno A list of the RGSet of the target data, the control data and the annotation data
#' @param thresh A positive float (>=0) indicating the threshold for an abberation.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param detail.regions Either NULL or a vector of gene names.
#' @param array_type Type of methylation array
#' @param showPlot boolean that determines, if the plots will be displayed
#'
#' @return Nothing. Will print the figures to the default plotting terminal.
multiSampleSeg2 <- function(mSetsAnno, thresh, array_type, colour.amplification, colour.loss, detail.regions, showPlot, gamma){

  #load and bin each sample in conumee
  x <- conumee2.0::CNV.bin(conumee2.0::CNV.fit(query = mSetsAnno$target_mset_loaded, ref = mSetsAnno$control_mset_loaded, mSetsAnno$anno_targets))

  genes <- mSetsAnno$anno_targets@probes@elementMetadata@listData[["genes"]]
  genes <- na.omit(genes)
  gene <- sample(genes, 1)

  binInfo <- mSetsAnno$anno_targets@bins@elementMetadata@listData[["genes"]]

  bins <- binInfo[grep(gene, binInfo)]

  toChange <- binInfo[which(binInfo %in% bins)]
  toChange2 <- which(binInfo %in% bins)
  allColsToChange <- integer()

  for (colIndex in toChange2) {
    colRange <- (colIndex - 3):(colIndex + 3)
    colRange <- colRange[colRange > 0 & colRange <= ncol(binInfo)]
    allColsToChange <- unique(c(allColsToChange, colRange))
  }

  x1 <- as.data.frame(x@bin$ratio)
  x1_changed <- x1


  fractions <- c(0.2, 0.4, 0.6, 0.8, 1)
  frac <- sample(fractions, 1)

  amplifications <- c(0.5, 1, 1.5, 2)
  strength <- sample(amplifications, 1)

  numColsToChange <- round(ncol(x1) * frac)

  colsToChange <- sample(ncol(x1), numColsToChange)
  x1_changed[allColsToChange, colsToChange] <- x1[allColsToChange, colsToChange] + strength


  target_ratios <- cbind(x@anno@bins@ranges@NAMES, x@anno@bins@ranges@start, x1_changed)

  names(target_ratios) <- c("Chrom", "Median.bp", names(mSetsAnno$target_mset_loaded@intensity))

  target_ratios$Chrom <- gsub("chr","",target_ratios$Chrom)
  target_ratios$Chrom <- substr(target_ratios$Chrom,1,nchar(target_ratios$Chrom)-4)
  target_ratios$Chrom <- sub("\\-", "", target_ratios$Chrom)

  target_ratios <- na.omit(as.data.frame(sapply(target_ratios, as.numeric)))

  #################### Segmentation #################################
  start <- Sys.time()
  seg_mpcf <- FastMultiPCF(target_ratios, gamma = gamma)
  end <- Sys.time()
  execution_time_multi <- as.numeric(as.POSIXct(start,origin = "1970-01-01")) - as.numeric(as.POSIXct(end,origin = "1970-01-01"))
  print(paste("Runtime multi:", execution_time_multi))

  require(dplyr)
  cumCNV <- CCNV(mSetsAnno, seg_mpcf, target_ratios, array_type, colour.amplification, colour.loss, detail.regions, array_type)
  cumFreq <- cumFreq(mSetsAnno, seg_mpcf, target_ratios, colour.amplification, colour.loss, thresh, array_type)

  if (showPlot == "TRUE") {
    #draw plots
    suppressWarnings(print(cumCNV))
    suppressWarnings(print(cumFreq))
  }

  return(list(seg_mpcf, gene, frac, strength))

}
