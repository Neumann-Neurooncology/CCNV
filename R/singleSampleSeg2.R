#' Segments the data using the conumee 2.0 package and visualizes DNA methylation data and generative cumulative plots
#' @param mSetsAnno A list of the RGSet of the target data, the control data and the annotation data
#' @param thresh A positive float (>=0) indicating the threshold for an aberration.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param detail.regions Either NULL or a vector of gene names.
#' @param array_type Type of methylation array used
#' @param showPlot boolean that determines, if the plots will be displayed
#'
#' @return Nothing. Will print the figures to the default plotting terminal.
singleSampleSeg2<- function(mSetsAnno, thresh, colour.amplification, colour.loss, array_type, showPlot){

  #x1 <- conumee2::CNV.bin(conumee2::CNV.fit(query = mSetsAnno$target_mset_loaded, ref = mSetsAnno$control_mset_loaded, anno = mSetsAnno$anno_targets))
  x1 <- conumee2::CNV.bin(conumee2::CNV.fit(query = mSetsAnno$target_mset_loaded, ref = mSetsAnno$control_mset_loaded, mSetsAnno$anno_targets))

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
    allColsToChange <- unique(c(allColsToChange, colRange))
  }
  test <- x1@bin$ratio
  x2 <- as.data.frame(x1@bin$ratio)
  x2_changed <- x2


  fractions <- c(0.2, 0.4, 0.6, 0.8, 1)
  frac <- sample(fractions, 1)

  amplifications <- c(0.5, 1, 1.5, 2)
  strength <- sample(amplifications, 1)

  numColsToChange <- round(ncol(x2) * frac)

  colsToChange <- sample(ncol(x2), numColsToChange)
  x2_changed[allColsToChange, colsToChange] <- x2[allColsToChange, colsToChange] + strength

  x3 <- x1
  for (i in seq_along(x3@bin$ratio)) {
    x3@bin$ratio[[i]][] <- x2_changed[[i]]
  }

  x <- conumee2::CNV.segment(x3)


  conSegData <- dplyr::bind_rows(x@seg$summary, .id = "column_label")
  segmentation_data <- as.data.frame(cbind(conSegData$chrom, conSegData$loc.start, conSegData$loc.end, conSegData$seg.mean, conSegData$ID))

  names(segmentation_data) <- c("chromosome", "start", "end","segmean", "sample")

  segmentation_data$chromosome <- gsub("chr","",segmentation_data$chromosome)
  segmentation_data$chromosome <- as.numeric(segmentation_data$chromosome)
  sample_no <- as.numeric(length(unique(segmentation_data$sample)))
  segmentation_data$start <- as.numeric(segmentation_data$start)
  segmentation_data$end <- as.numeric(segmentation_data$end)
  segmentation_data$segmean <- as.numeric(segmentation_data$segmean)
  segmentation_data <- as.data.frame(segmentation_data)

  overlayPlot <- overlayPlot(mSetsAnno, segmentation_data, colour.amplification, colour.loss, array_type)
  singleFreqPlot <- singleFrequencyPlot(mSetsAnno, segmentation_data, colour.amplification, colour.loss, thresh, array_type)

  if (showPlot == "TRUE"){
    #return all plots and data, even if summaryplot not functioning
    summaryplot <- function(x, threshold = thresh, overlayPlot, singleFreqPlot, segmentationData){
      tryCatch(
        expr = {
          message(CNV.summaryplot(x, threshold = thresh))
        },
        error = function(e){
          message('Summary plot could not be generated.')
          print(e)
        },
        finally = {
          suppressWarnings(print(overlayPlot))
          suppressWarnings(print(singleFreqPlot))
          return(segmentation_data)
        }
      )
    }

    summaryplot(x, threshold = thresh, overlayPlot, singleFreqPlot, segmentationData)
  }

  else {
    return(list(segmentation_data, gene, frac, strength))
  }


}
