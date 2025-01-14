#'Determines the array type of the user input dataframe
#'
#'@param dataFiles Dataframe with a column labelled "ArrayType" containing the arraytype (450k/EPIC/EPIC2) for each sample.
#'@return A string with either "450k" (only 450k samples), "EPIC" (only EPIC samples), "combined" (both 450k and EPIC samples),  EPICv2 (only EPICv2)
# overlap.1 (450k and EPIC), overlap.2  (EPIC and EPICv2) and overlap.3 (450k, EPIC, EPICv2).

get.ArrayType <- function(dataFiles) {
    ArrayType <- NULL
    types = unique(dataFiles$ArrayType)
    if ((all(c("EPIC", "450k") %in% types)) &
        (length(types) == 2)) {
        ArrayType <- "overlap.1"
    } else if ((all(c("EPIC", "EPICv2") %in% types)) &
               (length(types) == 2)) {
        ArrayType <- "overlap.2"
    } else if ((all(c("EPIC", "EPICv2", "450k") %in% types)) &
              (length(types) == 3)) {
        ArrayType <- "overlap.3"
    } else if (("450k" %in% types) & (length(types) == 1)) {
        ArrayType <- "450k"
    } else if (("EPIC" %in% types) & (length(types) == 1)) {
        ArrayType <- "EPIC"
    } else if (("EPICv2" %in% types) & (length(types) == 1)) {
        ArrayType <- "EPICv2"
    }else if (("mouse" %in% types) & (length(types) == 1)) {
      ArrayType <- "mouse"
    }
    return(ArrayType)
}

#' Determines the conumee version from the array type
#'
#' @param ArrayType The ArrayType as returned by get.ArrayType
#'
#' @return The conumee version (either 1 or 2)
get.ConumeeVersion <- function(ArrayType) {
    v <- 2
    if (ArrayType %in% c("450k", "EPIC")) {
        v <- 1
    }
    return(v)
}

#' Reads the specified methylation arrays into an RGSet. Note that any combination of 450k and EPIC arrays will be coerced into an RGSet of 450k type.
#' @param dataFiles Dataframe with a column batch name as requested by minfi for reading in experiments.
#' @param ArrayTye A string (either "450k", "EPIC", "combined" or "EPIC2")
#'
#' @return A list of the RGSet of the target data, the control data and the annotation data
read.RGSet <- function(dataFiles, ArrayType) {
    types = unique(dataFiles$ArrayType)
    
    #read data transform to RGChannelSet
    if (ArrayType == "450k") {
        rgset_450k <-
            minfi::read.metharray.exp(targets = dataFiles, force = TRUE)
        target_rgset <- rgset_450k
    }
    if (ArrayType == "EPIC") {
        rgset_EPIC <-
            minfi::read.metharray.exp(targets = dataFiles, force = TRUE)
        target_rgset <- rgset_EPIC
        
    }
    if (ArrayType %in% c("EPICv2", "mouse", "overlap.1", "overlap.2", "overlap.3")) {
      require(sesame)
      require(sesameData)
      
      if ( ArrayType == "overlap.1") {
        #separate into 450k and EPIC samples
        data_EPIC <-
          dataFiles[which(dataFiles$ArrayType == "EPIC"),]
        data_450k <-
          dataFiles[which(dataFiles$ArrayType == "450k"),]
        #read with sesame
        if (length(data_450k$Basename) == 1) {
          samplename <- basename(data_450k$Basename)
          sdf.450 <- list(name = sesame::openSesame(data_450k$Basename, prep = "QCDPB", func = NULL))
          names(sdf.450) <- samplename
        } else {
        sdf.450 <-
          sesame::openSesame(data_450k$Basename, prep = "QCDPB", func = NULL)
        }
        
        if (length(data_EPIC$Basename) == 1) {
          samplename <- basename(data_EPIC$Basename)
          sdf.EPIC <- list(name = sesame::openSesame(data_EPIC$Basename, prep = "QCDPB", func = NULL))
          names(sdf.EPIC) <- samplename
        } else {
          sdf.EPIC <-
            sesame::openSesame(data_EPIC$Basename, prep = "QCDPB", func = NULL)
        }
        
       #liftOver to smaller platform
        
        sdf.lift.EPIC.450 <- sesame::mLiftOver(sdf.EPIC, "HM450", impute = FALSE)
        
        sdf.all <- c(sdf.450, sdf.lift.EPIC.450)
      } else if (ArrayType == "overlap.2") {
        #separate into 450k and EPIC samples
        data_EPIC <-
          dataFiles[which(dataFiles$ArrayType == "EPIC"),]
        data_EPICv2 <-
          dataFiles[which(dataFiles$ArrayType == "EPICv2"),]
        #read with sesame
        if (length(data_EPIC$Basename) == 1) {
          samplename <- basename(data_EPIC$Basename)
          sdf.EPIC <- list(name = sesame::openSesame(data_EPIC$Basename, prep = "QCDPB", func = NULL))
          names(sdf.EPIC) <- samplename
        } else {
          sdf.EPIC <-
            sesame::openSesame(sdf.EPIC$Basename, prep = "QCDPB", func = NULL)
        }
        if (length(data_EPICv2$Basename) == 1) {
          samplename <- basename(data_EPICv2$Basename)
          sdf.EPICv2 <- list(name = sesame::openSesame(data_EPICv2$Basename, prep = "QCDPB", func = NULL))
          names(sdf.EPICv2) <- samplename
        } else {
          sdf.EPICv2 <-
            sesame::openSesame(data_EPICv2$Basename, prep = "QCDPB", func = NULL)
        }
        
        #liftOver to smaller platform
        
        sdf.lift.EPICv2.EPIC <- sesame::mLiftOver(sdf.EPICv2, "EPIC", impute = FALSE)
        
        sdf.all <- c(sdf.EPIC, sdf.lift.EPICv2.EPIC)
      } else if (ArrayType == "overlap.3") {
        #separate into 450k and EPIC samples
        data_450k <-
          dataFiles[which(dataFiles$ArrayType == "450k"),]
        data_EPIC <-
          dataFiles[which(dataFiles$ArrayType == "EPIC"),]
        data_EPICv2 <-
          dataFiles[which(dataFiles$ArrayType == "EPICv2"),]
        #read with sesame
        if (length(data_450k$Basename) == 1) {
          samplename <- basename(data_450k$Basename)
          sdf.450 <- list(name = sesame::openSesame(data_450k$Basename, prep = "QCDPB", func = NULL))
          names(sdf.450) <- samplename
        } else {
          sdf.450 <-
            sesame::openSesame(data_450k$Basename, prep = "QCDPB", func = NULL)
        }
        if (length(data_EPIC$Basename) == 1) {
          samplename <- basename(data_EPIC$Basename)
          sdf.EPIC <- list(name = sesame::openSesame(data_EPIC$Basename, prep = "QCDPB", func = NULL))
          names(sdf.EPIC) <- samplename
        } else {
          sdf.EPIC <-
            sesame::openSesame(data_EPIC$Basename, prep = "QCDPB", func = NULL)
        }
        if (length(data_EPICv2$Basename) == 1) {
          samplename <- basename(data_EPICv2$Basename)
          sdf.EPICv2 <- list(name = sesame::openSesame(data_EPICv2$Basename, prep = "QCDPB", func = NULL))
          names(sdf.EPICv2) <- samplename
        } else {
          sdf.EPICv2 <-
            sesame::openSesame(data_EPICv2$Basename, prep = "QCDPB", func = NULL)
        }
        #liftOver to smaller platform
        sdf.lift.EPIC.450 <- sesame::mLiftOver(sdf.EPIC, "HM450", impute = FALSE)
        sdf.lift.EPICv2.450 <- sesame::mLiftOver(sdf.EPICv2, "HM450", impute = FALSE)
        
        sdf.all <- c(sdf.450, sdf.lift.EPIC.450, sdf.lift.EPICv2.450)
      } else if (ArrayType == "mouse") {
        #read with sesame
        sdf.mouse <-
          sesame::openSesame(dataFiles$Basename, prep = "QCDPB", func = NULL)
        
        sdf.all <- sdf.mouse
      } else {
        #read with sesame
        sdf.EPICv2 <-
          sesame::openSesame(dataFiles$Basename, prep = "QCDPB", func = NULL)
        
        sdf.all <- sdf.EPICv2
      }

      target_rgset <- list("array_type" = ArrayType, "sdfs" = sdf.all)
    } 
    
    return(target_rgset)
}

#' Segments and visualizes DNA methylation data and generative cumulative plots
#' @param dataFiles A dataframe with columns ArrayType and Basename
#' @param segmentationMode specifying the segmentation mode. Allowed values are single/multi/all.
#' @param thresh A positive float (>=0) indicating the threshold for an abberation.
#' @param gamma A positive integer (>0) indicating the threshold for the multi-sample segmentation.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param detail.regions Either NULL or a vector of gene names.
#' @param conumee.version The conumee version to use.
#' @param controls control data to normalize against, if the publicly available should not be used or ArrayType = EPICv2 or mouse
#' @param showPlot boolean that determines, if the plots will be displayed
#'
#' @return Returns the segmentation values either as a dataframe or a list of two dataframes.
segment.Plot <-
    function(target_rgset,
             array_type,
             segmentationMode,
             gamma,
             thresh,
             colour.amplification,
             colour.loss,
             detail.regions,
             conumee.version,
             controls,
             showPlot){
      if(conumee.version == 1) {
        require(conumee)
        mSetsAnno <-  sampleBinContr(target_rgset, array_type, controls)
      } else {
        require(conumee2)
        mSetsAnno <-  sampleBinContr2(target_rgset, array_type, controls)
      }
      require(ggplot2)
        if (segmentationMode == "single") {
            if (conumee.version == 1) {
                singleSeg <- singleSampleSeg(mSetsAnno,
                                thresh,
                                colour.amplification,
                                colour.loss,
                                array_type,
                                showPlot)
            } else{
              singleSeg <- singleSampleSeg2(mSetsAnno,
                                 thresh,
                                 colour.amplification,
                                 colour.loss,
                                 array_type,
                                 showPlot)
            }
          return(singleSeg)
        } else if (segmentationMode == "multi") {
          if (conumee.version == 1) {
            multiSeg <- multiSampleSeg(mSetsAnno,
                           thresh,
                           array_type,
                           colour.amplification,
                           colour.loss,
                           detail.regions,
                           showPlot,
                             gamma)
          } else{
            multiSeg <- multiSampleSeg2(mSetsAnno,
                            thresh,
                            array_type,
                            colour.amplification,
                            colour.loss,
                            detail.regions,
                            showPlot,
                            gamma)
          }
          return(multiSeg)
        } else if (segmentationMode == "all") {
            if (conumee.version == 1) {
              
              singleSeg <- singleSampleSeg(mSetsAnno,
                                thresh,
                                colour.amplification,
                                colour.loss,
                                array_type,
                                showPlot)
              
              multiSeg <- multiSampleSeg(mSetsAnno,
                             thresh,
                             array_type,
                             colour.amplification,
                             colour.loss,
                             detail.regions,
                             showPlot,
                             gamma)
              
            } else{
              
              singleSeg <- singleSampleSeg2(mSetsAnno,
                                 thresh,
                                 colour.amplification,
                                 colour.loss,
                                 array_type,
                                 showPlot)
              
              multiSeg <- multiSampleSeg2(mSetsAnno,
                             thresh,
                             array_type,
                             colour.amplification,
                             colour.loss,
                             detail.regions,
                             showPlot,
                             gamma)
              
            }
          
            output <-
              list(
                "multiSeg" = multiSeg,
                "singleSeg" = singleSeg
              )
            return(output)
        }
    }

#' compute the cumulative CNV plots
#' @param dataFiles A dataframe with columns ArrayType and Basename
#' @param segmentationMode specifying the segmentation mode. Allowed values are single/multi/all.
#' @param thresh A positive float (>=0) indicating the threshold for an abberation.
#' @param gamma A positive integer (>0) indicating the threshold for the multi-sample segmentation.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param detail.regions Either NULL or a vector of gene names.
#' #conumee version new name
#' @param conumee.version The version of conumee to use (either 1 or 2). 1 is incompatible with mouse or EPICv2 arrays. NULL will set the version heuristically to 1 for 450k, EPIC and to 2 for Mouse and EPICv2
#' @param output determines the type of output. Can be either plot, data or all
#' @param controls control samples as a dataframe with columns ArrayType and Basename or NULL. If NULL, publicly available can be used (except for mouse or EPICv2 or mouse data)
#'
#' @return If the output is set to data or all, te segmentation values will be returned, else nothing will be returned and the figures are printed to the default plotting terminal.
#' @export
cumul.CNV <-
    function(dataFiles,
             segmentationMode = "all",
             thresh = 0.2,
             gamma = 5,
             colour.amplification = "red3",
             colour.loss = "blue4",
             detail.regions = NULL,
             conumee.version = NULL,
             output = "plot",
             controls = NULL) {
        # check user input
        stopifnot("dataFiles must be a dataframe" = typeof(dataFiles) == "list")
        stopifnot("Basename must be a column of input dataframe dataFiles" =
                      ("Basename" %in% colnames(dataFiles)))
        stopifnot(
            "ArrayType must be a column of input dataframe dataFiles" = ("ArrayType" %in% colnames(dataFiles))
        )
        stopifnot(
            "Parameter segmentationMode must be one of single/multi/all" = segmentationMode %in% c("multi", "all", "single")
        )
        stopifnot("Parameter thresh must be a float >=0" = (thresh >= 0) &&
                      (typeof(thresh) == "double"))
        stopifnot("Parameter gamma must be an integer >0" = (gamma > 0) &&
                      (typeof(gamma) == "double"))
        stopifnot(
            "Parameter colour.amplification must be a string" = typeof(colour.amplification) ==
                "character"
        )
        stopifnot("Parameter colour.loss must be a string" = typeof(colour.loss) ==
                      "character")
        stopifnot(
            "Parameter colour.loss must be a string" = (detail.regions == NULL) &
                (typeof(detail.regions) == "list")
        )
        stopifnot(
            "Parameter detail.regions must either be NULL or a vector of strings" =
                (
                    is.null(detail.regions) | typeof(detail.regions) == "character"
                )
        )
        stopifnot(
            "Parameter conumee.version must either be NULL, 1 or 2" = (conumee.version %in% c(1, 2) ||
                                                                           is.null(conumee.version))
        )
        stopifnot(
          "Parameter output must either be plot, data or all" = (output %in% c("plot", "data", "all"))
        )
        
        # determine type of input files
        array_type <- get.ArrayType(dataFiles)
        # determine the conumee version (if not specified by the user)
        if (is.null(conumee.version)) {
            conumee.version <- get.ConumeeVersion(array_type)
        }
        # determine whether to return the segmentation values
        if (output == "plot") {
          segVal = FALSE
        } else {
          segVal = TRUE
        }
        
        if (output == "data") {
        showPlot = "FALSE"
        } else {
          showPlot = "TRUE"
        }
        
        if(!is.null(controls)) {
          controlsRG <- read.RGSet(controls, array_type)
          if (array_type %in% c("EPICv2", "mouse")) {
            require(conumee2)
            control_mset <- conumee2::CNV.import(controlsRG$array_type, controlsRG$directory, controlsRG$sample_sheet)
          }
          else {
            control_mset <- minfi::preprocessIllumina(controlsRG)
          }
          controls <- control_mset
        }
        
        # read in RGSet
        target_rgset <- read.RGSet(dataFiles, array_type)
        
        # segment and plot
        SegementationValues <- segment.Plot(
            target_rgset,
            array_type,
            segmentationMode,
            gamma,
            thresh,
            colour.amplification,
            colour.loss,
            detail.regions,
            conumee.version,
            controls,
            showPlot)
        if (segVal == TRUE) {
          return(SegementationValues)
        } else {
          return(NULL)
        }
        
    }
