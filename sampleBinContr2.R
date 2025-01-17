#' binning target data and publicly available control data and mapping them to
#' the genome and generating annotation file for bins using conumee2
#'
#' @param target_rgset RGchannelSet of data to be binned
#' @param ArrayType character representing a that Array Type used for acquiring
#' data
#' @return A list containing three objects, the mapped MethylSet of the target
#' data, the mapped MethylSet of the control Data and the annotation file of
#' the bins.
sampleBinContr2 <- function(target_rgset, ArrayType, controls) {
  #generate bins with some good default values
  
  
  if (ArrayType == "mouse") {
    anno_targets <- conumee2::CNV.create_anno(array_type = ArrayType, chrXY = TRUE)
    target_mset_loaded <- conumee2::CNV.load(do.call(cbind, lapply(target_rgset$sdfs, totalIntensities)))
  }  else if (ArrayType == "EPICv2") {
    anno_targets <- conumee2::CNV.create_anno(array_type = ArrayType)
    target_mset_loaded <- conumee2::CNV.load(do.call(cbind, lapply(target_rgset$sdfs, totalIntensities)))
  } else if (ArrayType == "overlap.1") {
    anno_targets <- conumee2::CNV.create_anno(array_type = c("450k", "EPIC"))
    target_mset_loaded <- conumee2::CNV.load(do.call(cbind, lapply(target_rgset$sdfs, totalIntensities)))
  } else if (ArrayType == "overlap.2") {
    anno_targets <- conumee2::CNV.create_anno(array_type = c("EPIC", "EPICv2"))
    target_mset_loaded <- conumee2::CNV.load(do.call(cbind, lapply(target_rgset$sdfs, totalIntensities)))
  } else if (ArrayType == "overlap.3") {
    anno_targets <- conumee2::CNV.create_anno(array_type = c("450k", "EPIC", "EPICv2"))
    target_mset_loaded <- conumee2::CNV.load(do.call(cbind, lapply(target_rgset$sdfs, totalIntensities)))
  } else {
      anno_targets <- conumee2::CNV.create_anno(array_type = ArrayType)
      # Illumina normalisation
      target_mset <- minfi::preprocessIllumina(target_rgset)
      target_mset_loaded <- conumee2::CNV.load(target_mset) 
    }
  if (is.null(controls))  {
    #load controls based on ArrayType
    if (ArrayType == "overlap.1" || ArrayType == "450k" || ArrayType == "overlap.3") {
      control_mset <- minfiData::MsetEx
      control_mset_loaded <- conumee2::CNV.load(control_mset)
    }
    if (ArrayType == "EPIC" || ArrayType == "overlap.2") {
      control_mset <- minfiDataEPIC::MsetEPIC
      control_mset_loaded <- conumee2::CNV.load(control_mset)
    }
  } else {
    control_mset_loaded <- controls
  }
  
  
  output <-
    list(
      "target_mset_loaded" = target_mset_loaded,
      "control_mset_loaded" = control_mset_loaded,
      "anno_targets" = anno_targets
    )
  return(output)
}