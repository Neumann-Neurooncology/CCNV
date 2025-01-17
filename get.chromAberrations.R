#' Takes the weighted mean intensity per common segment for each sample or chromosome and determines whether a chromosome 
#' has am amplification (1), a loss (-1) or no change (0). This is loaded into a dataframe and returned. A similarity clustering 
#' is performed using the euclidean distance between the samples to reorder them, based on their similarity. Sample that show similar
#' segment boundaries (from the combined segmentation) can be processed segment wise (chrom.arm = FALSE) or chromosomal arm wise (chrom.arm = TRUE).
#' Samples segmented separately will be processed on the choromosomal level.
#' @param segmentation_data dataframe that has been returned from the function cumul.CNV
#' @param thresh a threshold to determine from when on a segment is considered an aberration
#' @param chrom boolean, whether the aberration should be calculated on the segment level (FALSE) or on the chromosomal level (TRUE)
#'
#' @return returns the frequency plot of the samples that were segmented togethher using the multi sample segmentation algorithm
#' @export
#'
get.chromAberrations <- function(segmentation_data, thresh = 0.2, chrom.arm = FALSE){
  require(dplyr)
  if(all(c("chrom", "start.pos", "end.pos", "n.probes") %in% names(segmentation_data)))
  {
    if (chrom.arm == TRUE) {
      anno <- conumee::CNV.create_anno(array_type = "EPIC")
      annopq <- anno@genome$pq
      
      pq_split_data <- data.frame()
      
      for (i in 1:length(unique(segmentation_data$chrom))){
        
        df_new <- segmentation_data[which(segmentation_data$chrom == i),]
        pq_split <- annopq[i]
        
        for (j in 1:nrow(df_new)) {
          row <- df_new[j, ]
          if (row$start.pos <= pq_split && row$end.pos >= pq_split) {
            new_row1 <- row
            new_row2 <- row
            
            new_row1$end.pos <- pq_split
            new_row1$chrom <- paste0(new_row1$chrom, "p")
            
            new_row2$start.pos <- pq_split + 1
            new_row2$chrom <- paste0(row$chrom, "q")
            
            # Add the two new rows to the pq_split_data
            pq_split_data <- rbind(pq_split_data, new_row1, new_row2)
          } else if (row$end.pos <= pq_split) {
            row$chrom <- paste0(row$chrom, "p")
            pq_split_data <- rbind(pq_split_data, row)
          } else {
            row$chrom <- paste0(row$chrom, "q")
            pq_split_data <- rbind(pq_split_data, row)
          }
        
        }
      }
      
      segmentation_data <- pq_split_data
      
      segmentation_data$n.bp <- segmentation_data$end.pos - segmentation_data$start.pos
      
      segmentation_data_filtered <- segmentation_data[,-which(names(segmentation_data) %in% c("start.pos", "end.pos", "n.probes"))]
      sample_names <- names(segmentation_data_filtered)
      sample_names <- sample_names[-1]
      chromosomal_weights <- segmentation_data_filtered %>% group_by(chrom) %>%  mutate(weight = n.bp / sum(n.bp),) 
      
      chromosomal_average_multi <- chromosomal_weights %>% group_by(chrom) %>% summarise(across(sample_names, ~weighted.mean(., w = weight)))
      chromNames <- chromosomal_average_multi$chrom
      chromosomal_average_multi <- chromosomal_average_multi[,-which(names(chromosomal_average_multi) %in% c("n.bp", "chrom")) ,drop=FALSE]

      chromosomal_average_multi[chromosomal_average_multi > 0 & chromosomal_average_multi <= thresh ] <- 1
      chromosomal_average_multi[chromosomal_average_multi <= -thresh ] <- -1
      chromosomal_average_multi[chromosomal_average_multi != 1 & chromosomal_average_multi != -1 ] <- 0
      
      
      euclideanDist <- hclust(dist(t(chromosomal_average_multi)))
      ordereuclideanDist <- euclideanDist$order
      chromosomal_average_multi <- chromosomal_average_multi[,ordereuclideanDist]
      output <- cbind(chromosome = chromNames, chromosomal_average_multi)
      
    } else {
      segmentation_data_filtered <- segmentation_data[,-which(names(segmentation_data) %in% c("chrom", "start.pos", "end.pos", "n.probes"))]
      segment_average_multi <- segmentation_data_filtered
      segment_average_multi[segment_average_multi > 0 & segment_average_multi <= thresh ] <- 1
      segment_average_multi[segment_average_multi <= -thresh ] <- -1
      segment_average_multi[segment_average_multi != 1 & segment_average_multi != -1 ] <- 0
      
      euclideanDist <- hclust(dist(t(segment_average_multi)))
      ordereuclideanDist <- euclideanDist$order
      segment_average_multi <- segment_average_multi[,ordereuclideanDist]
      output <- as.data.frame(cbind(segmentation_data[,c(1:3)], segment_average_multi))
    }
  }
  else if (all(c("chromosome", "start", "end", "segmean", "sample", "noise") %in% names(segmentation_data)))
    {
    anno <- conumee::CNV.create_anno(array_type = "EPIC")
    annopq <- anno@genome$pq
    
    pq_split_data <- data.frame()
    
    samples <- unique(segmentation_data$sample)
    
      
      for (x in 1: length(unique(segmentation_data$chromosome))){
        
        df_new <- segmentation_data[which(segmentation_data$chromosome == x),]
        pq_split <- annopq[x]
        
        for (j in 1:nrow(df_new)) {
          row <- df_new[j, ]
          if (row$start <= pq_split && row$end >= pq_split) {
            new_row1 <- row
            new_row2 <- row
            
            new_row1$end <- pq_split
            new_row1$chromosome <- paste0(new_row1$chromosome, "p")
            
            new_row2$start <- pq_split + 1
            new_row2$chromosome <- paste0(row$chromosome, "q")
            
            # Add the two new rows to the pq_split_data
            pq_split_data <- rbind(pq_split_data, new_row1, new_row2)
          } else if (row$end <= pq_split) {
            row$chromosome <- paste0(row$chromosome, "p")
            pq_split_data <- rbind(pq_split_data, row)
          } else {
            row$chromosome <- paste0(row$chromosome, "q")
            pq_split_data <- rbind(pq_split_data, row)
          }
        }
    }
    
    segmentation_data <- pq_split_data
    
    
    segmentation_data$n.bp <- segmentation_data$end - segmentation_data$start
    single_data_filtered <- segmentation_data[,-which(names(segmentation_data) %in% c("start", "end", "noise"))]
    
    single_data_weights <- single_data_filtered %>%   group_by(sample, chromosome) %>%   mutate(weight = n.bp / sum(n.bp)) 
    chromosomal_average_single <- single_data_weights %>% group_by(sample, chromosome) %>% summarise(segmean_weighted_mean = weighted.mean(segmean, weight))
    
    require(tidyverse)
    
    
    chromosomal_wide <- chromosomal_average_single %>% pivot_wider(names_from = sample, values_from = segmean_weighted_mean)
    chromNames <- chromosomal_wide[,1]
    chromosomal_average_single <- chromosomal_wide[,-1]
    
    chromosomal_average_single[chromosomal_average_single > 0 & chromosomal_average_single <= thresh ] <- 1
    chromosomal_average_single[chromosomal_average_single <= -thresh ] <- -1
    chromosomal_average_single[chromosomal_average_single != 1 & chromosomal_average_single != -1 ] <- 0
    
    euclideanDist <- hclust(dist(t(chromosomal_average_single)))
    ordereuclideanDist <- euclideanDist$order
    chromosomal_average_single <- chromosomal_average_single[,ordereuclideanDist]
    output <- cbind(chromosome = chromNames, chromosomal_average_single)
  } else {
    stop("Please verify that your input is a dataframe from the cumul.CNV function or is ordered accordingly.")
  }
  return(output)
}