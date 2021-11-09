library(TFregulomeR)

# tf's in question
tfs <- c("MM1_HSA_H1-hESC_USF1", "MM1_HSA_H1-hESC_USF2", "MM1_HSA_H1-hESC_ATF3")

exclusivePeak_result_output <- c()
commonPeak_result_output <- c()
exclusive_names <- c()
common_names <- c()

# a loop to generate the exclusive peaks and then common peaks
for (target_tf in tfs) {
  # target_tf is the tf in the x angle
  # remaining_tfs is a list of tfs without the current target_tf
  remaining_tfs <- tfs[!(tfs %in% target_tf)]
  
  # generates the peaks for your target_tf that do not include the peaks of the specified tfs 
  exclusivePeak_output <- exclusivePeaks(target_peak_id = target_tf,
                                         motif_only_for_target_peak = TRUE,
                                         excluded_peak_id = remaining_tfs,
                                         TFregulome_url = 'http://localhost/MethMotif/2021/')
  
  exclusivePeak_result <- exclusivePeakResult(exclusivePeaks = exclusivePeak_output,
                                              return_exclusive_peak_sites = TRUE)
  
  # generates a name for the list below
  exclusive_names <- c(exclusive_names, paste0("exclusivePeaks_", target_tf))
  # a list variable to store the results from this loop
  exclusivePeak_result_output <- c(exclusivePeak_result_output, exclusivePeak_result)
  
  # a loop to generate the common peaks for all possible combinations
  for (tf_y in remaining_tfs) {
    # common peaks function returns all shared peaks between target_tf (angle x) and tf_y (angle y)
    commonPeak_output <- commonPeaks(target_peak_id = target_tf,
                                     compared_peak_id = tf_y,
                                     TFregulome_url = 'http://localhost/MethMotif/2021/')
    
    commonPeak_result <- commonPeakResult(commonPeaks = commonPeak_output,
                                          return_common_peak_sites = TRUE)
    
    # generates a name for the list below
    common_names <- c(common_names, paste0("commonPeaks_", target_tf, "_", tf_y))
    # a list variable to store the results from this loop
    commonPeak_result_output <- c(commonPeak_result_output, commonPeak_result)
  }
}

# applying the names generated to the results lists 
names(exclusivePeak_result_output) <- exclusive_names
names(commonPeak_result_output) <- common_names

# example of the plotting process
plot(density(exclusivePeak_result_output$`exclusivePeaks_MM1_HSA_H1-hESC_USF1`$`MM1_HSA_H1-hESC_USF1_exclusive_peaks`$tag_fold_change))
points(density(commonPeak_result_output$`commonPeaks_MM1_HSA_H1-hESC_USF1_MM1_HSA_H1-hESC_USF2`$`MM1_HSA_H1-hESC_USF1_common_peaks`$tag_fold_change), type = "l", col = "red")
